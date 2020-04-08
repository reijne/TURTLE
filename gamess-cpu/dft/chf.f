c***********************************************************************
c
c  Chf.m
c
c  This file contains all the subroutines that relate to the CHF 
c  equations. That is right-hand-sides, left-hand-sides and 
c  perturbed Kohn-Sham matrices (no explicit derivative terms though).
c
c***********************************************************************
c
      subroutine ver_dft_chf(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/chf.m,v $
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
      subroutine den_pert_ao(rkstyp_sw,uamat,ubmat,
     &     nao,npert,
     &     bfn_val,bfng_val,
     &     drhodb,dgrhodb,t,
     &     npts,mxp,ider)
C     ******************************************************************
C     *Description:						       *
C     *Calculate the perturbed density and the perturbed grad density  *
c     *at a point                                                      *
C     ******************************************************************
      implicit none
c
c     Parameters
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
c     In variables
c
      integer npts, npert, mxp
      integer nao
      logical rkstyp_sw
      real*8 uamat(nao,nao,npert),ubmat(nao,nao,npert)
      real*8 bfn_val(mxp,*)
      real*8 bfng_val(mxp,nao,3)
      integer ider
c
C     Out variables
c
      real*8 drhodb(mxp,2,npert),dgrhodb(mxp,2,3,npert)
c
c     Work space
c
      real*8 t(npts)
c
c     Local variables
c
      integer ialpha, ibeta, ipert
      integer kbasi,lbasi

      integer ipt, k

      ialpha = 1
      ibeta  = 2

      call aclear_dp(drhodb,2*mxp*npert,0.0d0)
      if(ider.gt.0)call aclear_dp(dgrhodb,3*2*mxp*npert,0.0d0)
c     
c     Calculate drho/dRb  = sum(nu,sum(mu,Pb(mu,nu)chi(mu))chi(nu))
c     Calculate dgrho/dRb = sum(nu,sum(mu,Pb(mu,nu)chi(mu))dchi(nu))
c
      if (ider.le.0) then
         do ipert=1,npert
            do lbasi=1,nao
               do ipt=1,npts
                  t(ipt) = 0.0d0
               enddo
               do kbasi=1,nao
                  do ipt=1,npts
                     t(ipt)=t(ipt)
     &                     +uamat(kbasi,lbasi,ipert)*bfn_val(ipt,kbasi)
                  enddo
               enddo
               do ipt=1,npts
                  drhodb(ipt,1,ipert)=drhodb(ipt,1,ipert)
     &                               +t(ipt)*bfn_val(ipt,lbasi)
               enddo
            enddo
         enddo
      else
         do ipert=1,npert
            do lbasi=1,nao
               do ipt=1,npts
                  t(ipt) = 0.0d0
               enddo
               do kbasi=1,nao
                  do ipt=1,npts
                     t(ipt)=t(ipt)
     &                     +uamat(kbasi,lbasi,ipert)*bfn_val(ipt,kbasi)
                  enddo
               enddo
               do ipt=1,npts
                  drhodb(ipt,1,ipert)=drhodb(ipt,1,ipert)
     &                               +t(ipt)*bfn_val(ipt,lbasi)
                  dgrhodb(ipt,1,1,ipert)=dgrhodb(ipt,1,1,ipert)
     &                                  +t(ipt)*bfng_val(ipt,lbasi,1)
                  dgrhodb(ipt,1,2,ipert)=dgrhodb(ipt,1,2,ipert)
     &                                  +t(ipt)*bfng_val(ipt,lbasi,2)
                  dgrhodb(ipt,1,3,ipert)=dgrhodb(ipt,1,3,ipert)
     &                                  +t(ipt)*bfng_val(ipt,lbasi,3)
               enddo ! ipt
            enddo ! lbasi 
         enddo ! ipert
         do ipert=1,npert
            do k=1,3
               do lbasi=1,nao
                  do ipt=1,npts
                     t(ipt) = 0.0d0
                  enddo
                  do kbasi=1,nao
                     do ipt=1,npts
                        t(ipt)=t(ipt)
     &                     +uamat(kbasi,lbasi,ipert)*
     &                      bfng_val(ipt,kbasi,k)
                     enddo
                  enddo
                  do ipt=1,npts
                     dgrhodb(ipt,1,k,ipert)=dgrhodb(ipt,1,k,ipert)
     &                                     +t(ipt)*bfn_val(ipt,lbasi)
                  enddo ! ipt
               enddo ! lbasi 
            enddo ! k
         enddo ! ipert
      endif
c
      if (.not.rkstyp_sw) then
         if (ider.le.0) then
            do ipert=1,npert
               do lbasi=1,nao
                  do ipt=1,npts
                     t(ipt) = 0.0d0
                  enddo
                  do kbasi=1,nao
                     do ipt=1,npts
                        t(ipt)=t(ipt)
     &                        +ubmat(kbasi,lbasi,ipert)*
     &                         bfn_val(ipt,kbasi)
                     enddo
                  enddo
               enddo
               do ipt=1,npts
                  drhodb(ipt,2,ipert)=drhodb(ipt,2,ipert)
     &                               +t(ipt)*bfn_val(ipt,lbasi)
               enddo
            enddo
         else
            do ipert=1,npert
               do lbasi=1,nao
                  do ipt=1,npts
                     t(ipt) = 0.0d0
                  enddo
                  do kbasi=1,nao
                     do ipt=1,npts
                        t(ipt)=t(ipt)
     &                        +ubmat(kbasi,lbasi,ipert)*
     &                         bfn_val(ipt,kbasi)
                     enddo
                  enddo
               enddo
               do ipt=1,npts
                  drhodb(ipt,2,ipert)=drhodb(ipt,2,ipert)
     &                               +t(ipt)*bfn_val(ipt,lbasi)
                  dgrhodb(ipt,2,1,ipert)=dgrhodb(ipt,2,1,ipert)
     &                                  +t(ipt)*bfng_val(ipt,lbasi,1)
                  dgrhodb(ipt,2,2,ipert)=dgrhodb(ipt,2,2,ipert)
     &                                  +t(ipt)*bfng_val(ipt,lbasi,2)
                  dgrhodb(ipt,2,3,ipert)=dgrhodb(ipt,2,3,ipert)
     &                                 +t(ipt)*bfng_val(ipt,lbasi,3)
               enddo
            enddo
            do ipert=1,npert
               do k=1,3
                  do lbasi=1,nao
                     do ipt=1,npts
                        t(ipt) = 0.0d0
                     enddo
                     do kbasi=1,nao
                        do ipt=1,npts
                           t(ipt)=t(ipt)
     &                        +ubmat(kbasi,lbasi,ipert)*
     &                         bfng_val(ipt,kbasi,k)
                        enddo
                     enddo
                     do ipt=1,npts
                        dgrhodb(ipt,2,k,ipert)=dgrhodb(ipt,2,k,ipert)
     &                                        +t(ipt)*bfn_val(ipt,lbasi)
                     enddo ! ipt
                  enddo ! lbasi 
               enddo ! k
            enddo ! ipert
         endif
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine den_pert_ao_scr(rkstyp_sw,uamat,ubmat,
     &     active_bfn_list,n_active_bfn, 
     &     active_chf_pert,n_active_chf_prt,
     &     nao,npert,
     &     bfn_val,bfng_val,
     &     drhodb,dgrhodb,t,
     &     npts,mxp,ider,dentol)
C     ******************************************************************
C     *Description:						       *
C     *Calculate the perturbed density and the perturbed grad density  *
c     *at a point                                                      *
C     ******************************************************************
      implicit none

C     *Parameters	
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

C In variables	

c
      integer iky, ikyp, ilifq, mapie, ilifm, i4096
      common/mapper/iky(maxorb),ikyp(maxorb),ilifq(maxorb),
     +            mapie(maxorb),ilifm(maxorb),i4096(maxorb)
c

      integer npts, npert, mxp
      integer nao
      integer n_active_bfn
      integer active_bfn_list(n_active_bfn)
      integer n_active_chf_prt
      integer active_chf_pert(n_active_chf_prt)
      logical rkstyp_sw
      real*8 uamat(nao,nao,npert),ubmat(nao,nao,npert)
      real*8 bfn_val(mxp,*)
      real*8 bfng_val(mxp,nao,3)
      real*8 dentol
      
      integer ider
c
c     Out variables
c
      real*8 drhodb(mxp,2,npert),dgrhodb(mxp,2,3,npert)
c
c     Work space
c
      real*8 t(npts)
c
c     Local variables
c
      integer ialpha, ibeta, ipert, iprt
      integer kbasi,lbasi,k,kk,ll

      integer ipt

      ialpha = 1
      ibeta  = 2

      call aclear_dp(drhodb,2*mxp*npert,0.0d0)
      if(ider.gt.0)call aclear_dp(dgrhodb,3*2*mxp*npert,0.0d0)
C     
C     Calculate drho/dRb  = sum(nu,sum(mu,Pb(mu,nu)chi(mu))chi(nu))
C     Calculate dgrho/dRb = sum(nu,sum(mu,Pb(mu,nu)chi(mu))dchi(nu))
c
      if (ider.le.0) then
         do iprt=1,n_active_chf_prt
            ipert=active_chf_pert(iprt)
c           do l=1,3
c              iprt=3*(iatm-1)+l
c              ipert=lhs_atm_list(iatm)+l
               do lbasi=1,n_active_bfn
                  ll = active_bfn_list(lbasi)
                  do ipt=1,npts
                     t(ipt) = 0.0d0
                  enddo
                  do kbasi=1,n_active_bfn
                     kk = active_bfn_list(kbasi)
                     if (abs(uamat(kk,ll,ipert)).gt.dentol) then
                        do ipt=1,npts
                           t(ipt)=t(ipt)
     &                        +uamat(kk,ll,ipert)*bfn_val(ipt,kbasi)
                        enddo
                     endif
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,1,iprt)=drhodb(ipt,1,iprt)
     &                                 +t(ipt)*bfn_val(ipt,lbasi)
                  enddo
               enddo
c           enddo
         enddo
      else
         do iprt=1,n_active_chf_prt
            ipert=active_chf_pert(iprt)
c        do iatm=1,n_lhs_atm
c           do l=1,3
c              iprt=3*(iatm-1)+l
c              ipert=lhs_atm_list(iatm)+l
               do lbasi=1,n_active_bfn
                  ll = active_bfn_list(lbasi)
                  do ipt=1,npts
                     t(ipt) = 0.0d0
                  enddo
                  do kbasi=1,n_active_bfn
                     kk = active_bfn_list(kbasi)
                     if (abs(uamat(kk,ll,ipert)).gt.dentol) then
                        do ipt=1,npts
                           t(ipt)=t(ipt)
     &                     +uamat(kk,ll,ipert)*bfn_val(ipt,kbasi)
                        enddo
                     endif
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,1,iprt)=drhodb(ipt,1,iprt)
     &                                 +t(ipt)*bfn_val(ipt,lbasi)
                     dgrhodb(ipt,1,1,iprt)=dgrhodb(ipt,1,1,iprt)
     &                                    +t(ipt)*bfng_val(ipt,lbasi,1)
                     dgrhodb(ipt,1,2,iprt)=dgrhodb(ipt,1,2,iprt)
     &                                    +t(ipt)*bfng_val(ipt,lbasi,2)
                     dgrhodb(ipt,1,3,iprt)=dgrhodb(ipt,1,3,iprt)
     &                                    +t(ipt)*bfng_val(ipt,lbasi,3)
                 enddo ! ipt
              enddo ! lbasi 
c          enddo ! k
         enddo ! iatm
         do iprt=1,n_active_chf_prt
            ipert=active_chf_pert(iprt)
c        do iatm=1,n_lhs_atm
c           do l=1,3
c              iprt=3*(iatm-1)+l
c              ipert=lhs_atm_list(iatm)+l
               do k=1,3
                  do lbasi=1,n_active_bfn
                     ll = active_bfn_list(lbasi)
                     do ipt=1,npts
                        t(ipt) = 0.0d0
                     enddo
                     do kbasi=1,n_active_bfn
                        kk = active_bfn_list(kbasi)
                        if (abs(uamat(kk,ll,ipert)).gt.dentol) then
                           do ipt=1,npts
                              t(ipt)=t(ipt)
     &                        +uamat(kk,ll,ipert)*
     &                         bfng_val(ipt,kbasi,k)
                           enddo
                        endif
                     enddo
                     do ipt=1,npts
                        dgrhodb(ipt,1,k,iprt)=dgrhodb(ipt,1,k,iprt)
     &                                       +t(ipt)*bfn_val(ipt,lbasi)
                     enddo ! ipt
                  enddo ! lbasi 
               enddo ! k
c           enddo ! l
         enddo ! iatm
      endif
c
      if (.not.rkstyp_sw) then
         if (ider.le.0) then
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do l=1,3
c                 iprt=3*(iatm-1)+l
c                 ipert=lhs_atm_list(iatm)+l
                  do lbasi=1,n_active_bfn
                     ll = active_bfn_list(lbasi)
                     do ipt=1,npts
                        t(ipt) = 0.0d0
                     enddo
                     do kbasi=1,n_active_bfn
                        kk = active_bfn_list(kbasi)
                        if (abs(ubmat(kk,ll,ipert)).gt.dentol) then
                           do ipt=1,npts
                              t(ipt)=t(ipt)
     &                              +ubmat(kk,ll,ipert)*
     &                               bfn_val(ipt,kbasi)
                           enddo
                        endif
                     enddo
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,2,iprt)=drhodb(ipt,2,iprt)
     &                                 +t(ipt)*bfn_val(ipt,lbasi)
                  enddo
c              enddo
            enddo
         else
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do l=1,3
c                 iprt=3*(iatm-1)+l
c                 ipert=lhs_atm_list(iatm)+l
                  do lbasi=1,n_active_bfn
                     ll = active_bfn_list(lbasi)
                     do ipt=1,npts
                        t(ipt) = 0.0d0
                     enddo
                     do kbasi=1,n_active_bfn
                        kk = active_bfn_list(kbasi)
c                       dd = bdens(iky(max(ll,kk))+min(ll,kk))
                        if (abs(ubmat(kk,ll,ipert)).gt.dentol) then
                           do ipt=1,npts
                              t(ipt)=t(ipt)
     &                              +ubmat(kk,ll,ipert)*
     &                               bfn_val(ipt,kbasi)
                           enddo
                        endif
                     enddo
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,2,iprt)=drhodb(ipt,2,iprt)
     &                                 +t(ipt)*bfn_val(ipt,lbasi)
                     dgrhodb(ipt,2,1,iprt)=dgrhodb(ipt,2,1,iprt)
     &                                    +t(ipt)*bfng_val(ipt,lbasi,1)
                     dgrhodb(ipt,2,2,iprt)=dgrhodb(ipt,2,2,iprt)
     &                                    +t(ipt)*bfng_val(ipt,lbasi,2)
                     dgrhodb(ipt,2,3,iprt)=dgrhodb(ipt,2,3,iprt)
     &                                    +t(ipt)*bfng_val(ipt,lbasi,3)
                  enddo
c              enddo
            enddo
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do l=1,3
c                 iprt=3*(iatm-1)+l
c                 ipert=lhs_atm_list(iatm)+l
                  do k=1,3
                     do lbasi=1,n_active_bfn
                        ll = active_bfn_list(lbasi)
                        do ipt=1,npts
                           t(ipt) = 0.0d0
                        enddo
                        do kbasi=1,n_active_bfn
                           kk = active_bfn_list(kbasi)
c                          dd = bdens(iky(max(ll,kk))+min(ll,kk))
                           if (abs(ubmat(kk,ll,ipert)).gt.dentol) then
                              do ipt=1,npts
                                 t(ipt)=t(ipt)
     &                                 +ubmat(kk,ll,ipert)*
     &                                  bfng_val(ipt,kbasi,k)
                              enddo
                           endif
                        enddo
                        do ipt=1,npts
                           dgrhodb(ipt,2,k,iprt)=dgrhodb(ipt,2,k,iprt)
     &                                          +t(ipt)*
     &                                           bfn_val(ipt,lbasi)
                        enddo ! ipt
                     enddo ! lbasi 
                  enddo ! k
c              enddo ! l
            enddo ! iatm
         endif
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine den_pert_dksm_mo(rkstyp_sw,ider,npts,mxp,npert,nvec,
     &                         naocc,nbocc,uamat,ubmat,
     &                         amo_val,amo_grad,
     &                         bmo_val,bmo_grad,
     &                         drhodb,dgrhodb,t,dentol)
      implicit none
c
c     Compute the wavefunction part of the derivatives of rho and grho 
c     with respect to nuclear coordinates. For the left-hand-sides
c     these dependent on the occupied-virtual block of the U matrix.
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
c     Inputs
c
      logical rkstyp_sw
      integer ider
      integer npts
      integer mxp
      integer npert
      integer nvec
      integer naocc
      integer nbocc
      real*8 uamat(nvec*(nvec+1)/2,npert)
      real*8 ubmat(nvec*(nvec+1)/2,npert)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 bmo_val(npts,nvec)
      real*8 bmo_grad(npts,nvec,3)
      real*8 dentol
c
c     Outputs
c
      real*8 drhodb(mxp,2,npert)
      real*8 dgrhodb(mxp,2,3,npert)
c
c     Workspace
c
      real*8 t(npts)
c
c     Local variables
c
      integer i, j, k, iki, maxj
      integer iprt
      integer ipt
      real*8 uval
c
c     Code
c
      call aclear_dp(drhodb,mxp*2*npert,0.0d0)
      if (ider.ge.1) then
         call aclear_dp(dgrhodb,mxp*2*3*npert,0.0d0)
      endif
c
      if (ider.le.0) then
         do iprt=1,npert
            do i=1,nvec
               iki = iky(i)
               call aclear_dp(t,npts,0.0d0)
               maxj=min(i,naocc)
               do j=1,maxj
                  uval=uamat(iki+j,iprt)
                  if (abs(uval).gt.dentol) then
                     do ipt=1,npts
                        t(ipt)=t(ipt)+amo_val(ipt,j)*uamat(iki+j,iprt)
                     enddo
                  endif
               enddo
               do ipt=1,npts
                  drhodb(ipt,1,iprt)=drhodb(ipt,1,iprt)
     &                              +2*amo_val(ipt,i)*t(ipt)
               enddo
            enddo
         enddo
      else if (ider.le.1) then
         do iprt=1,npert
            do i=1,naocc
               iki = iky(i)
               do ipt=1,npts
                  t(ipt)=2.0d0*amo_val(ipt,i)*uamat(iki+i,iprt)
               enddo
               do j=1,i-1
                  do ipt=1,npts
                     t(ipt)=t(ipt)+amo_val(ipt,j)*uamat(iki+j,iprt)
                  enddo
               enddo
               do j=i+1,naocc
                  iki = iky(j)
                  do ipt=1,npts
                     t(ipt)=t(ipt)+amo_val(ipt,j)*uamat(iki+i,iprt)
                  enddo
               enddo
               do ipt=1,npts
                  drhodb(ipt,1,iprt)=drhodb(ipt,1,iprt)
     &                              +amo_val(ipt,i)*t(ipt)
                  t(ipt)=2.0d0*t(ipt)
               enddo
               do k=1,3
                  do ipt=1,npts
                     dgrhodb(ipt,1,k,iprt)=dgrhodb(ipt,1,k,iprt)
     &                                    +amo_grad(ipt,i,k)*
     &                                     t(ipt)
                  enddo
               enddo
            enddo
            do i=naocc+1,nvec
               iki = iky(i)
               do ipt=1,npts
                  t(ipt)=0.0d0
               enddo
               do j=1,naocc
                  do ipt=1,npts
                     t(ipt)=t(ipt)+amo_val(ipt,j)*uamat(iki+j,iprt)
                  enddo
               enddo
               do ipt=1,npts
                  t(ipt)=2.0d0*t(ipt)
               enddo
               do ipt=1,npts
                  drhodb(ipt,1,iprt)=drhodb(ipt,1,iprt)
     &                              +amo_val(ipt,i)*t(ipt)
               enddo
               do k=1,3
                  do ipt=1,npts
                     dgrhodb(ipt,1,k,iprt)=dgrhodb(ipt,1,k,iprt)
     &                                    +amo_grad(ipt,i,k)*
     &                                     t(ipt)
                  enddo
               enddo
            enddo
            do j=1,naocc
               do ipt=1,npts
                  t(ipt)=0.0d0
               enddo
               do i=naocc+1,nvec
                  iki=iky(i)
                  do ipt=1,npts
                     t(ipt)=t(ipt)+amo_val(ipt,i)*
     &                             uamat(iki+j,iprt)
                  enddo
               enddo
               do ipt=1,npts
                  t(ipt)=2*t(ipt)
               enddo
               do k=1,3
                  do ipt=1,npts
                     dgrhodb(ipt,1,k,iprt)=dgrhodb(ipt,1,k,iprt)
     &                                    +amo_grad(ipt,j,k)*
     &                                     t(ipt)
                  enddo
               enddo
            enddo
         enddo
      endif
      if (.not.rkstyp_sw) then
         if (ider.le.0) then
            do iprt=1,npert
               do i=1,nvec
                  iki=iky(i)
                  do ipt=1,npts
                     t(ipt)=0.0d0
                  enddo
                  maxj=min(i,nbocc)
                  do j=1,maxj
                     do ipt=1,npts
                        t(ipt)=t(ipt)+bmo_val(ipt,j)*ubmat(iki+j,iprt)
                     enddo
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,2,iprt)=drhodb(ipt,2,iprt)
     &                                 +2*bmo_val(ipt,i)*t(ipt)
                  enddo
               enddo
            enddo
         else if (ider.le.1) then
            do iprt=1,npert
               do i=1,nbocc
                  do ipt=1,npts
                     t(ipt)=2.0d0*bmo_val(ipt,i)*ubmat(iki+i,iprt)
                  enddo
                  do j=1,i-1
                     do ipt=1,npts
                        t(ipt)=t(ipt)+bmo_val(ipt,j)*ubmat(iki+j,iprt)
                     enddo
                  enddo
                  do j=i+1,nbocc
                     iki=iky(j)
                     do ipt=1,npts
                        t(ipt)=t(ipt)+bmo_val(ipt,j)*ubmat(iki+i,iprt)
                     enddo
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,2,iprt)=drhodb(ipt,2,iprt)
     &                                 +bmo_val(ipt,i)*t(ipt)
                     t(ipt)=2.0d0*t(ipt)
                  enddo
                  do k=1,3
                     do ipt=1,npts
                        dgrhodb(ipt,2,k,iprt)=dgrhodb(ipt,2,k,iprt)
     &                                       +bmo_grad(ipt,i,k)*
     &                                        t(ipt)
                     enddo
                  enddo
               enddo
               do i=nbocc+1,nvec
                  iki = iky(i)
                  do ipt=1,npts
                     t(ipt)=0.0d0
                  enddo
                  do j=1,nbocc
                     do ipt=1,npts
                        t(ipt)=t(ipt)+bmo_val(ipt,j)*ubmat(iki+j,iprt)
                     enddo
                  enddo
                  do ipt=1,npts
                     t(ipt)=2.0d0*t(ipt)
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,2,iprt)=drhodb(ipt,2,iprt)
     &                                 +bmo_val(ipt,i)*t(ipt)
                  enddo
                  do k=1,3
                     do ipt=1,npts
                        dgrhodb(ipt,2,k,iprt)=dgrhodb(ipt,2,k,iprt)
     &                                       +bmo_grad(ipt,i,k)*
     &                                        t(ipt)
                     enddo
                  enddo
               enddo
               do j=1,nbocc
                  do ipt=1,npts
                     t(ipt)=0.0d0
                  enddo
                  do i=nbocc+1,nvec
                     iki=iky(i)
                     do ipt=1,npts
                        t(ipt)=t(ipt)+bmo_val(ipt,i)*
     &                                ubmat(iki+j,iprt)
                     enddo
                  enddo
                  do ipt=1,npts
                     t(ipt)=2*t(ipt)
                  enddo
                  do k=1,3
                     do ipt=1,npts
                        dgrhodb(ipt,2,k,iprt)=dgrhodb(ipt,2,k,iprt)
     &                                       +bmo_grad(ipt,j,k)*
     &                                        t(ipt)
                     enddo
                  enddo
               enddo
            enddo
         endif
      endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine den_pert_dksm_mo_scr(rkstyp_sw,ider,npts,mxp,npert,
     &                         nvec,naocc,nbocc,
     &                         active_chf_pert,n_active_chf_prt,
     &                         uamat,ubmat,
     &                         amo_val,amo_grad,
     &                         bmo_val,bmo_grad,
     &                         drhodb,dgrhodb,t,dentol)
      implicit none
c
c     Compute the wavefunction part of the derivatives of rho and grho 
c     with respect to nuclear coordinates. For the left-hand-sides
c     these dependent on the occupied-virtual block of the U matrix.
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
c     Inputs
c
      logical rkstyp_sw
      integer ider
      integer npts
      integer mxp
      integer npert
      integer nvec
      integer naocc
      integer nbocc
      integer n_active_chf_prt
      integer active_chf_pert(n_active_chf_prt)
      real*8 uamat(nvec*(nvec+1)/2,npert)
      real*8 ubmat(nvec*(nvec+1)/2,npert)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 bmo_val(npts,nvec)
      real*8 bmo_grad(npts,nvec,3)
      real*8 dentol
c
c     Outputs
c
      real*8 drhodb(mxp,2,npert)
      real*8 dgrhodb(mxp,2,3,npert)
c
c     Workspace
c
      real*8 t(npts)
c
c     Local variables
c
      integer i, j, k, iki, maxj
      integer ipert
      integer iprt 
c     integer iatm, ic
      integer ipt
      real*8 uval
c
c     Code
c
      call aclear_dp(drhodb,mxp*2*npert,0.0d0)
      if (ider.ge.1) then
         call aclear_dp(dgrhodb,mxp*2*3*npert,0.0d0)
      endif
c
      if (ider.le.0) then
         do iprt=1,n_active_chf_prt
            ipert=active_chf_pert(iprt)
c        do iatm=1,n_lhs_atm
c           do ic=1,3
c              ipert=lhs_atm_list(iatm)+ic
c              iprt=3*(iatm-1)+ic
               do i=1,nvec
                  iki = iky(i)
                  call aclear_dp(t,npts,0.0d0)
                  maxj=min(i,naocc)
                  do j=1,maxj
                     uval=uamat(iki+j,ipert)
                     if (abs(uval).gt.dentol) then
                        do ipt=1,npts
                           t(ipt)=t(ipt)+amo_val(ipt,j)*uval
                        enddo
                     endif
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,1,iprt)=drhodb(ipt,1,iprt)
     &                                 +2*amo_val(ipt,i)*t(ipt)
                  enddo
               enddo
c           enddo
         enddo
      else if (ider.le.1) then
         do iprt=1,n_active_chf_prt
            ipert=active_chf_pert(iprt)
c        do iatm=1,n_lhs_atm
c           do ic=1,3
c              ipert=lhs_atm_list(iatm)+ic
c              iprt=3*(iatm-1)+ic
               do i=1,naocc
                  iki = iky(i)
                  uval=uamat(iki+i,ipert)
                  if (abs(uval).gt.dentol) then
                     do ipt=1,npts
                        t(ipt)=2.0d0*amo_val(ipt,i)*uval
                     enddo
                  else
                     do ipt=1,npts
                        t(ipt)=0.0d0
                     enddo
                  endif
                  do j=1,i-1
                     uval=uamat(iki+j,ipert)
                     if (abs(uval).gt.dentol) then
                        do ipt=1,npts
                           t(ipt)=t(ipt)+amo_val(ipt,j)*uval
                        enddo
                     endif
                  enddo
                  do j=i+1,naocc
                     iki = iky(j)
                     uval=uamat(iki+i,ipert)
                     if (abs(uval).gt.dentol) then
                        do ipt=1,npts
                           t(ipt)=t(ipt)+amo_val(ipt,j)*uval
                        enddo
                     endif
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,1,iprt)=drhodb(ipt,1,iprt)
     &                                 +amo_val(ipt,i)*t(ipt)
                     t(ipt)=2.0d0*t(ipt)
                  enddo
                  do k=1,3
                     do ipt=1,npts
                        dgrhodb(ipt,1,k,iprt)=dgrhodb(ipt,1,k,iprt)
     &                                       +amo_grad(ipt,i,k)*
     &                                        t(ipt)
                     enddo
                  enddo
               enddo
               do i=naocc+1,nvec
                  iki = iky(i)
                  do ipt=1,npts
                     t(ipt)=0.0d0
                  enddo
                  do j=1,naocc
                     uval=uamat(iki+j,ipert)
                     if (abs(uval).gt.dentol) then
                        do ipt=1,npts
                           t(ipt)=t(ipt)+amo_val(ipt,j)*uval
                        enddo
                     endif
                  enddo
                  do ipt=1,npts
                     t(ipt)=2.0d0*t(ipt)
                  enddo
                  do ipt=1,npts
                     drhodb(ipt,1,iprt)=drhodb(ipt,1,iprt)
     &                                 +amo_val(ipt,i)*t(ipt)
                  enddo
                  do k=1,3
                     do ipt=1,npts
                        dgrhodb(ipt,1,k,iprt)=dgrhodb(ipt,1,k,iprt)
     &                                       +amo_grad(ipt,i,k)*
     &                                        t(ipt)
                     enddo
                  enddo
               enddo
               do j=1,naocc
                  do ipt=1,npts
                     t(ipt)=0.0d0
                  enddo
                  do i=naocc+1,nvec
                     iki=iky(i)
                     uval=uamat(iki+j,ipert)
                     if (abs(uval).gt.dentol) then
                        do ipt=1,npts
                           t(ipt)=t(ipt)+amo_val(ipt,i)*uval
                        enddo
                     endif
                  enddo
                  do ipt=1,npts
                     t(ipt)=2*t(ipt)
                  enddo
                  do k=1,3
                     do ipt=1,npts
                        dgrhodb(ipt,1,k,iprt)=dgrhodb(ipt,1,k,iprt)
     &                                       +amo_grad(ipt,j,k)*
     &                                        t(ipt)
                     enddo
                  enddo
               enddo
c           enddo
         enddo
      endif
      if (.not.rkstyp_sw) then
         if (ider.le.0) then
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do ic=1,3
c                 ipert=lhs_atm_list(iatm)+ic
c                 iprt=3*(iatm-1)+ic
                  do i=1,nvec
                     iki=iky(i)
                     do ipt=1,npts
                        t(ipt)=0.0d0
                     enddo
                     maxj=min(i,nbocc)
                     do j=1,maxj
                        uval=ubmat(iki+j,ipert)
                        if (abs(uval).gt.dentol) then
                           do ipt=1,npts
                              t(ipt)=t(ipt)+bmo_val(ipt,j)*uval
                           enddo
                        endif
                     enddo
                     do ipt=1,npts
                        drhodb(ipt,2,iprt)=drhodb(ipt,2,iprt)
     &                                    +2*bmo_val(ipt,i)*t(ipt)
                     enddo
                  enddo
c              enddo
            enddo
         else if (ider.le.1) then
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do ic=1,3
c                 ipert=lhs_atm_list(iatm)+ic
c                 iprt=3*(iatm-1)+ic
                  do i=1,nbocc
                     uval=ubmat(iki+i,ipert)
                     if (abs(uval).gt.dentol) then
                        do ipt=1,npts
                           t(ipt)=2.0d0*bmo_val(ipt,i)*uval
                        enddo
                     else
                        do ipt=1,npts
                           t(ipt)=0.0d0
                        enddo
                     endif
                     do j=1,i-1
                        uval=ubmat(iki+j,ipert)
                        if (abs(uval).gt.dentol) then
                           do ipt=1,npts
                              t(ipt)=t(ipt)+bmo_val(ipt,j)*uval
                           enddo
                        endif
                     enddo
                     do j=i+1,nbocc
                        iki=iky(j)
                        uval=ubmat(iki+i,ipert)
                        if (abs(uval).gt.dentol) then
                           do ipt=1,npts
                              t(ipt)=t(ipt)+bmo_val(ipt,j)*uval
                           enddo
                        endif
                     enddo
                     do ipt=1,npts
                        drhodb(ipt,2,iprt)=drhodb(ipt,2,iprt)
     &                                    +bmo_val(ipt,i)*t(ipt)
                        t(ipt)=2.0d0*t(ipt)
                     enddo
                     do k=1,3
                        do ipt=1,npts
                           dgrhodb(ipt,2,k,iprt)=dgrhodb(ipt,2,k,iprt)
     &                                          +bmo_grad(ipt,i,k)*
     &                                           t(ipt)
                        enddo
                     enddo
                  enddo
                  do i=nbocc+1,nvec
                     iki = iky(i)
                     do ipt=1,npts
                        t(ipt)=0.0d0
                     enddo
                     do j=1,nbocc
                        uval=ubmat(iki+j,ipert)
                        if (abs(uval).gt.dentol) then
                           do ipt=1,npts
                              t(ipt)=t(ipt)+bmo_val(ipt,j)*uval
                           enddo
                        endif
                     enddo
                     do ipt=1,npts
                        t(ipt)=2.0d0*t(ipt)
                     enddo
                     do ipt=1,npts
                        drhodb(ipt,2,iprt)=drhodb(ipt,2,iprt)
     &                                    +bmo_val(ipt,i)*t(ipt)
                     enddo
                     do k=1,3
                        do ipt=1,npts
                           dgrhodb(ipt,2,k,iprt)=dgrhodb(ipt,2,k,iprt)
     &                                          +bmo_grad(ipt,i,k)*
     &                                           t(ipt)
                        enddo
                     enddo
                  enddo
                  do j=1,nbocc
                     do ipt=1,npts
                        t(ipt)=0.0d0
                     enddo
                     do i=nbocc+1,nvec
                        iki=iky(i)
                        uval=ubmat(iki+j,ipert)
                        if (abs(uval).gt.dentol) then
                           do ipt=1,npts
                              t(ipt)=t(ipt)+bmo_val(ipt,i)*uval
                           enddo
                        endif
                     enddo
                     do ipt=1,npts
                        t(ipt)=2*t(ipt)
                     enddo
                     do k=1,3
                        do ipt=1,npts
                           dgrhodb(ipt,2,k,iprt)=dgrhodb(ipt,2,k,iprt)
     &                                          +bmo_grad(ipt,j,k)*
     &                                           t(ipt)
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
         endif
      endif
      end
c
c-----------------------------------------------------------------------
c
      subroutine gamma_pert(rkstyp_sw,npts,npert,grho,dgrhodb,dgammadb,
     &                      mxp)
      implicit none
c
c     Computes the perturbed gamma's on the grid. I.e. it computes the
c     various inproducts between (grad rho) and the perturbed 
c     (grad rho)'s.
c
c     In variables:
c
      logical rkstyp_sw ! .true. if closed shell
      integer mxp
      integer npts
      integer npert
      real*8 grho(mxp,2,3)
      real*8 dgrhodb(mxp,2,3,npert)
c
c     Out variables:
c
      real*8 dgammadb(mxp,3,npert)
c
c     Parameters:
c
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
c     Local variables:
c
      integer ipert, ipt
c
c     Code:
c
      if (rkstyp_sw) then
         do ipert = 1, npert
            do ipt = 1, npts
c              dgammadb(ipt,igaa,ipert) = 2.0d0*
               dgammadb(ipt,igaa,ipert) = 
     &            (grho(ipt,1,1)*dgrhodb(ipt,1,1,ipert)
     &            +grho(ipt,1,2)*dgrhodb(ipt,1,2,ipert)
     &            +grho(ipt,1,3)*dgrhodb(ipt,1,3,ipert))
            enddo
         enddo
      else
         do ipert = 1, npert
            do ipt = 1, npts
               dgammadb(ipt,igaa,ipert) = 2.0d0*
     &            (grho(ipt,1,1)*dgrhodb(ipt,1,1,ipert)
     &            +grho(ipt,1,2)*dgrhodb(ipt,1,2,ipert)
     &            +grho(ipt,1,3)*dgrhodb(ipt,1,3,ipert))
               dgammadb(ipt,igbb,ipert) = 2.0d0*
     &            (grho(ipt,2,1)*dgrhodb(ipt,2,1,ipert)
     &            +grho(ipt,2,2)*dgrhodb(ipt,2,2,ipert)
     &            +grho(ipt,2,3)*dgrhodb(ipt,2,3,ipert))
               dgammadb(ipt,igab,ipert) = 
     &            (grho(ipt,1,1)*dgrhodb(ipt,2,1,ipert)
     &            +grho(ipt,1,2)*dgrhodb(ipt,2,2,ipert)
     &            +grho(ipt,1,3)*dgrhodb(ipt,2,3,ipert)
     &            +grho(ipt,2,1)*dgrhodb(ipt,1,1,ipert)
     &            +grho(ipt,2,2)*dgrhodb(ipt,1,2,ipert)
     &            +grho(ipt,2,3)*dgrhodb(ipt,1,3,ipert))
            enddo
         enddo
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine lhs_dft_ao(rkstyp_sw,gradcorr_sw,nao,npts,npert,wt,
     &     grho,drhodb,dgrhodb,dgammadb,bfn_val,bfng_val,xc_dvpt,
     &     xc_hpt,xc_dhpt,ga,gb,mxp,wrrp,wrrpd,wrrpdn,drmn,dgrn,dvn,
     &     ddr,dhnd)
      implicit none
c
c     Computes the DFT terms that contribute to the Left-Hand-Side
c     of the CPHF equations.
c
c     For performance analysis we introduce the following variables:
c     N  = the number of atoms
c     Np = a*N = the number of grid points (a = +- thousands)
c     Nb = b*N = the number of basis functions (b = 1 to 2 dozen)
c
c     Parameters
c
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
c     In variables
c
      logical rkstyp_sw ! .true. if closed shell
      logical gradcorr_sw
      integer nao, npts, npert, mxp
      real*8 wt(mxp)
      real*8 grho(mxp,2,3)
      real*8 drhodb(mxp,2,npert)
      real*8 dgrhodb(mxp,2,3,npert)
      real*8 dgammadb(mxp,3,npert)
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
c     Out variables
c
      real*8 ga(nao,nao,npert)
      real*8 gb(nao,nao,npert)
c
c     Local variables
c
      integer ipert, ipt, imu, inu
      real*8 rhoamn, rhobmn, rhoamnb, rhobmnb

c     real*8 wrr, wrg, wrg2, wgg, wg
c     real*8 wrrdr, wrgdr, wrgdg, wggdg
c     real*8 grgb(nao), dgrgb(nao)
c     real*8 wrrdrn(nao), wrgdrn(nao), wrgdgn(nao), wggdgn(nao)
c     real*8 wrrn(nao), rbn(nao)
      real*8 wrrp(npts), wrrpd(npts,npert), wrrpdn(npts,nao,npert)
      real*8 drmn(npts,nao,nao), dgrn(npts,nao,npert), dvn(npts,nao)
      real*8 ddr(npts,npert), dhnd(npts,nao,npert)
c
c     Code
c
      if (gradcorr_sw) then
         if (rkstyp_sw) then
c
c           The following code section costs (FLOPs)
c           3*N*Nb*Nb*Na*(13+13+22)
c           = 144*a*b*b*N**4
c
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c
c                       rhoamn = 
c    &                    (grho(ipt,1,1)*bfng_val(ipt,inu,1)+
c    &                     grho(ipt,1,2)*bfng_val(ipt,inu,2)+
c    &                     grho(ipt,1,3)*bfng_val(ipt,inu,3))*
c    &                    bfn_val(ipt,imu)+
c    &                    (grho(ipt,1,1)*bfng_val(ipt,imu,1)+
c    &                     grho(ipt,1,2)*bfng_val(ipt,imu,2)+
c    &                     grho(ipt,1,3)*bfng_val(ipt,imu,3))*
c    &                    bfn_val(ipt,inu)
c
c                       rhoamnb = 
c    &                    (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,inu,1)+
c    &                     dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,inu,2)+
c    &                     dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,inu,3))*
c    &                    bfn_val(ipt,imu)+
c    &                    (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,imu,1)+
c    &                     dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,imu,2)+
c    &                     dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,imu,3))*
c    &                    bfn_val(ipt,inu)
c
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                  + wt(ipt)*xc_hpt(ipt,irara)*
c    &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
c    &                    drhodb(ipt,1,ipert)
c    &                  + wt(ipt)*xc_dhpt(ipt,iragaa)*
c    &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
c    &                    dgammadb(ipt,1,ipert)
c    &                  + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*rhoamn*
c    &                      drhodb(ipt,1,ipert)
c    &                  +   wt(ipt)*xc_dhpt(ipt,igaagaa)*rhoamn*
c    &                      dgammadb(ipt,1,ipert)
c    &                  +   wt(ipt)*xc_dvpt(ipt,igaa)*rhoamnb
c
c                    enddo
c                 enddo
c              enddo
c           enddo
c
c           The following section of code costs (FLOPs)
c           Nb*Np*3 + Nb*Nb*Np*13 + 3*N*Nb*Np*(5+6) + 3*N*Np*5
c           + 3*N*Nb*Nb*Np*8
c           = 3*a*b*N**2 + 15*a*N**2
c           + 13*a*b*b*N**3 + 33*a*b*N**3
c           + 24*a*b*b*N**4
c           = saving of max 6 (realized on methane 2)
c
            do inu = 1, nao
               do ipt = 1, npts
                  dvn(ipt,inu) = wt(ipt)*xc_dvpt(ipt,igaa)*
     &                           bfn_val(ipt,inu)
               enddo
            enddo
            do inu = 1, nao
               do imu = 1, nao
                  do ipt = 1, npts
                     drmn(ipt,imu,inu) = 
     &                   (grho(ipt,1,1)*bfng_val(ipt,inu,1)+
     &                    grho(ipt,1,2)*bfng_val(ipt,inu,2)+
     &                    grho(ipt,1,3)*bfng_val(ipt,inu,3))*
     &                   bfn_val(ipt,imu)+
     &                   (grho(ipt,1,1)*bfng_val(ipt,imu,1)+
     &                    grho(ipt,1,2)*bfng_val(ipt,imu,2)+
     &                    grho(ipt,1,3)*bfng_val(ipt,imu,3))*
     &                   bfn_val(ipt,inu)
                  enddo
               enddo
            enddo
            do ipert = 1, npert
               do inu = 1, nao
                  do ipt = 1, npts
                     dgrn(ipt,inu,ipert) = 
     &                   (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,inu,1)+
     &                    dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,inu,2)+
     &                    dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,inu,3))
                     dhnd(ipt,inu,ipert) = wt(ipt)*(
     &                    xc_hpt(ipt,irara)*bfn_val(ipt,inu)*
     &                    drhodb(ipt,1,ipert)
     &                  + xc_dhpt(ipt,iragaa)*bfn_val(ipt,inu)*
     &                    dgammadb(ipt,1,ipert))
                  enddo
               enddo
            enddo
            do ipert = 1, npert
               do ipt = 1, npts
                  ddr(ipt,ipert) = wt(ipt)*(
     &                2*xc_dhpt(ipt,iragaa)*drhodb(ipt,1,ipert)
     &              + xc_dhpt(ipt,igaagaa)*dgammadb(ipt,1,ipert))
               enddo
            enddo
            do ipert = 1, npert
               do inu = 1, nao
                  do imu = 1, nao
                     do ipt = 1, npts
                        ga(imu,inu,ipert) = ga(imu,inu,ipert)
     &                  + dhnd(ipt,inu,ipert)*bfn_val(ipt,imu)
     &                  + ddr(ipt,ipert)*drmn(ipt,imu,inu)
     &                  + dvn(ipt,imu)*dgrn(ipt,inu,ipert)
     &                  + dvn(ipt,inu)*dgrn(ipt,imu,ipert)
                     enddo
                  enddo
               enddo
            enddo
         else ! rkstyp_sw
            do ipert = 1, npert
               do inu = 1, nao
                  do imu = 1, nao
                     do ipt = 1, npts
c
                        rhoamn = 
     &                    (grho(ipt,1,1)*bfng_val(ipt,inu,1)+
     &                     grho(ipt,1,2)*bfng_val(ipt,inu,2)+
     &                     grho(ipt,1,3)*bfng_val(ipt,inu,3))*
     &                    bfn_val(ipt,imu)+
     &                    (grho(ipt,1,1)*bfng_val(ipt,imu,1)+
     &                     grho(ipt,1,2)*bfng_val(ipt,imu,2)+
     &                     grho(ipt,1,3)*bfng_val(ipt,imu,3))*
     &                    bfn_val(ipt,inu)
c
                        rhobmn = 
     &                    (grho(ipt,2,1)*bfng_val(ipt,inu,1)+
     &                     grho(ipt,2,2)*bfng_val(ipt,inu,2)+
     &                     grho(ipt,2,3)*bfng_val(ipt,inu,3))*
     &                    bfn_val(ipt,imu)+
     &                    (grho(ipt,2,1)*bfng_val(ipt,imu,1)+
     &                     grho(ipt,2,2)*bfng_val(ipt,imu,2)+
     &                     grho(ipt,2,3)*bfng_val(ipt,imu,3))*
     &                    bfn_val(ipt,inu)
c
                        rhoamnb = 
     &                    (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,inu,1)+
     &                     dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,inu,2)+
     &                     dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,inu,3))*
     &                    bfn_val(ipt,imu)+
     &                    (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,imu,1)+
     &                     dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,imu,2)+
     &                     dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,imu,3))*
     &                    bfn_val(ipt,inu)
c
                        rhobmnb = 
     &                    (dgrhodb(ipt,2,1,ipert)*bfng_val(ipt,inu,1)+
     &                     dgrhodb(ipt,2,2,ipert)*bfng_val(ipt,inu,2)+
     &                     dgrhodb(ipt,2,3,ipert)*bfng_val(ipt,inu,3))*
     &                    bfn_val(ipt,imu)+
     &                    (dgrhodb(ipt,2,1,ipert)*bfng_val(ipt,imu,1)+
     &                     dgrhodb(ipt,2,2,ipert)*bfng_val(ipt,imu,2)+
     &                     dgrhodb(ipt,2,3,ipert)*bfng_val(ipt,imu,3))*
     &                    bfn_val(ipt,inu)
c
                        ga(imu,inu,ipert) = ga(imu,inu,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irara)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    drhodb(ipt,1,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irarb)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    drhodb(ipt,2,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    dgammadb(ipt,1,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,iragab)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    dgammadb(ipt,2,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,iragbb)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    dgammadb(ipt,3,ipert)
c
     &                  + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*rhoamn*
     &                      drhodb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,iragab)*rhobmn*
     &                      drhodb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,irbgaa)*rhoamn*
     &                      drhodb(ipt,2,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,irbgab)*rhobmn*
     &                      drhodb(ipt,2,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igaagaa)*rhoamn*
     &                      dgammadb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igaagab)*rhobmn*
     &                      dgammadb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igaagab)*rhoamn*
     &                      dgammadb(ipt,2,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igabgab)*rhobmn*
     &                      dgammadb(ipt,2,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igaagbb)*rhoamn*
     &                      dgammadb(ipt,3,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igabgbb)*rhobmn*
     &                      dgammadb(ipt,3,ipert)
c
     &                  + 2*wt(ipt)*xc_dvpt(ipt,igaa)*rhoamnb
     &                  +   wt(ipt)*xc_dvpt(ipt,igab)*rhobmnb
c
                        gb(imu,inu,ipert) = gb(imu,inu,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irarb)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    drhodb(ipt,1,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    drhodb(ipt,2,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgaa)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    dgammadb(ipt,1,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgab)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    dgammadb(ipt,2,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgbb)*
     &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                    dgammadb(ipt,3,ipert)
c
     &                  + 2*wt(ipt)*xc_dhpt(ipt,iragbb)*rhobmn*
     &                      drhodb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,iragab)*rhoamn*
     &                      drhodb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,irbgbb)*rhobmn*
     &                      drhodb(ipt,2,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,irbgab)*rhoamn*
     &                      drhodb(ipt,2,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igaagbb)*rhobmn*
     &                      dgammadb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igaagab)*rhoamn*
     &                      dgammadb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igabgbb)*rhobmn*
     &                      dgammadb(ipt,2,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igabgab)*rhoamn*
     &                      dgammadb(ipt,2,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*rhobmn*
     &                      dgammadb(ipt,3,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igabgbb)*rhoamn*
     &                      dgammadb(ipt,3,ipert)
c
     &                  + 2*wt(ipt)*xc_dvpt(ipt,igbb)*rhobmnb
     &                  +   wt(ipt)*xc_dvpt(ipt,igab)*rhoamnb
c
                     enddo
                  enddo
               enddo
            enddo
         endif ! rkstyp_sw
      else ! gradcorr_sw
         if (rkstyp_sw) then
            if (.true.) then
c
c           The following piece of code costs (FLOPs)
c           Np + Np*Nb + Np*3*N*Nb + Np*3*N*Nb*Nb*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + a*b*N**2 + a*N
c           = saving max 2.5 (realized on methane: 4 s out of 127 s)
c
c           do ipt = 1, npts
c              wrr = wt(ipt)*xc_hpt(ipt,irara)
c              do inu = 1, nao
c                 wrrn(inu) = wrr*bfn_val(ipt,inu)
c              enddo
c              do ipert = 1, npert
c                 do inu = 1, nao
c                    rbn(inu) = bfn_val(ipt,inu)*drhodb(ipt,1,ipert)
c                 enddo
c                 do inu = 1, nao
c                    do imu = 1, nao
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrn(inu)*rbn(imu)
c                    enddo
c                 enddo
c              enddo
c           enddo
c
c           The following section of code costs (FLOPs)
c           Np + 3*N*Np + 3*N*Nb*Np + 3*N*Nb*Nb*Np*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + 3*a*N**2 + a*N
c           = saving max 2.5 (realized on methane: 1.5)
c
            do ipt = 1, npts
               wrrp(ipt) = wt(ipt)*xc_hpt(ipt,irara)
            enddo
            do ipert = 1, npert
               do ipt = 1, npts
                  wrrpd(ipt,ipert) = wrrp(ipt)*drhodb(ipt,1,ipert)
               enddo
            enddo
            do ipert = 1, npert
               do inu = 1, nao
                  do ipt = 1, npts
                     wrrpdn(ipt,inu,ipert) = wrrpd(ipt,ipert)*
     &                                       bfn_val(ipt,inu)
                  enddo
               enddo
            enddo
            do ipert = 1, npert
               do inu = 1, nao
                  do imu = 1, nao
                     do ipt = 1, npts
                        ga(imu,inu,ipert) = ga(imu,inu,ipert)
     &                    + wrrpdn(ipt,inu,ipert)*bfn_val(ipt,imu)
                     enddo
                  enddo
               enddo
            enddo
            else
c
c           The following piece of code costs (FLOPs)
c           3*N*Nb*Nb*Np*5 = 15*a*b*b*N**4
c
            do ipert = 1, npert
               do inu = 1, nao
                  do imu = 1, nao
                     do ipt = 1, npts
                        ga(imu,inu,ipert) = ga(imu,inu,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irara)*
     &                      bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                      drhodb(ipt,1,ipert)
                     enddo
                  enddo
               enddo
            enddo
c
c           The cost of the following bit of code is (FLOPs)
c           Np + Np*Nb + Np*Nb*3*N + 3*N*Nb*Nb*Np*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + a*b*N**2 + a*N
c
c           do ipt = 1, npts
c              wrr = wt(ipt)*xc_hpt(ipt,irara)
c              do inu = 1, nao
c                 wrrn(ipt,inu) = wrr*bfn_val(ipt,inu)
c                 do ipert = 1, npert
c                    drn(ipt,inu,ipert) = bfn_val(ipt,inu)*
c    &                                    drhodb(ipt,1,ipert)
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrn(ipt,inu)*drn(ipt,imu,ipert)
c                    enddo
c                 enddo
c              enddo
c           enddo
            endif
         else ! rkstyp_sw
            do ipert = 1, npert
               do inu = 1, nao
                  do imu = 1, nao
                     do ipt = 1, npts
                        ga(imu,inu,ipert) = ga(imu,inu,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irara)*
     &                      bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                      drhodb(ipt,1,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irarb)*
     &                      bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                      drhodb(ipt,2,ipert)
                        gb(imu,inu,ipert) = gb(imu,inu,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                      bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                      drhodb(ipt,2,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irarb)*
     &                      bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                      drhodb(ipt,1,ipert)
                     enddo
                  enddo
               enddo
            enddo
         endif ! rkstyp_sw
      endif ! gradcorr_sw
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine lhs_dft_ao_scr(rkstyp_sw,gradcorr_sw,
     &     active_bfn_list,n_active_bfn,
     &     active_chf_pert,n_active_chf_prt,
     &     nao,npts,npert,wt,grho,
     &     drhodb,dgrhodb,dgammadb,bfn_val,bfng_val,xc_dvpt,
     &     xc_hpt,xc_dhpt,ga,gb,mxp,wrrp,wrrpd,wrrpdn,drmn,dgrn,dvn,
     &     ddr,dhnd)
      implicit none
c
c     Computes the DFT terms that contribute to the Left-Hand-Side
c     of the CPHF equations.
c
c     For performance analysis we introduce the following variables:
c     N  = the number of atoms
c     Np = a*N = the number of grid points (a = +- thousands)
c     Nb = b*N = the number of basis functions (b = 1 to 2 dozen)
c
c     Parameters
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
c     In variables
c
      logical rkstyp_sw ! .true. if closed shell
      logical gradcorr_sw
      integer n_active_bfn
      integer active_bfn_list(n_active_bfn)
      integer n_active_chf_prt
      integer active_chf_pert(n_active_chf_prt)
      integer nao, npts, npert, mxp
      real*8 wt(mxp)
      real*8 grho(mxp,2,3)
      real*8 drhodb(mxp,2,npert)
      real*8 dgrhodb(mxp,2,3,npert)
      real*8 dgammadb(mxp,3,npert)
      real*8 bfn_val(mxp,nao)
      real*8 bfng_val(mxp,nao,3)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
c     Out variables
c
      real*8 ga(nao,nao,npert)
      real*8 gb(nao,nao,npert)
c
c     Local variables
c
      integer ipert, iprt, ipt, imu, inu, im, in
      real*8 rhoamn, rhobmn, rhoamnb, rhobmnb

c     real*8 wrr, wrg, wrg2, wgg, wg
c     real*8 wrrdr, wrgdr, wrgdg, wggdg
c     real*8 grgb(n_active_bfn), dgrgb(n_active_bfn)
c     real*8 wrrdrn(n_active_bfn), wrgdrn(n_active_bfn)
c     real*8 wrgdgn(n_active_bfn), wggdgn(n_active_bfn)
c     real*8 wrrn(n_active_bfn), rbn(n_active_bfn)
      real*8 wrrp(npts), wrrpd(npts,npert)
      real*8 wrrpdn(npts,n_active_bfn,npert)
      real*8 drmn(npts,n_active_bfn,n_active_bfn)
      real*8 dgrn(npts,n_active_bfn,npert)
      real*8 dvn(npts,n_active_bfn)
      real*8 ddr(npts,npert), dhnd(npts,n_active_bfn,npert)
c
c     Code
c
      if (gradcorr_sw) then
         if (rkstyp_sw) then
c
c           The following code section costs (FLOPs)
c           3*N*Nb*Nb*Na*(13+13+22)
c           = 144*a*b*b*N**4
c
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c
c                       rhoamn = 
c    &                    (grho(ipt,1,1)*bfng_val(ipt,inu,1)+
c    &                     grho(ipt,1,2)*bfng_val(ipt,inu,2)+
c    &                     grho(ipt,1,3)*bfng_val(ipt,inu,3))*
c    &                    bfn_val(ipt,imu)+
c    &                    (grho(ipt,1,1)*bfng_val(ipt,imu,1)+
c    &                     grho(ipt,1,2)*bfng_val(ipt,imu,2)+
c    &                     grho(ipt,1,3)*bfng_val(ipt,imu,3))*
c    &                    bfn_val(ipt,inu)
c
c                       rhoamnb = 
c    &                    (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,inu,1)+
c    &                     dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,inu,2)+
c    &                     dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,inu,3))*
c    &                    bfn_val(ipt,imu)+
c    &                    (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,imu,1)+
c    &                     dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,imu,2)+
c    &                     dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,imu,3))*
c    &                    bfn_val(ipt,inu)
c
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                  + wt(ipt)*xc_hpt(ipt,irara)*
c    &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
c    &                    drhodb(ipt,1,ipert)
c    &                  + wt(ipt)*xc_dhpt(ipt,iragaa)*
c    &                    bfn_val(ipt,imu)*bfn_val(ipt,inu)*
c    &                    dgammadb(ipt,1,ipert)
c    &                  + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*rhoamn*
c    &                      drhodb(ipt,1,ipert)
c    &                  +   wt(ipt)*xc_dhpt(ipt,igaagaa)*rhoamn*
c    &                      dgammadb(ipt,1,ipert)
c    &                  +   wt(ipt)*xc_dvpt(ipt,igaa)*rhoamnb
c
c                    enddo
c                 enddo
c              enddo
c           enddo
c
c           The following section of code costs (FLOPs)
c           Nb*Np*3 + Nb*Nb*Np*13 + 3*N*Nb*Np*(5+6) + 3*N*Np*5
c           + 3*N*Nb*Nb*Np*8
c           = 3*a*b*N**2 + 15*a*N**2
c           + 13*a*b*b*N**3 + 33*a*b*N**3
c           + 24*a*b*b*N**4
c           = saving of max 6 (realized on methane 2)
c
            do inu = 1, n_active_bfn
               do ipt = 1, npts
                  dvn(ipt,inu) = wt(ipt)*xc_dvpt(ipt,igaa)*
     &                           bfn_val(ipt,inu)
               enddo
            enddo
            do inu = 1, n_active_bfn
               do imu = 1, n_active_bfn
                  do ipt = 1, npts
                     drmn(ipt,imu,inu) = 
     &                   (grho(ipt,1,1)*bfng_val(ipt,inu,1)+
     &                    grho(ipt,1,2)*bfng_val(ipt,inu,2)+
     &                    grho(ipt,1,3)*bfng_val(ipt,inu,3))*
     &                   bfn_val(ipt,imu)+
     &                   (grho(ipt,1,1)*bfng_val(ipt,imu,1)+
     &                    grho(ipt,1,2)*bfng_val(ipt,imu,2)+
     &                    grho(ipt,1,3)*bfng_val(ipt,imu,3))*
     &                   bfn_val(ipt,inu)
                  enddo
               enddo
            enddo
            do ipert = 1, n_active_chf_prt
               do inu = 1, n_active_bfn
                  do ipt = 1, npts
                     dgrn(ipt,inu,ipert) = 
     &                   (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,inu,1)+
     &                    dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,inu,2)+
     &                    dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,inu,3))
                     dhnd(ipt,inu,ipert) = wt(ipt)*(
     &                    xc_hpt(ipt,irara)*bfn_val(ipt,inu)*
     &                    drhodb(ipt,1,ipert)
     &                  + xc_dhpt(ipt,iragaa)*bfn_val(ipt,inu)*
     &                    dgammadb(ipt,1,ipert))
                  enddo
               enddo
            enddo
            do ipert = 1, n_active_chf_prt
               do ipt = 1, npts
                  ddr(ipt,ipert) = wt(ipt)*(
     &                2*xc_dhpt(ipt,iragaa)*drhodb(ipt,1,ipert)
     &              + xc_dhpt(ipt,igaagaa)*dgammadb(ipt,1,ipert))
               enddo
            enddo
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm = 1, n_lhs_atm
c              do l = 1, 3
c                 iprt=3*(iatm-1)+l
c                 ipert=lhs_atm_list(iatm)+l
                  do inu = 1, n_active_bfn
                     in = active_bfn_list(inu)
                     do imu = 1, n_active_bfn
                        im = active_bfn_list(imu)
                        do ipt = 1, npts
                           ga(im,in,ipert) = ga(im,in,ipert)
     &                     + dhnd(ipt,inu,iprt)*bfn_val(ipt,imu)
     &                     + ddr(ipt,iprt)*drmn(ipt,imu,inu)
     &                     + dvn(ipt,imu)*dgrn(ipt,inu,iprt)
     &                     + dvn(ipt,inu)*dgrn(ipt,imu,iprt)
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
         else ! rkstyp_sw
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm = 1, n_lhs_atm
c              do l = 1, 3
c                 iprt=3*(iatm-1)+l
c                 ipert=lhs_atm_list(iatm)+l
                  do inu = 1, n_active_bfn
                     in = active_bfn_list(inu)
                     do imu = 1, n_active_bfn
                        im = active_bfn_list(imu)
                        do ipt = 1, npts
c
                           rhoamn = 
     &                       (grho(ipt,1,1)*bfng_val(ipt,inu,1)+
     &                        grho(ipt,1,2)*bfng_val(ipt,inu,2)+
     &                        grho(ipt,1,3)*bfng_val(ipt,inu,3))*
     &                       bfn_val(ipt,imu)+
     &                       (grho(ipt,1,1)*bfng_val(ipt,imu,1)+
     &                        grho(ipt,1,2)*bfng_val(ipt,imu,2)+
     &                        grho(ipt,1,3)*bfng_val(ipt,imu,3))*
     &                       bfn_val(ipt,inu)
c
                           rhobmn = 
     &                       (grho(ipt,2,1)*bfng_val(ipt,inu,1)+
     &                        grho(ipt,2,2)*bfng_val(ipt,inu,2)+
     &                        grho(ipt,2,3)*bfng_val(ipt,inu,3))*
     &                       bfn_val(ipt,imu)+
     &                       (grho(ipt,2,1)*bfng_val(ipt,imu,1)+
     &                        grho(ipt,2,2)*bfng_val(ipt,imu,2)+
     &                        grho(ipt,2,3)*bfng_val(ipt,imu,3))*
     &                       bfn_val(ipt,inu)
c
                           rhoamnb = 
     &                       (dgrhodb(ipt,1,1,iprt)*bfng_val(ipt,inu,1)+
     &                        dgrhodb(ipt,1,2,iprt)*bfng_val(ipt,inu,2)+
     &                        dgrhodb(ipt,1,3,iprt)*bfng_val(ipt,inu,3))
     &                       *bfn_val(ipt,imu)+
     &                       (dgrhodb(ipt,1,1,iprt)*bfng_val(ipt,imu,1)+
     &                        dgrhodb(ipt,1,2,iprt)*bfng_val(ipt,imu,2)+
     &                        dgrhodb(ipt,1,3,iprt)*bfng_val(ipt,imu,3))
     &                       *bfn_val(ipt,inu)
c
                           rhobmnb = 
     &                       (dgrhodb(ipt,2,1,iprt)*bfng_val(ipt,inu,1)+
     &                        dgrhodb(ipt,2,2,iprt)*bfng_val(ipt,inu,2)+
     &                        dgrhodb(ipt,2,3,iprt)*bfng_val(ipt,inu,3))
     &                       *bfn_val(ipt,imu)+
     &                       (dgrhodb(ipt,2,1,iprt)*bfng_val(ipt,imu,1)+
     &                        dgrhodb(ipt,2,2,iprt)*bfng_val(ipt,imu,2)+
     &                        dgrhodb(ipt,2,3,iprt)*bfng_val(ipt,imu,3))
     &                       *bfn_val(ipt,inu)
c
                           ga(im,in,ipert) = ga(im,in,ipert)
     &                     + wt(ipt)*xc_hpt(ipt,irara)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       drhodb(ipt,1,iprt)
     &                     + wt(ipt)*xc_hpt(ipt,irarb)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       drhodb(ipt,2,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       dgammadb(ipt,1,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,iragab)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       dgammadb(ipt,2,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,iragbb)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       dgammadb(ipt,3,iprt)
c
     &                     + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*rhoamn*
     &                         drhodb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,iragab)*rhobmn*
     &                         drhodb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,irbgaa)*rhoamn*
     &                         drhodb(ipt,2,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,irbgab)*rhobmn*
     &                         drhodb(ipt,2,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igaagaa)*rhoamn*
     &                         dgammadb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igaagab)*rhobmn*
     &                         dgammadb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igaagab)*rhoamn*
     &                         dgammadb(ipt,2,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igabgab)*rhobmn*
     &                         dgammadb(ipt,2,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igaagbb)*rhoamn*
     &                         dgammadb(ipt,3,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igabgbb)*rhobmn*
     &                         dgammadb(ipt,3,iprt)
c
     &                     + 2*wt(ipt)*xc_dvpt(ipt,igaa)*rhoamnb
     &                     +   wt(ipt)*xc_dvpt(ipt,igab)*rhobmnb
c
                           gb(im,in,ipert) = gb(im,in,ipert)
     &                     + wt(ipt)*xc_hpt(ipt,irarb)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       drhodb(ipt,1,iprt)
     &                     + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       drhodb(ipt,2,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgaa)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       dgammadb(ipt,1,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgab)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       dgammadb(ipt,2,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgbb)*
     &                       bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                       dgammadb(ipt,3,iprt)
c
     &                     + 2*wt(ipt)*xc_dhpt(ipt,iragbb)*rhobmn*
     &                         drhodb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,iragab)*rhoamn*
     &                         drhodb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,irbgbb)*rhobmn*
     &                         drhodb(ipt,2,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,irbgab)*rhoamn*
     &                         drhodb(ipt,2,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igaagbb)*rhobmn*
     &                         dgammadb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igaagab)*rhoamn*
     &                         dgammadb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igabgbb)*rhobmn*
     &                         dgammadb(ipt,2,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igabgab)*rhoamn*
     &                         dgammadb(ipt,2,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*rhobmn*
     &                         dgammadb(ipt,3,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igabgbb)*rhoamn*
     &                         dgammadb(ipt,3,iprt)
c
     &                     + 2*wt(ipt)*xc_dvpt(ipt,igbb)*rhobmnb
     &                     +   wt(ipt)*xc_dvpt(ipt,igab)*rhoamnb
c
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
         endif ! rkstyp_sw
      else ! gradcorr_sw
         if (rkstyp_sw) then
            if (.true.) then
c
c           The following piece of code costs (FLOPs)
c           Np + Np*Nb + Np*3*N*Nb + Np*3*N*Nb*Nb*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + a*b*N**2 + a*N
c           = saving max 2.5 (realized on methane: 4 s out of 127 s)
c
c           do ipt = 1, npts
c              wrr = wt(ipt)*xc_hpt(ipt,irara)
c              do inu = 1, nao
c                 wrrn(inu) = wrr*bfn_val(ipt,inu)
c              enddo
c              do ipert = 1, npert
c                 do inu = 1, nao
c                    rbn(inu) = bfn_val(ipt,inu)*drhodb(ipt,1,ipert)
c                 enddo
c                 do inu = 1, nao
c                    do imu = 1, nao
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrn(inu)*rbn(imu)
c                    enddo
c                 enddo
c              enddo
c           enddo
c
c           The following section of code costs (FLOPs)
c           Np + 3*N*Np + 3*N*Nb*Np + 3*N*Nb*Nb*Np*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + 3*a*N**2 + a*N
c           = saving max 2.5 (realized on methane: 1.5)
c
            do ipt = 1, npts
               wrrp(ipt) = wt(ipt)*xc_hpt(ipt,irara)
            enddo
            do ipert = 1, n_active_chf_prt
               do ipt = 1, npts
                  wrrpd(ipt,ipert) = wrrp(ipt)*drhodb(ipt,1,ipert)
               enddo
            enddo
            do ipert = 1, n_active_chf_prt
               do inu = 1, n_active_bfn
                  do ipt = 1, npts
                     wrrpdn(ipt,inu,ipert) = wrrpd(ipt,ipert)*
     &                                       bfn_val(ipt,inu)
                  enddo
               enddo
            enddo
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm = 1, n_lhs_atm
c              do l = 1, 3
c                 iprt=3*(iatm-1)+l
c                 ipert=lhs_atm_list(iatm)+l
                  do inu = 1, n_active_bfn
                     in = active_bfn_list(inu)
                     do imu = 1, n_active_bfn
                        im = active_bfn_list(imu)
                        do ipt = 1, npts
                           ga(im,in,ipert) = ga(im,in,ipert)
     &                       + wrrpdn(ipt,inu,iprt)*bfn_val(ipt,imu)
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
            else
c
c           The following piece of code costs (FLOPs)
c           3*N*Nb*Nb*Np*5 = 15*a*b*b*N**4
c
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm = 1, n_lhs_atm
c              do l = 1, 3
c                 iprt=3*(iatm-1)+l
c                 ipert=lhs_atm_list(iatm)+l
                  do inu = 1, n_active_bfn
                     in = active_bfn_list(inu)
                     do imu = 1, n_active_bfn
                        im = active_bfn_list(imu)
                        do ipt = 1, npts
                           ga(im,in,ipert) = ga(im,in,ipert)
     &                       + wt(ipt)*xc_hpt(ipt,irara)*
     &                         bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                         drhodb(ipt,1,iprt)
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
c
c           The cost of the following bit of code is (FLOPs)
c           Np + Np*Nb + Np*Nb*3*N + 3*N*Nb*Nb*Np*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + a*b*N**2 + a*N
c
c           do ipt = 1, npts
c              wrr = wt(ipt)*xc_hpt(ipt,irara)
c              do inu = 1, nao
c                 wrrn(ipt,inu) = wrr*bfn_val(ipt,inu)
c                 do ipert = 1, npert
c                    drn(ipt,inu,ipert) = bfn_val(ipt,inu)*
c    &                                    drhodb(ipt,1,ipert)
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrn(ipt,inu)*drn(ipt,imu,ipert)
c                    enddo
c                 enddo
c              enddo
c           enddo
            endif
         else ! rkstyp_sw
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm = 1, n_lhs_atm
c              do l = 1, 3
c                 iprt=3*(iatm-1)+l
c                 ipert=lhs_atm_list(iatm)+l
                  do inu = 1, n_active_bfn
                     in = active_bfn_list(inu)
                     do imu = 1, n_active_bfn
                        im = active_bfn_list(imu)
                        do ipt = 1, npts
                           ga(im,in,ipert) = ga(im,in,ipert)
     &                       + wt(ipt)*xc_hpt(ipt,irara)*
     &                         bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                         drhodb(ipt,1,iprt)
     &                       + wt(ipt)*xc_hpt(ipt,irarb)*
     &                         bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                         drhodb(ipt,2,iprt)
                           gb(im,in,ipert) = gb(im,in,ipert)
     &                       + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                         bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                         drhodb(ipt,2,iprt)
     &                       + wt(ipt)*xc_hpt(ipt,irarb)*
     &                         bfn_val(ipt,imu)*bfn_val(ipt,inu)*
     &                         drhodb(ipt,1,iprt)
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
         endif ! rkstyp_sw
      endif ! gradcorr_sw
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_lhs_mo(npts,nvec,naocc,npert,dentol,wt,xc_hpt,
     &           amo_val,uamat,ga,drhodb,t,t2,t3,twt)
      implicit none
c
c     Compute the left-hand-side in MO-basis for the local density
c     closed shell case.
c
c     Inputs
c
      integer npts
      integer npert
      integer nvec
      integer naocc
      real*8 wt(npts)
      real*8 xc_hpt(npts)
      real*8 uamat(naocc,nvec-naocc,npert)
      real*8 amo_val(npts,nvec)
      real*8 dentol
c
c     Outputs
c
      real*8 ga(naocc,nvec-naocc,npert)
c
c     Workspace
c
      real*8 drhodb(npts)
      real*8 t(npts)
      real*8 t2(npts)
      real*8 t3(npts)
      real*8 twt(npts)
c
c     Local variables
c
      integer nvirt
      integer inu, imu
      integer ipt, iprt
      integer i, j
      real*8 uval
c
c     Code
c
      nvirt = nvec-naocc
      do ipt=1,npts
         twt(ipt)=2.0d0*wt(ipt)*xc_hpt(ipt)
      enddo
      do iprt=1,npert
         call aclear_dp(drhodb,npts,0.0d0)
         do i=1,nvirt
            call aclear_dp(t,npts,0.0d0)
            do j=1,naocc
               uval=uamat(j,i,iprt)
               if (abs(uval).gt.dentol) then
                  do ipt=1,npts
                     t(ipt)=t(ipt)+amo_val(ipt,j)*uval
                  enddo
               endif
            enddo
            do ipt=1,npts
               drhodb(ipt)=drhodb(ipt)
     &                    +amo_val(ipt,naocc+i)*t(ipt)
            enddo
         enddo
         do ipt=1,npts
            t2(ipt)=twt(ipt)*drhodb(ipt)
         enddo
         do inu = 1, nvirt
            do ipt=1,npts
               t3(ipt)=t2(ipt)*amo_val(ipt,naocc+inu)
            enddo
            do imu = 1, naocc
               do ipt = 1, npts
                  ga(imu,inu,iprt) = ga(imu,inu,iprt)
     &              + t3(ipt)*amo_val(ipt,imu)
               enddo
            enddo
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_lhs_mo_scr(npts,nvec,naocc,npert,
     &           active_chf_pert,n_active_chf_prt,
     &           dentol,wt,xc_hpt,
     &           amo_val,uamat,ga,drhodb,t,t2,t3,twt)
      implicit none
c
c     Compute the left-hand-side in MO-basis for the local density
c     closed shell case.
c
c     Inputs
c
      integer npts
      integer npert
      integer nvec
      integer naocc
      integer n_active_chf_prt
      integer active_chf_pert(n_active_chf_prt)
      real*8 wt(npts)
      real*8 xc_hpt(npts)
      real*8 uamat(naocc,nvec-naocc,npert)
      real*8 amo_val(npts,nvec)
      real*8 dentol
c
c     Outputs
c
      real*8 ga(naocc,nvec-naocc,npert)
c
c     Workspace
c
      real*8 drhodb(npts)
      real*8 t(npts)
      real*8 t2(npts)
      real*8 t3(npts)
      real*8 twt(npts)
c
c     Local variables
c
      integer nvirt
      integer inu, imu
      integer ipt, ipert, iprt
      integer i, j
      real*8 uval
c
c     Code
c
      nvirt = nvec-naocc
      do ipt=1,npts
         twt(ipt)=2.0d0*wt(ipt)*xc_hpt(ipt)
      enddo
      do iprt=1,n_active_chf_prt
         ipert=active_chf_pert(iprt)
c     do iatm=1,n_lhs_atm
c        do ic=1,3
c           ipert=lhs_atm_list(iatm)+ic
            call aclear_dp(drhodb,npts,0.0d0)
            do i=1,nvirt
               call aclear_dp(t,npts,0.0d0)
               do j=1,naocc
                  uval=uamat(j,i,ipert)
                  if (abs(uval).gt.dentol) then
                     do ipt=1,npts
                        t(ipt)=t(ipt)+amo_val(ipt,j)*uval
                     enddo
                  endif
               enddo
               do ipt=1,npts
                  drhodb(ipt)=drhodb(ipt)
     &                       +amo_val(ipt,naocc+i)*t(ipt)
               enddo
            enddo
            do ipt=1,npts
               t2(ipt)=twt(ipt)*drhodb(ipt)
            enddo
            do inu = 1, nvirt
               do ipt=1,npts
                  t3(ipt)=t2(ipt)*amo_val(ipt,naocc+inu)
               enddo
               do imu = 1, naocc
                  do ipt = 1, npts
                     ga(imu,inu,ipert) = ga(imu,inu,ipert)
     &                 + t3(ipt)*amo_val(ipt,imu)
                  enddo
               enddo
            enddo
c        enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_lhs_mo_gga(mxp,npts,npert,nvec,naocc,dentol,
     &           wt,xc_dvpt,xc_hpt,xc_dhpt,
     &           uamat,amo_val,amo_grad,grho,drhodb,dgrhodb,t,dt,
     &           rhoan,rhoanb,dgammadb,ga)
      implicit none
c
c     Compute the left-hand-side in MO-basis for a GGA functional
c     in a closed shell calculation.
c
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
c     Inputs
c
      integer mxp
      integer npts
      integer npert
      integer nvec
      integer naocc
      real*8 wt(npts)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
      real*8 uamat(naocc,nvec-naocc,npert)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 grho(mxp,2,3)
      real*8 dentol
c
c     Outputs
c
      real*8 ga(naocc,nvec-naocc,npert)
c
c     Workspace
c
      real*8 drhodb(npts)
      real*8 dgrhodb(npts,3)
      real*8 t(npts)
      real*8 dt(npts,3)
      real*8 rhoan(npts,nvec)
      real*8 rhoanb(npts,nvec)
      real*8 dgammadb(npts)
c
c     Local variables
c
      integer nvirt
      integer iprt, ipt, inu, imu
      integer k
      real*8 uval
c
c     Code
c
      nvirt = nvec-naocc
c
      call aclear_dp(rhoan,npts*nvec,0.0d0)
      do k=1,3
         do inu=1,nvec
            do ipt=1,npts
               rhoan(ipt,inu)=rhoan(ipt,inu)
     &                       +grho(ipt,1,k)*amo_grad(ipt,inu,k)
            enddo
         enddo
      enddo
      do inu=1,nvec
         do ipt=1,npts
            rhoan(ipt,inu)=0.5d0*rhoan(ipt,inu)
         enddo
      enddo
c
      do iprt=1,npert
c
c        Construct drhodb and dgrhodb
c
         call aclear_dp(drhodb,npts,0.0d0)
         call aclear_dp(dgrhodb,3*npts,0.0d0)
         do inu=1,nvirt
            call aclear_dp(t,npts,0.0d0)
            call aclear_dp(dt,3*npts,0.0d0)
            do imu=1,naocc
               uval=uamat(imu,inu,iprt)
               if (abs(uval).gt.dentol) then
                  do ipt=1,npts
                     t(ipt)=t(ipt)+amo_val(ipt,imu)*uval
                  enddo
                  do k=1,3
                     do ipt=1,npts
                        dt(ipt,k)=dt(ipt,k)+amo_grad(ipt,imu,k)*uval
                     enddo
                  enddo
               endif
            enddo
            do ipt=1,npts
               drhodb(ipt)=drhodb(ipt)
     &                    +amo_val(ipt,naocc+inu)*t(ipt)
            enddo
            do k=1,3
               do ipt=1,npts
                  dgrhodb(ipt,k)=dgrhodb(ipt,k)
     &                          +amo_grad(ipt,naocc+inu,k)*
     &                           t(ipt)
     &                          +amo_val(ipt,naocc+inu)*
     &                           dt(ipt,k)
               enddo
            enddo
         enddo
         do ipt=1,npts
            drhodb(ipt)=2.0d0*drhodb(ipt)
         enddo
         do k=1,3
            do ipt=1,npts
               dgrhodb(ipt,k)=2.0d0*dgrhodb(ipt,k)
            enddo
         enddo
c
c        Construct dgammadb and rhoanb
c
         call aclear_dp(dgammadb,npts,0.0d0)
         do k=1,3
            do ipt=1,npts
               dgammadb(ipt)=dgammadb(ipt)
     &                      +grho(ipt,1,k)*dgrhodb(ipt,k)
            enddo
         enddo
c        do ipt=1,npts
c           dgammadb(ipt)=2.0d0*dgammadb(ipt)
c        enddo
         call aclear_dp(rhoanb,npts*nvec,0.0d0)
         do k=1,3
            do inu=1,nvec
               do ipt=1,npts
                  rhoanb(ipt,inu)=rhoanb(ipt,inu)
     &                           +dgrhodb(ipt,k)*amo_grad(ipt,inu,k)
               enddo
            enddo
         enddo
c
c        Construct the final expression
c
         do ipt=1,npts
            dt(ipt,1) = wt(ipt)*(xc_hpt(ipt,irara)*drhodb(ipt)
     &                +          xc_dhpt(ipt,iragaa)*dgammadb(ipt))
            dt(ipt,2) = wt(ipt)*(xc_dhpt(ipt,iragaa)*drhodb(ipt)*2.0d0
     &                +          xc_dhpt(ipt,igaagaa)*dgammadb(ipt))
            dt(ipt,3) = wt(ipt)*xc_dvpt(ipt,igaa)
         enddo
         do inu=1,nvirt
            do ipt=1,npts
               t(ipt)=dt(ipt,1)*amo_val(ipt,naocc+inu)
     &               +dt(ipt,2)*rhoan(ipt,naocc+inu)
     &               +dt(ipt,3)*rhoanb(ipt,naocc+inu)
            enddo
            do imu=1,naocc
               do ipt=1,npts
                  ga(imu,inu,iprt)=ga(imu,inu,iprt)
     &                            +t(ipt)*amo_val(ipt,imu)
               enddo
            enddo
         enddo
         if (.false.) then
c
c           This section uses more FLOPs but therefore memory access
c           is more regular.
c
c           This section led to a cost of:
c               102.60 secs on Methane (6-31G*).
c              1143.75 secs on Ethane (6-31G*).
c
            do inu=1,nvirt
               do ipt=1,npts
                  t(ipt)=dt(ipt,2)*amo_val(ipt,naocc+inu)
               enddo
               do imu=1,naocc
                  do ipt=1,npts
                     ga(imu,inu,iprt)=ga(imu,inu,iprt)
     &                               +t(ipt)*rhoan(ipt,imu)
                  enddo
               enddo
               do ipt=1,npts
                  t(ipt)=dt(ipt,3)*amo_val(ipt,naocc+inu)
               enddo
               do imu=1,naocc
                  do ipt=1,npts
                     ga(imu,inu,iprt)=ga(imu,inu,iprt)
     &                               +t(ipt)*rhoanb(ipt,imu)
                  enddo
               enddo
            enddo
         else
c
c           This section uses less FLOPs but therefore memory access
c           requires more jumps
c
c           This section led to a cost of:
c                84.05 secs on Methane (6-31G*).
c               862.62 secs on Ethane (6-31G*).
c
            do imu=1,naocc
               do ipt=1,npts
                  t(ipt)=dt(ipt,2)*rhoan(ipt,imu)
     &                  +dt(ipt,3)*rhoanb(ipt,imu)
               enddo
               do inu=1,nvirt
                  do ipt=1,npts
                     ga(imu,inu,iprt)=ga(imu,inu,iprt)
     &                               +t(ipt)*amo_val(ipt,naocc+inu)
                  enddo
               enddo
            enddo
         endif
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_lhs_mo_gga_scr(mxp,npts,npert,nvec,naocc,
     &           active_chf_pert,n_active_chf_prt,
     &           dentol,wt,xc_dvpt,xc_hpt,xc_dhpt,
     &           uamat,amo_val,amo_grad,grho,drhodb,dgrhodb,t,dt,
     &           rhoan,rhoanb,dgammadb,ga)
      implicit none
c
c     Compute the left-hand-side in MO-basis for a GGA functional
c     in a closed shell calculation.
c
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
c     Inputs
c
      integer mxp
      integer npts
      integer npert
      integer nvec
      integer naocc
      integer n_active_chf_prt
      integer active_chf_pert(n_active_chf_prt)
      real*8 wt(npts)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
      real*8 uamat(naocc,nvec-naocc,npert)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 grho(mxp,2,3)
      real*8 dentol
c
c     Outputs
c
      real*8 ga(naocc,nvec-naocc,npert)
c
c     Workspace
c
      real*8 drhodb(npts)
      real*8 dgrhodb(npts,3)
      real*8 t(npts)
      real*8 dt(npts,3)
      real*8 rhoan(npts,nvec)
      real*8 rhoanb(npts,nvec)
      real*8 dgammadb(npts)
c
c     Local variables
c
      integer nvirt
      integer ipert, iprt, ipt, inu, imu
      integer k
      real*8 uval
c
c     Code
c
      nvirt = nvec-naocc
c
      call aclear_dp(rhoan,npts*nvec,0.0d0)
      do k=1,3
         do inu=1,nvec
            do ipt=1,npts
               rhoan(ipt,inu)=rhoan(ipt,inu)
     &                       +grho(ipt,1,k)*amo_grad(ipt,inu,k)
            enddo
         enddo
      enddo
      do inu=1,nvec
         do ipt=1,npts
            rhoan(ipt,inu)=0.5d0*rhoan(ipt,inu)
         enddo
      enddo
c
      do iprt=1,n_active_chf_prt
         ipert=active_chf_pert(iprt)
c     do iatm=1,n_lhs_atm
c        do ic=1,3
c           ipert=lhs_atm_list(iatm)+ic
c
c           Construct drhodb and dgrhodb
c
            call aclear_dp(drhodb,npts,0.0d0)
            call aclear_dp(dgrhodb,3*npts,0.0d0)
            do inu=1,nvirt
               call aclear_dp(t,npts,0.0d0)
               call aclear_dp(dt,3*npts,0.0d0)
               do imu=1,naocc
                  uval=uamat(imu,inu,ipert)
                  if (abs(uval).gt.dentol) then
                     do ipt=1,npts
                        t(ipt)=t(ipt)+amo_val(ipt,imu)*uval
                     enddo
                     do k=1,3
                        do ipt=1,npts
                           dt(ipt,k)=dt(ipt,k)+amo_grad(ipt,imu,k)*uval
                        enddo
                     enddo
                  endif
               enddo
               do ipt=1,npts
                  drhodb(ipt)=drhodb(ipt)
     &                       +amo_val(ipt,naocc+inu)*t(ipt)
               enddo
               do k=1,3
                  do ipt=1,npts
                     dgrhodb(ipt,k)=dgrhodb(ipt,k)
     &                             +amo_grad(ipt,naocc+inu,k)*
     &                              t(ipt)
     &                             +amo_val(ipt,naocc+inu)*
     &                              dt(ipt,k)
                  enddo
               enddo
            enddo
            do ipt=1,npts
               drhodb(ipt)=2.0d0*drhodb(ipt)
            enddo
            do k=1,3
               do ipt=1,npts
                  dgrhodb(ipt,k)=2.0d0*dgrhodb(ipt,k)
               enddo
            enddo
c
c           Construct dgammadb and rhoanb
c
            call aclear_dp(dgammadb,npts,0.0d0)
            do k=1,3
               do ipt=1,npts
                  dgammadb(ipt)=dgammadb(ipt)
     &                         +grho(ipt,1,k)*dgrhodb(ipt,k)
               enddo
            enddo
c           do ipt=1,npts
c              dgammadb(ipt)=2.0d0*dgammadb(ipt)
c           enddo
            call aclear_dp(rhoanb,npts*nvec,0.0d0)
            do k=1,3
               do inu=1,nvec
                  do ipt=1,npts
                     rhoanb(ipt,inu)=rhoanb(ipt,inu)
     &                              +dgrhodb(ipt,k)*amo_grad(ipt,inu,k)
                  enddo
               enddo
            enddo
c
c           Construct the final expression
c
            do ipt=1,npts
               dt(ipt,1) = wt(ipt)*(xc_hpt(ipt,irara)*drhodb(ipt)
     &                   +          xc_dhpt(ipt,iragaa)*dgammadb(ipt))
               dt(ipt,2) = wt(ipt)*(xc_dhpt(ipt,iragaa)*drhodb(ipt)*2
     &                   +          xc_dhpt(ipt,igaagaa)*dgammadb(ipt))
               dt(ipt,3) = wt(ipt)*xc_dvpt(ipt,igaa)
            enddo
            do inu=1,nvirt
               do ipt=1,npts
                  t(ipt)=dt(ipt,1)*amo_val(ipt,naocc+inu)
     &                  +dt(ipt,2)*rhoan(ipt,naocc+inu)
     &                  +dt(ipt,3)*rhoanb(ipt,naocc+inu)
               enddo
               do imu=1,naocc
                  do ipt=1,npts
                     ga(imu,inu,ipert)=ga(imu,inu,ipert)
     &                                +t(ipt)*amo_val(ipt,imu)
                  enddo
               enddo
            enddo
            if (.false.) then
c
c              This section uses more FLOPs but therefore memory access
c              is more regular.
c
c              This section led to a cost of:
c                  102.60 secs on Methane (6-31G*).
c                 1143.75 secs on Ethane (6-31G*).
c
               do inu=1,nvirt
                  do ipt=1,npts
                     t(ipt)=dt(ipt,2)*amo_val(ipt,naocc+inu)
                  enddo
                  do imu=1,naocc
                     do ipt=1,npts
                        ga(imu,inu,ipert)=ga(imu,inu,ipert)
     &                                   +t(ipt)*rhoan(ipt,imu)
                     enddo
                  enddo
                  do ipt=1,npts
                     t(ipt)=dt(ipt,3)*amo_val(ipt,naocc+inu)
                  enddo
                  do imu=1,naocc
                     do ipt=1,npts
                        ga(imu,inu,ipert)=ga(imu,inu,ipert)
     &                                   +t(ipt)*rhoanb(ipt,imu)
                     enddo
                  enddo
               enddo
            else
c
c              This section uses less FLOPs but therefore memory access
c              requires more jumps
c
c              This section led to a cost of:
c                   84.05 secs on Methane (6-31G*).
c                  862.62 secs on Ethane (6-31G*).
c
               do imu=1,naocc
                  do ipt=1,npts
                     t(ipt)=dt(ipt,2)*rhoan(ipt,imu)
     &                     +dt(ipt,3)*rhoanb(ipt,imu)
                  enddo
                  do inu=1,nvirt
                     do ipt=1,npts
                        ga(imu,inu,ipert)=ga(imu,inu,ipert)
     &                                   +t(ipt)*amo_val(ipt,naocc+inu)
                     enddo
                  enddo
               enddo
            endif
c        enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_rhs_mo(npts,nvec,naocc,npert,dentol,wt,xc_hpt,
     &           amo_val,samat,ga,drhodb,t,t2,t3,twt)
      implicit none
c
c     Compute the right-hand-side in MO-basis for the local density
c     closed shell case.
c
c     Inputs
c
      integer npts
      integer npert
      integer nvec
      integer naocc
      real*8 wt(npts)
      real*8 xc_hpt(npts)
      real*8 samat(nvec*(nvec+1)/2,npert)
      real*8 amo_val(npts,nvec)
      real*8 dentol
c
c     Outputs
c
      real*8 ga(naocc,nvec-naocc,npert)
c
c     Workspace
c
      real*8 drhodb(npts)
      real*8 t(npts)
      real*8 t2(npts)
      real*8 t3(npts)
      real*8 twt(npts)
c
c     Local variables
c
      integer nvirt
      integer inu, imu
      integer ipt, iprt
      integer i, j, n
      real*8 sval
c
c     Code
c
      nvirt = nvec-naocc
      do ipt=1,npts
         twt(ipt)=2.0d0*wt(ipt)*xc_hpt(ipt)
      enddo
      do iprt=1,npert
         call aclear_dp(drhodb,npts,0.0d0)
         n = 0
         do i=1,naocc
            call aclear_dp(t,npts,0.0d0)
            do j=1,i-1
               n=n+1
               sval=samat(n,iprt)
               if (abs(sval).gt.dentol) then
                  do ipt=1,npts
                     t(ipt)=t(ipt)-amo_val(ipt,j)*sval
                  enddo
               endif
            enddo
            n=n+1
            sval=samat(n,iprt)
            if (abs(sval).gt.dentol) then
               sval=0.5d0*sval
               do ipt=1,npts
                  t(ipt)=t(ipt)-amo_val(ipt,i)*sval
               enddo
            endif
            do ipt=1,npts
               drhodb(ipt)=drhodb(ipt)
     &                    +amo_val(ipt,i)*t(ipt)
            enddo
         enddo
         do ipt=1,npts
            t2(ipt)=twt(ipt)*drhodb(ipt)
         enddo
         do inu = 1, nvirt
            do ipt=1,npts
               t3(ipt)=t2(ipt)*amo_val(ipt,naocc+inu)
            enddo
            do imu = 1, naocc
               do ipt = 1, npts
                  ga(imu,inu,iprt) = ga(imu,inu,iprt)
     &              + t3(ipt)*amo_val(ipt,imu)
               enddo
            enddo
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_rhs_mo_scr(npts,nvec,naocc,npert,
     &           active_chf_pert,n_active_chf_prt,
     &           dentol,wt,xc_hpt,
     &           amo_val,samat,ga,drhodb,t,t2,t3,twt)
      implicit none
c
c     Compute the left-hand-side in MO-basis for the local density
c     closed shell case.
c
c     Inputs
c
      integer npts
      integer npert
      integer nvec
      integer naocc
      integer n_active_chf_prt
      integer active_chf_pert(n_active_chf_prt)
      real*8 wt(npts)
      real*8 xc_hpt(npts)
      real*8 samat(nvec*(nvec+1)/2,npert)
      real*8 amo_val(npts,nvec)
      real*8 dentol
c
c     Outputs
c
      real*8 ga(naocc,nvec-naocc,npert)
c
c     Workspace
c
      real*8 drhodb(npts)
      real*8 t(npts)
      real*8 t2(npts)
      real*8 t3(npts)
      real*8 twt(npts)
c
c     Local variables
c
      integer nvirt
      integer inu, imu
      integer ipt, ipert, iprt
      integer i, j, n
      real*8 sval
c
c     Code
c
      nvirt = nvec-naocc
      do ipt=1,npts
         twt(ipt)=2.0d0*wt(ipt)*xc_hpt(ipt)
      enddo
      do iprt=1,n_active_chf_prt
         ipert=active_chf_pert(iprt)
c     do iatm=1,n_lhs_atm
c        do ic=1,3
c           ipert=lhs_atm_list(iatm)+ic
            call aclear_dp(drhodb,npts,0.0d0)
            n=0
            do i=1,naocc
               call aclear_dp(t,npts,0.0d0)
               do j=1,i-1
                  n=n+1
                  sval=samat(n,ipert)
                  if (abs(sval).gt.dentol) then
                     do ipt=1,npts
                        t(ipt)=t(ipt)-amo_val(ipt,j)*sval
                     enddo
                  endif
               enddo
               n=n+1
               sval=samat(n,ipert)
               if (abs(sval).gt.dentol) then
                  sval=0.5d0*sval
                  do ipt=1,npts
                     t(ipt)=t(ipt)-amo_val(ipt,i)*sval
                  enddo
               endif
               do ipt=1,npts
                  drhodb(ipt)=drhodb(ipt)
     &                       +amo_val(ipt,i)*t(ipt)
               enddo
            enddo
            do ipt=1,npts
               t2(ipt)=twt(ipt)*drhodb(ipt)
            enddo
            do inu = 1, nvirt
               do ipt=1,npts
                  t3(ipt)=t2(ipt)*amo_val(ipt,naocc+inu)
               enddo
               do imu = 1, naocc
                  do ipt = 1, npts
                     ga(imu,inu,ipert) = ga(imu,inu,ipert)
     &                 + t3(ipt)*amo_val(ipt,imu)
                  enddo
               enddo
            enddo
c        enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_rhs_mo_gga(mxp,npts,npert,nvec,naocc,dentol,
     &           wt,xc_dvpt,xc_hpt,xc_dhpt,
     &           samat,amo_val,amo_grad,grho,drhodb,dgrhodb,t,dt,
     &           rhoan,rhoanb,dgammadb,ga)
      implicit none
c
c     Compute the left-hand-side in MO-basis for a GGA functional
c     in a closed shell calculation.
c
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
c     Inputs
c
      integer mxp
      integer npts
      integer npert
      integer nvec
      integer naocc
      real*8 wt(npts)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
      real*8 samat(nvec*(nvec+1)/2,npert)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 grho(mxp,2,3)
      real*8 dentol
c
c     Outputs
c
      real*8 ga(naocc,nvec-naocc,npert)
c
c     Workspace
c
      real*8 drhodb(npts)
      real*8 dgrhodb(npts,3)
      real*8 t(npts)
      real*8 dt(npts,3)
      real*8 rhoan(npts,nvec)
      real*8 rhoanb(npts,nvec)
      real*8 dgammadb(npts)
c
c     Local variables
c
      integer nvirt
      integer iprt, ipt, inu, imu
      integer k, n, n1
      real*8 sval
c
c     Code
c
      nvirt = nvec-naocc
c
      call aclear_dp(rhoan,npts*nvec,0.0d0)
      do k=1,3
         do inu=1,nvec
            do ipt=1,npts
               rhoan(ipt,inu)=rhoan(ipt,inu)
     &                       +grho(ipt,1,k)*amo_grad(ipt,inu,k)
            enddo
         enddo
      enddo
cDEBUG
      do inu=1,nvec
         do ipt=1,npts
            rhoan(ipt,inu)=0.5d0*rhoan(ipt,inu)
         enddo
      enddo
cDEBUG
c
      do iprt=1,npert
c
c        Construct drhodb and dgrhodb
c
         call aclear_dp(drhodb,npts,0.0d0)
         call aclear_dp(dgrhodb,3*npts,0.0d0)
         n=0
         do inu=1,naocc
            call aclear_dp(t,npts,0.0d0)
            do imu=1,inu
               n=n+1
               sval=samat(n,iprt)
               if (abs(sval).gt.dentol) then
                  do ipt=1,npts
                     t(ipt)=t(ipt)-amo_val(ipt,imu)*sval
                  enddo
               endif
            enddo
            n1=n
            do imu=inu+1,naocc
               n1=n1+imu-1
               sval=samat(n1,iprt)
               if (abs(sval).gt.dentol) then
                  do ipt=1,npts
                     t(ipt)=t(ipt)-amo_val(ipt,imu)*sval
                  enddo
               endif
            enddo
            do ipt=1,npts
               drhodb(ipt)=drhodb(ipt)
     &                    +amo_val(ipt,inu)*t(ipt)
               t(ipt)=2.0d0*t(ipt)
            enddo
            do k=1,3
               do ipt=1,npts
                  dgrhodb(ipt,k)=dgrhodb(ipt,k)
     &                          +amo_grad(ipt,inu,k)*t(ipt)
               enddo
            enddo
         enddo
c        do ipt=1,npts
c           drhodb(ipt)=2.0d0*drhodb(ipt)
c        enddo
c        do k=1,3
c           do ipt=1,npts
c              dgrhodb(ipt,k)=2.0d0*dgrhodb(ipt,k)
c           enddo
c        enddo
c
c        Construct dgammadb and rhoanb
c
         do ipt=1,npts
c           dgammadb(ipt)=2.0d0*
            dgammadb(ipt)=
     &         (grho(ipt,1,1)*dgrhodb(ipt,1)
     &         +grho(ipt,1,2)*dgrhodb(ipt,2)
     &         +grho(ipt,1,3)*dgrhodb(ipt,3))
         enddo
c
         call aclear_dp(rhoanb,npts*nvec,0.0d0)
         do k=1,3
            do inu=1,nvec
               do ipt=1,npts
                  rhoanb(ipt,inu)=rhoanb(ipt,inu)
     &                           +dgrhodb(ipt,k)*amo_grad(ipt,inu,k)
               enddo
            enddo
         enddo
c
c        Construct the final expression
c
         do ipt=1,npts
            dt(ipt,1) = wt(ipt)*(xc_hpt(ipt,irara)*drhodb(ipt)
     &                +          xc_dhpt(ipt,iragaa)*dgammadb(ipt))
            dt(ipt,2) = wt(ipt)*(xc_dhpt(ipt,iragaa)*drhodb(ipt)*2.0d0
     &                +          xc_dhpt(ipt,igaagaa)*dgammadb(ipt))
            dt(ipt,3) = wt(ipt)*xc_dvpt(ipt,igaa)
         enddo
         do inu=1,nvirt
            do ipt=1,npts
               t(ipt)=dt(ipt,1)*amo_val(ipt,naocc+inu)
     &               +dt(ipt,2)*rhoan(ipt,naocc+inu)
     &               +dt(ipt,3)*rhoanb(ipt,naocc+inu)
            enddo
            do imu=1,naocc
               do ipt=1,npts
                  ga(imu,inu,iprt)=ga(imu,inu,iprt)
     &                            +t(ipt)*amo_val(ipt,imu)
               enddo
            enddo
         enddo
         if (.false.) then
c
c           This section uses more FLOPs but therefore memory access
c           is more regular.
c
c           This section led to a cost of:
c               102.60 secs on Methane (6-31G*).
c              1143.75 secs on Ethane (6-31G*).
c
            do inu=1,nvirt
               do ipt=1,npts
                  t(ipt)=dt(ipt,2)*amo_val(ipt,naocc+inu)
               enddo
               do imu=1,naocc
                  do ipt=1,npts
                     ga(imu,inu,iprt)=ga(imu,inu,iprt)
     &                               +t(ipt)*rhoan(ipt,imu)
                  enddo
               enddo
               do ipt=1,npts
                  t(ipt)=dt(ipt,3)*amo_val(ipt,naocc+inu)
               enddo
               do imu=1,naocc
                  do ipt=1,npts
                     ga(imu,inu,iprt)=ga(imu,inu,iprt)
     &                               +t(ipt)*rhoanb(ipt,imu)
                  enddo
               enddo
            enddo
         else
c
c           This section uses less FLOPs but therefore memory access
c           requires more jumps
c
c           This section led to a cost of:
c                84.05 secs on Methane (6-31G*).
c               862.62 secs on Ethane (6-31G*).
c
            do imu=1,naocc
               do ipt=1,npts
                  t(ipt)=dt(ipt,2)*rhoan(ipt,imu)
     &                  +dt(ipt,3)*rhoanb(ipt,imu)
               enddo
               do inu=1,nvirt
                  do ipt=1,npts
                     ga(imu,inu,iprt)=ga(imu,inu,iprt)
     &                               +t(ipt)*amo_val(ipt,naocc+inu)
                  enddo
               enddo
            enddo
         endif
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rks_rhs_mo_gga_scr(mxp,npts,npert,nvec,naocc,
     &           active_chf_pert,n_active_chf_prt,
     &           dentol,wt,xc_dvpt,xc_hpt,xc_dhpt,
     &           samat,amo_val,amo_grad,grho,drhodb,dgrhodb,t,dt,
     &           rhoan,rhoanb,dgammadb,ga)
      implicit none
c
c     Compute the left-hand-side in MO-basis for a GGA functional
c     in a closed shell calculation.
c
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
c     Inputs
c
      integer mxp
      integer npts
      integer npert
      integer nvec
      integer naocc
      integer n_active_chf_prt
      integer active_chf_pert(n_active_chf_prt)
      real*8 wt(npts)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
      real*8 samat(nvec*(nvec+1)/2,npert)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 grho(mxp,2,3)
      real*8 dentol
c
c     Outputs
c
      real*8 ga(naocc,nvec-naocc,npert)
c
c     Workspace
c
      real*8 drhodb(npts)
      real*8 dgrhodb(npts,3)
      real*8 t(npts)
      real*8 dt(npts,3)
      real*8 rhoan(npts,nvec)
      real*8 rhoanb(npts,nvec)
      real*8 dgammadb(npts)
c
c     Local variables
c
      integer nvirt
      integer ipert, iprt, ipt, inu, imu
      integer k, n, n1
      real*8 sval
c
c     Code
c
      nvirt = nvec-naocc
c
      call aclear_dp(rhoan,npts*nvec,0.0d0)
      do k=1,3
         do inu=1,nvec
            do ipt=1,npts
               rhoan(ipt,inu)=rhoan(ipt,inu)
     &                       +grho(ipt,1,k)*amo_grad(ipt,inu,k)
            enddo
         enddo
      enddo
cDEBUG
      do inu=1,nvec
         do ipt=1,npts
            rhoan(ipt,inu)=0.5d0*rhoan(ipt,inu)
         enddo
      enddo
cDEBUG
c
      do iprt=1,n_active_chf_prt
         ipert=active_chf_pert(iprt)
c     do iatm=1,n_lhs_atm
c        do ic=1,3
c           ipert=lhs_atm_list(iatm)+ic
c
c           Construct drhodb and dgrhodb
c
            call aclear_dp(drhodb,npts,0.0d0)
            call aclear_dp(dgrhodb,3*npts,0.0d0)
            n=0
            do inu=1,naocc
               call aclear_dp(t,npts,0.0d0)
               do imu=1,inu
                  n=n+1
                  sval=samat(n,ipert)
                  if (abs(sval).gt.dentol) then
                     do ipt=1,npts
                        t(ipt)=t(ipt)-amo_val(ipt,imu)*sval
                     enddo
                  endif
               enddo
               n1=n
               do imu=inu+1,naocc
                  n1=n1+imu-1
                  sval=samat(n1,ipert)
                  if (abs(sval).gt.dentol) then
                     do ipt=1,npts
                        t(ipt)=t(ipt)-amo_val(ipt,imu)*sval
                     enddo
                  endif
               enddo
               do ipt=1,npts
                  drhodb(ipt)=drhodb(ipt)
     &                       +amo_val(ipt,inu)*t(ipt)
                  t(ipt)=2.0d0*t(ipt)
               enddo
               do k=1,3
                  do ipt=1,npts
                     dgrhodb(ipt,k)=dgrhodb(ipt,k)
     &                             +amo_grad(ipt,inu,k)*t(ipt)
                  enddo
               enddo
            enddo
c           do ipt=1,npts
c              drhodb(ipt)=2.0d0*drhodb(ipt)
c           enddo
c           do k=1,3
c              do ipt=1,npts
c                 dgrhodb(ipt,k)=2.0d0*dgrhodb(ipt,k)
c              enddo
c           enddo
c
c           Construct dgammadb and rhoanb
c
            do ipt=1,npts
c              dgammadb(ipt)=2.0d0*
               dgammadb(ipt)=
     &            (grho(ipt,1,1)*dgrhodb(ipt,1)
     &            +grho(ipt,1,2)*dgrhodb(ipt,2)
     &            +grho(ipt,1,3)*dgrhodb(ipt,3))
            enddo
c
            call aclear_dp(rhoanb,npts*nvec,0.0d0)
            do k=1,3
               do inu=1,nvec
                  do ipt=1,npts
                     rhoanb(ipt,inu)=rhoanb(ipt,inu)
     &                              +dgrhodb(ipt,k)*amo_grad(ipt,inu,k)
                  enddo
               enddo
            enddo
c
c           Construct the final expression
c
            do ipt=1,npts
               dt(ipt,1) = wt(ipt)*(xc_hpt(ipt,irara)*drhodb(ipt)
     &                   +          xc_dhpt(ipt,iragaa)*dgammadb(ipt))
               dt(ipt,2) = wt(ipt)*(xc_dhpt(ipt,iragaa)*drhodb(ipt)*2
     &                   +          xc_dhpt(ipt,igaagaa)*dgammadb(ipt))
               dt(ipt,3) = wt(ipt)*xc_dvpt(ipt,igaa)
            enddo
            do inu=1,nvirt
               do ipt=1,npts
                  t(ipt)=dt(ipt,1)*amo_val(ipt,naocc+inu)
     &                  +dt(ipt,2)*rhoan(ipt,naocc+inu)
     &                  +dt(ipt,3)*rhoanb(ipt,naocc+inu)
               enddo
               do imu=1,naocc
                  do ipt=1,npts
                     ga(imu,inu,ipert)=ga(imu,inu,ipert)
     &                                +t(ipt)*amo_val(ipt,imu)
                  enddo
               enddo
            enddo
            if (.false.) then
c
c              This section uses more FLOPs but therefore memory access
c              is more regular.
c
c              This section led to a cost of:
c                  102.60 secs on Methane (6-31G*).
c                 1143.75 secs on Ethane (6-31G*).
c
               do inu=1,nvirt
                  do ipt=1,npts
                     t(ipt)=dt(ipt,2)*amo_val(ipt,naocc+inu)
                  enddo
                  do imu=1,naocc
                     do ipt=1,npts
                        ga(imu,inu,ipert)=ga(imu,inu,ipert)
     &                                   +t(ipt)*rhoan(ipt,imu)
                     enddo
                  enddo
                  do ipt=1,npts
                     t(ipt)=dt(ipt,3)*amo_val(ipt,naocc+inu)
                  enddo
                  do imu=1,naocc
                     do ipt=1,npts
                        ga(imu,inu,ipert)=ga(imu,inu,ipert)
     &                                   +t(ipt)*rhoanb(ipt,imu)
                     enddo
                  enddo
               enddo
            else
c
c              This section uses less FLOPs but therefore memory access
c              requires more jumps
c
c              This section led to a cost of:
c                   84.05 secs on Methane (6-31G*).
c                  862.62 secs on Ethane (6-31G*).
c
               do imu=1,naocc
                  do ipt=1,npts
                     t(ipt)=dt(ipt,2)*rhoan(ipt,imu)
     &                     +dt(ipt,3)*rhoanb(ipt,imu)
                  enddo
                  do inu=1,nvirt
                     do ipt=1,npts
                        ga(imu,inu,ipert)=ga(imu,inu,ipert)
     &                                   +t(ipt)*amo_val(ipt,naocc+inu)
                     enddo
                  enddo
               enddo
            endif
c        enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine dksm_dft_mo(rkstyp_sw,gradcorr_sw,npts,mxp,nvec,
     &     naocc,nbocc,npert,wt,
     &     grho,drhodb,dgrhodb,dgammadb,
     &     amo_val,amo_grad,
     &     bmo_val,bmo_grad,
     &     xc_dvpt,xc_hpt,xc_dhpt,ga,gb)
      implicit none
c
c     Computes the DFT wavefunction terms that contribute to the 
c     derivative Kohn-Sham matrices.
c
c     For performance analysis we introduce the following variables:
c     N  = the number of atoms
c     Np = a*N = the number of grid points (a = +- thousands)
c     Nb = b*N = the number of basis functions (b = 1 to 2 dozen)
c
c     Parameters
c
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
c     In variables
c
      logical rkstyp_sw ! .true. if closed shell
      logical gradcorr_sw
      integer npts, npert, mxp
      integer nvec
      integer naocc
      integer nbocc
      real*8 wt(mxp)
      real*8 grho(mxp,2,3)
      real*8 drhodb(mxp,2,npert)
      real*8 dgrhodb(mxp,2,3,npert)
      real*8 dgammadb(mxp,3,npert)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 bmo_val(npts,nvec)
      real*8 bmo_grad(npts,nvec,3)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
c     Out variables
c
      real*8 ga(nvec*(nvec+1)/2,npert)
      real*8 gb(nvec*(nvec+1)/2,npert)
c
c     Local variables
c
      integer ipert, ipt, imu, inu, n
      real*8 rhoamn, rhobmn, rhoamnb, rhobmnb

c     real*8 wrr, wrg, wrg2, wgg, wg
c     real*8 wrrdr, wrgdr, wrgdg, wggdg
c     real*8 grgb(nao), dgrgb(nao)
c     real*8 wrrdrn(nao), wrgdrn(nao), wrgdgn(nao), wggdgn(nao)
c     real*8 wrrn(nao), rbn(nao)
c     real*8 wrrp(npts), wrrpd(npts,npert), wrrpdn(npts,nao,npert)
c     real*8 drmn(npts,nao,nao), dgrn(npts,nao,npert), dvn(npts,nao)
c     real*8 ddr(npts,npert), dhnd(npts,nao,npert)
c
c     Code
c
cDEBUG
c     write(*,*)'*** rkstyp_sw   = ',rkstyp_sw
c     write(*,*)'*** gradcorr_sw = ',gradcorr_sw
c     write(*,*)'*** npts        = ',npts       
c     write(*,*)'*** mxp         = ',mxp        
c     write(*,*)'*** nvec        = ',nvec       
c     write(*,*)'*** naocc       = ',naocc      
c     write(*,*)'*** nbocc       = ',nbocc      
c     write(*,*)'*** npert       = ',npert      
c     write(*,*)'*** grho        = ',grho(1,1,1)*0.5d0
c     write(*,*)'*** grho        = ',grho(5,1,3)*0.5d0
c     write(*,*)'*** drhodb      = ',drhodb(1,1,1)
c     write(*,*)'*** drhodb      = ',drhodb(5,1,3)
c     write(*,*)'*** dgrhodb     = ',dgrhodb(1,1,1,1)
c     write(*,*)'*** dgrhodb     = ',dgrhodb(5,1,2,3)
c     write(*,*)'*** dgammadb    = ',dgammadb(1,igaa,1)
c     write(*,*)'*** dgammadb    = ',dgammadb(5,igaa,3)
c     write(*,*)'*** amo_val     = ',amo_val(1,1)
c     write(*,*)'*** amo_val     = ',amo_val(5,8)
c     write(*,*)'*** amo_grad    = ',amo_grad(1,1,1)
c     write(*,*)'*** amo_grad    = ',amo_grad(5,1,3)
c     write(*,*)'*** xc_dvpt     = ',xc_dvpt(1,igaa)
c     write(*,*)'*** xc_hpt      = ',xc_hpt(1,irara)
c     write(*,*)'*** xc_dhpt     = ',xc_dhpt(1,iragaa)
c     write(*,*)'*** xc_dhpt     = ',xc_dhpt(1,igaagaa)
cDEBUG
      if (gradcorr_sw) then
         if (rkstyp_sw) then
c
c           The following code section costs (FLOPs)
c           3*N*Nb*Nb*Na*(13+13+22)
c           = 144*a*b*b*N**4
c
            do ipert = 1, npert
               n=0
               do inu = 1, naocc
                  do imu = 1, inu
                     n=n+1
                     do ipt = 1, npts
 
                        rhoamn = 0.5d0*
     &                   ((grho(ipt,1,1)*amo_grad(ipt,inu,1)+
     &                     grho(ipt,1,2)*amo_grad(ipt,inu,2)+
     &                     grho(ipt,1,3)*amo_grad(ipt,inu,3))*
     &                    amo_val(ipt,imu)+
     &                    (grho(ipt,1,1)*amo_grad(ipt,imu,1)+
     &                     grho(ipt,1,2)*amo_grad(ipt,imu,2)+
     &                     grho(ipt,1,3)*amo_grad(ipt,imu,3))*
     &                    amo_val(ipt,inu))
c
                        rhoamnb = 
     &                    (dgrhodb(ipt,1,1,ipert)*amo_grad(ipt,inu,1)+
     &                     dgrhodb(ipt,1,2,ipert)*amo_grad(ipt,inu,2)+
     &                     dgrhodb(ipt,1,3,ipert)*amo_grad(ipt,inu,3))*
     &                    amo_val(ipt,imu)+
     &                    (dgrhodb(ipt,1,1,ipert)*amo_grad(ipt,imu,1)+
     &                     dgrhodb(ipt,1,2,ipert)*amo_grad(ipt,imu,2)+
     &                     dgrhodb(ipt,1,3,ipert)*amo_grad(ipt,imu,3))*
     &                    amo_val(ipt,inu)
c
                        ga(n,ipert) = ga(n,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irara)*
     &                    amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                    drhodb(ipt,1,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                    dgammadb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*rhoamn*
     &                      drhodb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igaagaa)*rhoamn*
     &                      dgammadb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dvpt(ipt,igaa)*rhoamnb
 
                     enddo
                  enddo
               enddo
            enddo
c
c           The following section of code costs (FLOPs)
c           Nb*Np*3 + Nb*Nb*Np*13 + 3*N*Nb*Np*(5+6) + 3*N*Np*5
c           + 3*N*Nb*Nb*Np*8
c           = 3*a*b*N**2 + 15*a*N**2
c           + 13*a*b*b*N**3 + 33*a*b*N**3
c           + 24*a*b*b*N**4
c           = saving of max 6 (realized on methane 2)
c
c           do inu = 1, nao
c              do ipt = 1, npts
c                 dvn(ipt,inu) = wt(ipt)*xc_dvpt(ipt,igaa)*
c    &                           bfn_val(ipt,inu)
c              enddo
c           enddo
c           do inu = 1, nao
c              do imu = 1, nao
c                 do ipt = 1, npts
c                    drmn(ipt,imu,inu) = 
c    &                   (grho(ipt,1,1)*bfng_val(ipt,inu,1)+
c    &                    grho(ipt,1,2)*bfng_val(ipt,inu,2)+
c    &                    grho(ipt,1,3)*bfng_val(ipt,inu,3))*
c    &                   bfn_val(ipt,imu)+
c    &                   (grho(ipt,1,1)*bfng_val(ipt,imu,1)+
c    &                    grho(ipt,1,2)*bfng_val(ipt,imu,2)+
c    &                    grho(ipt,1,3)*bfng_val(ipt,imu,3))*
c    &                   bfn_val(ipt,inu)
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do ipt = 1, npts
c                    dgrn(ipt,inu,ipert) = 
c    &                   (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,inu,1)+
c    &                    dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,inu,2)+
c    &                    dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,inu,3))
c                    dhnd(ipt,inu,ipert) = wt(ipt)*(
c    &                    xc_hpt(ipt,irara)*bfn_val(ipt,inu)*
c    &                    drhodb(ipt,1,ipert)
c    &                  + xc_dhpt(ipt,iragaa)*bfn_val(ipt,inu)*
c    &                    dgammadb(ipt,1,ipert))
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do ipt = 1, npts
c                 ddr(ipt,ipert) = wt(ipt)*(
c    &                2*xc_dhpt(ipt,iragaa)*drhodb(ipt,1,ipert)
c    &              + xc_dhpt(ipt,igaagaa)*dgammadb(ipt,1,ipert))
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                  + dhnd(ipt,inu,ipert)*bfn_val(ipt,imu)
c    &                  + ddr(ipt,ipert)*drmn(ipt,imu,inu)
c    &                  + dvn(ipt,imu)*dgrn(ipt,inu,ipert)
c    &                  + dvn(ipt,inu)*dgrn(ipt,imu,ipert)
c                    enddo
c                 enddo
c              enddo
c           enddo
         else ! rkstyp_sw
            do ipert = 1, npert
               n=0
               do inu = 1, naocc
                  do imu = 1, inu
                     n=n+1
                     do ipt = 1, npts
c
                        rhoamn = 
     &                    (grho(ipt,1,1)*amo_grad(ipt,inu,1)+
     &                     grho(ipt,1,2)*amo_grad(ipt,inu,2)+
     &                     grho(ipt,1,3)*amo_grad(ipt,inu,3))*
     &                    amo_val(ipt,imu)+
     &                    (grho(ipt,1,1)*amo_grad(ipt,imu,1)+
     &                     grho(ipt,1,2)*amo_grad(ipt,imu,2)+
     &                     grho(ipt,1,3)*amo_grad(ipt,imu,3))*
     &                    amo_val(ipt,inu)
c
                        rhobmn = 
     &                    (grho(ipt,2,1)*amo_grad(ipt,inu,1)+
     &                     grho(ipt,2,2)*amo_grad(ipt,inu,2)+
     &                     grho(ipt,2,3)*amo_grad(ipt,inu,3))*
     &                    amo_val(ipt,imu)+
     &                    (grho(ipt,2,1)*amo_grad(ipt,imu,1)+
     &                     grho(ipt,2,2)*amo_grad(ipt,imu,2)+
     &                     grho(ipt,2,3)*amo_grad(ipt,imu,3))*
     &                    amo_val(ipt,inu)
c
                        rhoamnb = 
     &                    (dgrhodb(ipt,1,1,ipert)*amo_grad(ipt,inu,1)+
     &                     dgrhodb(ipt,1,2,ipert)*amo_grad(ipt,inu,2)+
     &                     dgrhodb(ipt,1,3,ipert)*amo_grad(ipt,inu,3))*
     &                    amo_val(ipt,imu)+
     &                    (dgrhodb(ipt,1,1,ipert)*amo_grad(ipt,imu,1)+
     &                     dgrhodb(ipt,1,2,ipert)*amo_grad(ipt,imu,2)+
     &                     dgrhodb(ipt,1,3,ipert)*amo_grad(ipt,imu,3))*
     &                    amo_val(ipt,inu)
c
                        rhobmnb = 
     &                    (dgrhodb(ipt,2,1,ipert)*amo_grad(ipt,inu,1)+
     &                     dgrhodb(ipt,2,2,ipert)*amo_grad(ipt,inu,2)+
     &                     dgrhodb(ipt,2,3,ipert)*amo_grad(ipt,inu,3))*
     &                    amo_val(ipt,imu)+
     &                    (dgrhodb(ipt,2,1,ipert)*amo_grad(ipt,imu,1)+
     &                     dgrhodb(ipt,2,2,ipert)*amo_grad(ipt,imu,2)+
     &                     dgrhodb(ipt,2,3,ipert)*amo_grad(ipt,imu,3))*
     &                    amo_val(ipt,inu)
c
                        ga(n,ipert) = ga(n,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irara)*
     &                    amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                    drhodb(ipt,1,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irarb)*
     &                    amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                    drhodb(ipt,2,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                    amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                    dgammadb(ipt,1,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,iragab)*
     &                    amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                    dgammadb(ipt,2,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,iragbb)*
     &                    amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                    dgammadb(ipt,3,ipert)
c
     &                  + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*rhoamn*
     &                      drhodb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,iragab)*rhobmn*
     &                      drhodb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,irbgaa)*rhoamn*
     &                      drhodb(ipt,2,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,irbgab)*rhobmn*
     &                      drhodb(ipt,2,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igaagaa)*rhoamn*
     &                      dgammadb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igaagab)*rhobmn*
     &                      dgammadb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igaagab)*rhoamn*
     &                      dgammadb(ipt,2,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igabgab)*rhobmn*
     &                      dgammadb(ipt,2,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igaagbb)*rhoamn*
     &                      dgammadb(ipt,3,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igabgbb)*rhobmn*
     &                      dgammadb(ipt,3,ipert)
c
     &                  + 2*wt(ipt)*xc_dvpt(ipt,igaa)*rhoamnb
     &                  +   wt(ipt)*xc_dvpt(ipt,igab)*rhobmnb
c
                     enddo
                  enddo
               enddo
            enddo
            do ipert = 1, npert
               n=0
               do inu = 1, nbocc
                  do imu = 1, inu
                     n=n+1
                     do ipt = 1, npts
c
                        rhoamn = 
     &                    (grho(ipt,1,1)*bmo_grad(ipt,inu,1)+
     &                     grho(ipt,1,2)*bmo_grad(ipt,inu,2)+
     &                     grho(ipt,1,3)*bmo_grad(ipt,inu,3))*
     &                    bmo_val(ipt,imu)+
     &                    (grho(ipt,1,1)*bmo_grad(ipt,imu,1)+
     &                     grho(ipt,1,2)*bmo_grad(ipt,imu,2)+
     &                     grho(ipt,1,3)*bmo_grad(ipt,imu,3))*
     &                    bmo_val(ipt,inu)
c
                        rhobmn = 
     &                    (grho(ipt,2,1)*bmo_grad(ipt,inu,1)+
     &                     grho(ipt,2,2)*bmo_grad(ipt,inu,2)+
     &                     grho(ipt,2,3)*bmo_grad(ipt,inu,3))*
     &                    bmo_val(ipt,imu)+
     &                    (grho(ipt,2,1)*bmo_grad(ipt,imu,1)+
     &                     grho(ipt,2,2)*bmo_grad(ipt,imu,2)+
     &                     grho(ipt,2,3)*bmo_grad(ipt,imu,3))*
     &                    bmo_val(ipt,inu)
c
                        rhoamnb = 
     &                    (dgrhodb(ipt,1,1,ipert)*bmo_grad(ipt,inu,1)+
     &                     dgrhodb(ipt,1,2,ipert)*bmo_grad(ipt,inu,2)+
     &                     dgrhodb(ipt,1,3,ipert)*bmo_grad(ipt,inu,3))*
     &                    bmo_val(ipt,imu)+
     &                    (dgrhodb(ipt,1,1,ipert)*bmo_grad(ipt,imu,1)+
     &                     dgrhodb(ipt,1,2,ipert)*bmo_grad(ipt,imu,2)+
     &                     dgrhodb(ipt,1,3,ipert)*bmo_grad(ipt,imu,3))*
     &                    bmo_val(ipt,inu)
c
                        rhobmnb = 
     &                    (dgrhodb(ipt,2,1,ipert)*bmo_grad(ipt,inu,1)+
     &                     dgrhodb(ipt,2,2,ipert)*bmo_grad(ipt,inu,2)+
     &                     dgrhodb(ipt,2,3,ipert)*bmo_grad(ipt,inu,3))*
     &                    bmo_val(ipt,imu)+
     &                    (dgrhodb(ipt,2,1,ipert)*bmo_grad(ipt,imu,1)+
     &                     dgrhodb(ipt,2,2,ipert)*bmo_grad(ipt,imu,2)+
     &                     dgrhodb(ipt,2,3,ipert)*bmo_grad(ipt,imu,3))*
     &                    bmo_val(ipt,inu)
c
                        gb(n,ipert) = gb(n,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irarb)*
     &                    bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                    drhodb(ipt,1,ipert)
     &                  + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                    bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                    drhodb(ipt,2,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgaa)*
     &                    bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                    dgammadb(ipt,1,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgab)*
     &                    bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                    dgammadb(ipt,2,ipert)
     &                  + wt(ipt)*xc_dhpt(ipt,irbgbb)*
     &                    bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                    dgammadb(ipt,3,ipert)
c
     &                  + 2*wt(ipt)*xc_dhpt(ipt,iragbb)*rhobmn*
     &                      drhodb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,iragab)*rhoamn*
     &                      drhodb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,irbgbb)*rhobmn*
     &                      drhodb(ipt,2,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,irbgab)*rhoamn*
     &                      drhodb(ipt,2,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igaagbb)*rhobmn*
     &                      dgammadb(ipt,1,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igaagab)*rhoamn*
     &                      dgammadb(ipt,1,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igabgbb)*rhobmn*
     &                      dgammadb(ipt,2,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igabgab)*rhoamn*
     &                      dgammadb(ipt,2,ipert)
     &                  + 2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*rhobmn*
     &                      dgammadb(ipt,3,ipert)
     &                  +   wt(ipt)*xc_dhpt(ipt,igabgbb)*rhoamn*
     &                      dgammadb(ipt,3,ipert)
c
     &                  + 2*wt(ipt)*xc_dvpt(ipt,igbb)*rhobmnb
     &                  +   wt(ipt)*xc_dvpt(ipt,igab)*rhoamnb
c
                     enddo
                  enddo
               enddo
            enddo
         endif ! rkstyp_sw
      else ! gradcorr_sw
         if (rkstyp_sw) then
c           if (.true.) then
c
c           The following piece of code costs (FLOPs)
c           Np + Np*Nb + Np*3*N*Nb + Np*3*N*Nb*Nb*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + a*b*N**2 + a*N
c           = saving max 2.5 (realized on methane: 4 s out of 127 s)
c
c           do ipt = 1, npts
c              wrr = wt(ipt)*xc_hpt(ipt,irara)
c              do inu = 1, nao
c                 wrrn(inu) = wrr*bfn_val(ipt,inu)
c              enddo
c              do ipert = 1, npert
c                 do inu = 1, nao
c                    rbn(inu) = bfn_val(ipt,inu)*drhodb(ipt,1,ipert)
c                 enddo
c                 do inu = 1, nao
c                    do imu = 1, nao
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrn(inu)*rbn(imu)
c                    enddo
c                 enddo
c              enddo
c           enddo
c
c           The following section of code costs (FLOPs)
c           Np + 3*N*Np + 3*N*Nb*Np + 3*N*Nb*Nb*Np*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + 3*a*N**2 + a*N
c           = saving max 2.5 (realized on methane: 1.5)
c
c           do ipt = 1, npts
c              wrrp(ipt) = wt(ipt)*xc_hpt(ipt,irara)
c           enddo
c           do ipert = 1, npert
c              do ipt = 1, npts
c                 wrrpd(ipt,ipert) = wrrp(ipt)*drhodb(ipt,1,ipert)
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do ipt = 1, npts
c                    wrrpdn(ipt,inu,ipert) = wrrpd(ipt,ipert)*
c    &                                       bfn_val(ipt,inu)
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrpdn(ipt,inu,ipert)*bfn_val(ipt,imu)
c                    enddo
c                 enddo
c              enddo
c           enddo
c           else
c
c           The following piece of code costs (FLOPs)
c           3*N*Nb*Nb*Np*5 = 15*a*b*b*N**4
c
            do ipert = 1, npert
               n=0
               do inu = 1, naocc
                  do imu = 1, inu
                     n=n+1
                     do ipt = 1, npts
                        ga(n,ipert) = ga(n,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irara)*
     &                      amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                      drhodb(ipt,1,ipert)
                     enddo
                  enddo
               enddo
            enddo
c
c           The cost of the following bit of code is (FLOPs)
c           Np + Np*Nb + Np*Nb*3*N + 3*N*Nb*Nb*Np*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + a*b*N**2 + a*N
c
c           do ipt = 1, npts
c              wrr = wt(ipt)*xc_hpt(ipt,irara)
c              do inu = 1, nao
c                 wrrn(ipt,inu) = wrr*bfn_val(ipt,inu)
c                 do ipert = 1, npert
c                    drn(ipt,inu,ipert) = bfn_val(ipt,inu)*
c    &                                    drhodb(ipt,1,ipert)
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrn(ipt,inu)*drn(ipt,imu,ipert)
c                    enddo
c                 enddo
c              enddo
c           enddo
c           endif
         else ! rkstyp_sw
            do ipert = 1, npert
               n=0
               do inu = 1, naocc
                  do imu = 1, inu
                     n=n+1
                     do ipt = 1, npts
                        ga(n,ipert) = ga(n,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irara)*
     &                      amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                      drhodb(ipt,1,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irarb)*
     &                      amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                      drhodb(ipt,2,ipert)
                     enddo
                  enddo
               enddo
            enddo
            do ipert = 1, npert
               n=0
               do inu = 1, nbocc
                  do imu = 1, nbocc
                     n=n+1
                     do ipt = 1, npts
                        gb(n,ipert) = gb(n,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                      bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                      drhodb(ipt,2,ipert)
     &                    + wt(ipt)*xc_hpt(ipt,irarb)*
     &                      bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                      drhodb(ipt,1,ipert)
                     enddo
                  enddo
               enddo
            enddo
         endif ! rkstyp_sw
      endif ! gradcorr_sw
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dksm_dft_mo_scr(rkstyp_sw,gradcorr_sw,npts,mxp,nvec,
     &     naocc,nbocc,npert,
     &     active_chf_pert,n_active_chf_prt,
     &     wt,grho,drhodb,dgrhodb,dgammadb,
     &     amo_val,amo_grad,
     &     bmo_val,bmo_grad,
     &     xc_dvpt,xc_hpt,xc_dhpt,ga,gb)
      implicit none
c
c     Computes the DFT wavefunction terms that contribute to the 
c     derivative Kohn-Sham matrices.
c
c     For performance analysis we introduce the following variables:
c     N  = the number of atoms
c     Np = a*N = the number of grid points (a = +- thousands)
c     Nb = b*N = the number of basis functions (b = 1 to 2 dozen)
c
c     Parameters
c
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
c     In variables
c
      logical rkstyp_sw ! .true. if closed shell
      logical gradcorr_sw
      integer npts, npert, mxp
      integer nvec
      integer naocc
      integer nbocc
      integer n_active_chf_prt
      integer active_chf_pert(n_active_chf_prt)
      real*8 wt(mxp)
      real*8 grho(mxp,2,3)
      real*8 drhodb(mxp,2,npert)
      real*8 dgrhodb(mxp,2,3,npert)
      real*8 dgammadb(mxp,3,npert)
      real*8 amo_val(npts,nvec)
      real*8 amo_grad(npts,nvec,3)
      real*8 bmo_val(npts,nvec)
      real*8 bmo_grad(npts,nvec,3)
      real*8 xc_dvpt(mxp,3)
      real*8 xc_hpt(mxp,3)
      real*8 xc_dhpt(mxp,12)
c
c     Out variables
c
      real*8 ga(nvec*(nvec+1)/2,npert)
      real*8 gb(nvec*(nvec+1)/2,npert)
c
c     Local variables
c
      integer ipert, iprt, ipt, imu, inu, n
c     integer iatm, ic
      real*8 rhoamn, rhobmn, rhoamnb, rhobmnb

c     real*8 wrr, wrg, wrg2, wgg, wg
c     real*8 wrrdr, wrgdr, wrgdg, wggdg
c     real*8 grgb(nao), dgrgb(nao)
c     real*8 wrrdrn(nao), wrgdrn(nao), wrgdgn(nao), wggdgn(nao)
c     real*8 wrrn(nao), rbn(nao)
c     real*8 wrrp(npts), wrrpd(npts,npert), wrrpdn(npts,nao,npert)
c     real*8 drmn(npts,nao,nao), dgrn(npts,nao,npert), dvn(npts,nao)
c     real*8 ddr(npts,npert), dhnd(npts,nao,npert)
c
c     Code
c
      if (gradcorr_sw) then
         if (rkstyp_sw) then
c
c           The following code section costs (FLOPs)
c           3*N*Nb*Nb*Na*(13+13+22)
c           = 144*a*b*b*N**4
c
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do ic=1,3
c                 ipert=lhs_atm_list(iatm)+ic
c                 iprt=3*(iatm-1)+ic
                  n=0
                  do inu = 1, naocc
                     do imu = 1, inu
                        n=n+1
                        do ipt = 1, npts
 
                           rhoamn = 0.5d0*
     &                      ((grho(ipt,1,1)*amo_grad(ipt,inu,1)+
     &                        grho(ipt,1,2)*amo_grad(ipt,inu,2)+
     &                        grho(ipt,1,3)*amo_grad(ipt,inu,3))*
     &                       amo_val(ipt,imu)+
     &                       (grho(ipt,1,1)*amo_grad(ipt,imu,1)+
     &                        grho(ipt,1,2)*amo_grad(ipt,imu,2)+
     &                        grho(ipt,1,3)*amo_grad(ipt,imu,3))*
     &                       amo_val(ipt,inu))
c
                           rhoamnb = 
     &                       (dgrhodb(ipt,1,1,iprt)*
     &                        amo_grad(ipt,inu,1)+
     &                        dgrhodb(ipt,1,2,iprt)*
     &                        amo_grad(ipt,inu,2)+
     &                        dgrhodb(ipt,1,3,iprt)*
     &                        amo_grad(ipt,inu,3))*
     &                       amo_val(ipt,imu)+
     &                       (dgrhodb(ipt,1,1,iprt)*
     &                        amo_grad(ipt,imu,1)+
     &                        dgrhodb(ipt,1,2,iprt)*
     &                        amo_grad(ipt,imu,2)+
     &                        dgrhodb(ipt,1,3,iprt)*
     &                        amo_grad(ipt,imu,3))*
     &                       amo_val(ipt,inu)
c
                           ga(n,ipert) = ga(n,ipert)
     &                     + wt(ipt)*xc_hpt(ipt,irara)*
     &                       amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                       drhodb(ipt,1,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                       amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                       dgammadb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*rhoamn*
     &                         drhodb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igaagaa)*rhoamn*
     &                         dgammadb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dvpt(ipt,igaa)*rhoamnb
 
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
c
c           The following section of code costs (FLOPs)
c           Nb*Np*3 + Nb*Nb*Np*13 + 3*N*Nb*Np*(5+6) + 3*N*Np*5
c           + 3*N*Nb*Nb*Np*8
c           = 3*a*b*N**2 + 15*a*N**2
c           + 13*a*b*b*N**3 + 33*a*b*N**3
c           + 24*a*b*b*N**4
c           = saving of max 6 (realized on methane 2)
c
c           do inu = 1, nao
c              do ipt = 1, npts
c                 dvn(ipt,inu) = wt(ipt)*xc_dvpt(ipt,igaa)*
c    &                           bfn_val(ipt,inu)
c              enddo
c           enddo
c           do inu = 1, nao
c              do imu = 1, nao
c                 do ipt = 1, npts
c                    drmn(ipt,imu,inu) = 
c    &                   (grho(ipt,1,1)*bfng_val(ipt,inu,1)+
c    &                    grho(ipt,1,2)*bfng_val(ipt,inu,2)+
c    &                    grho(ipt,1,3)*bfng_val(ipt,inu,3))*
c    &                   bfn_val(ipt,imu)+
c    &                   (grho(ipt,1,1)*bfng_val(ipt,imu,1)+
c    &                    grho(ipt,1,2)*bfng_val(ipt,imu,2)+
c    &                    grho(ipt,1,3)*bfng_val(ipt,imu,3))*
c    &                   bfn_val(ipt,inu)
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do ipt = 1, npts
c                    dgrn(ipt,inu,ipert) = 
c    &                   (dgrhodb(ipt,1,1,ipert)*bfng_val(ipt,inu,1)+
c    &                    dgrhodb(ipt,1,2,ipert)*bfng_val(ipt,inu,2)+
c    &                    dgrhodb(ipt,1,3,ipert)*bfng_val(ipt,inu,3))
c                    dhnd(ipt,inu,ipert) = wt(ipt)*(
c    &                    xc_hpt(ipt,irara)*bfn_val(ipt,inu)*
c    &                    drhodb(ipt,1,ipert)
c    &                  + xc_dhpt(ipt,iragaa)*bfn_val(ipt,inu)*
c    &                    dgammadb(ipt,1,ipert))
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do ipt = 1, npts
c                 ddr(ipt,ipert) = wt(ipt)*(
c    &                2*xc_dhpt(ipt,iragaa)*drhodb(ipt,1,ipert)
c    &              + xc_dhpt(ipt,igaagaa)*dgammadb(ipt,1,ipert))
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                  + dhnd(ipt,inu,ipert)*bfn_val(ipt,imu)
c    &                  + ddr(ipt,ipert)*drmn(ipt,imu,inu)
c    &                  + dvn(ipt,imu)*dgrn(ipt,inu,ipert)
c    &                  + dvn(ipt,inu)*dgrn(ipt,imu,ipert)
c                    enddo
c                 enddo
c              enddo
c           enddo
         else ! rkstyp_sw
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do ic=1,3
c                 ipert=lhs_atm_list(iatm)+ic
c                 iprt=3*(iatm-1)+ic
                  n=0
                  do inu = 1, naocc
                     do imu = 1, inu
                        n=n+1
                        do ipt = 1, npts
c
                           rhoamn = 
     &                       (grho(ipt,1,1)*amo_grad(ipt,inu,1)+
     &                        grho(ipt,1,2)*amo_grad(ipt,inu,2)+
     &                        grho(ipt,1,3)*amo_grad(ipt,inu,3))*
     &                       amo_val(ipt,imu)+
     &                       (grho(ipt,1,1)*amo_grad(ipt,imu,1)+
     &                        grho(ipt,1,2)*amo_grad(ipt,imu,2)+
     &                        grho(ipt,1,3)*amo_grad(ipt,imu,3))*
     &                       amo_val(ipt,inu)
c
                           rhobmn = 
     &                       (grho(ipt,2,1)*amo_grad(ipt,inu,1)+
     &                        grho(ipt,2,2)*amo_grad(ipt,inu,2)+
     &                        grho(ipt,2,3)*amo_grad(ipt,inu,3))*
     &                       amo_val(ipt,imu)+
     &                       (grho(ipt,2,1)*amo_grad(ipt,imu,1)+
     &                        grho(ipt,2,2)*amo_grad(ipt,imu,2)+
     &                        grho(ipt,2,3)*amo_grad(ipt,imu,3))*
     &                       amo_val(ipt,inu)
c
                           rhoamnb = 
     &                       (dgrhodb(ipt,1,1,iprt)*
     &                        amo_grad(ipt,inu,1)+
     &                        dgrhodb(ipt,1,2,iprt)*
     &                        amo_grad(ipt,inu,2)+
     &                        dgrhodb(ipt,1,3,iprt)*
     &                        amo_grad(ipt,inu,3))*
     &                       amo_val(ipt,imu)+
     &                       (dgrhodb(ipt,1,1,iprt)*
     &                        amo_grad(ipt,imu,1)+
     &                        dgrhodb(ipt,1,2,iprt)*
     &                        amo_grad(ipt,imu,2)+
     &                        dgrhodb(ipt,1,3,iprt)*
     &                        amo_grad(ipt,imu,3))*
     &                       amo_val(ipt,inu)
c
                           rhobmnb = 
     &                       (dgrhodb(ipt,2,1,iprt)*
     &                        amo_grad(ipt,inu,1)+
     &                        dgrhodb(ipt,2,2,iprt)*
     &                        amo_grad(ipt,inu,2)+
     &                        dgrhodb(ipt,2,3,iprt)*
     &                        amo_grad(ipt,inu,3))*
     &                       amo_val(ipt,imu)+
     &                       (dgrhodb(ipt,2,1,iprt)*
     &                        amo_grad(ipt,imu,1)+
     &                        dgrhodb(ipt,2,2,iprt)*
     &                        amo_grad(ipt,imu,2)+
     &                        dgrhodb(ipt,2,3,iprt)*
     &                        amo_grad(ipt,imu,3))*
     &                       amo_val(ipt,inu)
c
                           ga(n,ipert) = ga(n,ipert)
     &                     + wt(ipt)*xc_hpt(ipt,irara)*
     &                       amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                       drhodb(ipt,1,iprt)
     &                     + wt(ipt)*xc_hpt(ipt,irarb)*
     &                       amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                       drhodb(ipt,2,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,iragaa)*
     &                       amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                       dgammadb(ipt,1,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,iragab)*
     &                       amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                       dgammadb(ipt,2,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,iragbb)*
     &                       amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                       dgammadb(ipt,3,iprt)
c
     &                     + 2*wt(ipt)*xc_dhpt(ipt,iragaa)*rhoamn*
     &                         drhodb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,iragab)*rhobmn*
     &                         drhodb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,irbgaa)*rhoamn*
     &                         drhodb(ipt,2,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,irbgab)*rhobmn*
     &                         drhodb(ipt,2,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igaagaa)*rhoamn*
     &                         dgammadb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igaagab)*rhobmn*
     &                         dgammadb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igaagab)*rhoamn*
     &                         dgammadb(ipt,2,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igabgab)*rhobmn*
     &                         dgammadb(ipt,2,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igaagbb)*rhoamn*
     &                         dgammadb(ipt,3,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igabgbb)*rhobmn*
     &                         dgammadb(ipt,3,iprt)
c
     &                     + 2*wt(ipt)*xc_dvpt(ipt,igaa)*rhoamnb
     &                     +   wt(ipt)*xc_dvpt(ipt,igab)*rhobmnb
c
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do ic=1,3
c                 ipert=lhs_atm_list(iatm)+ic
c                 iprt=3*(iatm-1)+ic
                  n=0
                  do inu = 1, nbocc
                     do imu = 1, inu
                        n=n+1
                        do ipt = 1, npts
c
                           rhoamn = 
     &                       (grho(ipt,1,1)*bmo_grad(ipt,inu,1)+
     &                        grho(ipt,1,2)*bmo_grad(ipt,inu,2)+
     &                        grho(ipt,1,3)*bmo_grad(ipt,inu,3))*
     &                       bmo_val(ipt,imu)+
     &                       (grho(ipt,1,1)*bmo_grad(ipt,imu,1)+
     &                        grho(ipt,1,2)*bmo_grad(ipt,imu,2)+
     &                        grho(ipt,1,3)*bmo_grad(ipt,imu,3))*
     &                       bmo_val(ipt,inu)
c
                           rhobmn = 
     &                       (grho(ipt,2,1)*bmo_grad(ipt,inu,1)+
     &                        grho(ipt,2,2)*bmo_grad(ipt,inu,2)+
     &                        grho(ipt,2,3)*bmo_grad(ipt,inu,3))*
     &                       bmo_val(ipt,imu)+
     &                       (grho(ipt,2,1)*bmo_grad(ipt,imu,1)+
     &                        grho(ipt,2,2)*bmo_grad(ipt,imu,2)+
     &                        grho(ipt,2,3)*bmo_grad(ipt,imu,3))*
     &                       bmo_val(ipt,inu)
c
                           rhoamnb = 
     &                       (dgrhodb(ipt,1,1,iprt)*
     &                        bmo_grad(ipt,inu,1)+
     &                        dgrhodb(ipt,1,2,iprt)*
     &                        bmo_grad(ipt,inu,2)+
     &                        dgrhodb(ipt,1,3,iprt)*
     &                        bmo_grad(ipt,inu,3))*
     &                       bmo_val(ipt,imu)+
     &                       (dgrhodb(ipt,1,1,iprt)*
     &                        bmo_grad(ipt,imu,1)+
     &                        dgrhodb(ipt,1,2,iprt)*
     &                        bmo_grad(ipt,imu,2)+
     &                        dgrhodb(ipt,1,3,iprt)*
     &                        bmo_grad(ipt,imu,3))*
     &                       bmo_val(ipt,inu)
c
                           rhobmnb = 
     &                       (dgrhodb(ipt,2,1,iprt)*
     &                        bmo_grad(ipt,inu,1)+
     &                        dgrhodb(ipt,2,2,iprt)*
     &                        bmo_grad(ipt,inu,2)+
     &                        dgrhodb(ipt,2,3,iprt)*
     &                        bmo_grad(ipt,inu,3))*
     &                       bmo_val(ipt,imu)+
     &                       (dgrhodb(ipt,2,1,iprt)*
     &                        bmo_grad(ipt,imu,1)+
     &                        dgrhodb(ipt,2,2,iprt)*
     &                        bmo_grad(ipt,imu,2)+
     &                        dgrhodb(ipt,2,3,iprt)*
     &                        bmo_grad(ipt,imu,3))*
     &                       bmo_val(ipt,inu)
c
                           gb(n,ipert) = gb(n,ipert)
     &                     + wt(ipt)*xc_hpt(ipt,irarb)*
     &                       bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                       drhodb(ipt,1,iprt)
     &                     + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                       bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                       drhodb(ipt,2,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgaa)*
     &                       bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                       dgammadb(ipt,1,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgab)*
     &                       bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                       dgammadb(ipt,2,iprt)
     &                     + wt(ipt)*xc_dhpt(ipt,irbgbb)*
     &                       bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                       dgammadb(ipt,3,iprt)
c
     &                     + 2*wt(ipt)*xc_dhpt(ipt,iragbb)*rhobmn*
     &                         drhodb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,iragab)*rhoamn*
     &                         drhodb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,irbgbb)*rhobmn*
     &                         drhodb(ipt,2,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,irbgab)*rhoamn*
     &                         drhodb(ipt,2,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igaagbb)*rhobmn*
     &                         dgammadb(ipt,1,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igaagab)*rhoamn*
     &                         dgammadb(ipt,1,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igabgbb)*rhobmn*
     &                         dgammadb(ipt,2,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igabgab)*rhoamn*
     &                         dgammadb(ipt,2,iprt)
     &                     + 2*wt(ipt)*xc_dhpt(ipt,igbbgbb)*rhobmn*
     &                         dgammadb(ipt,3,iprt)
     &                     +   wt(ipt)*xc_dhpt(ipt,igabgbb)*rhoamn*
     &                         dgammadb(ipt,3,iprt)
c
     &                     + 2*wt(ipt)*xc_dvpt(ipt,igbb)*rhobmnb
     &                     +   wt(ipt)*xc_dvpt(ipt,igab)*rhoamnb
c
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
         endif ! rkstyp_sw
      else ! gradcorr_sw
         if (rkstyp_sw) then
c           if (.true.) then
c
c           The following piece of code costs (FLOPs)
c           Np + Np*Nb + Np*3*N*Nb + Np*3*N*Nb*Nb*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + a*b*N**2 + a*N
c           = saving max 2.5 (realized on methane: 4 s out of 127 s)
c
c           do ipt = 1, npts
c              wrr = wt(ipt)*xc_hpt(ipt,irara)
c              do inu = 1, nao
c                 wrrn(inu) = wrr*bfn_val(ipt,inu)
c              enddo
c              do ipert = 1, npert
c                 do inu = 1, nao
c                    rbn(inu) = bfn_val(ipt,inu)*drhodb(ipt,1,ipert)
c                 enddo
c                 do inu = 1, nao
c                    do imu = 1, nao
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrn(inu)*rbn(imu)
c                    enddo
c                 enddo
c              enddo
c           enddo
c
c           The following section of code costs (FLOPs)
c           Np + 3*N*Np + 3*N*Nb*Np + 3*N*Nb*Nb*Np*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + 3*a*N**2 + a*N
c           = saving max 2.5 (realized on methane: 1.5)
c
c           do ipt = 1, npts
c              wrrp(ipt) = wt(ipt)*xc_hpt(ipt,irara)
c           enddo
c           do ipert = 1, npert
c              do ipt = 1, npts
c                 wrrpd(ipt,ipert) = wrrp(ipt)*drhodb(ipt,1,ipert)
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do ipt = 1, npts
c                    wrrpdn(ipt,inu,ipert) = wrrpd(ipt,ipert)*
c    &                                       bfn_val(ipt,inu)
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrpdn(ipt,inu,ipert)*bfn_val(ipt,imu)
c                    enddo
c                 enddo
c              enddo
c           enddo
c           else
c
c           The following piece of code costs (FLOPs)
c           3*N*Nb*Nb*Np*5 = 15*a*b*b*N**4
c
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do ic=1,3
c                 ipert=lhs_atm_list(iatm)+ic
c                 iprt=3*(iatm-1)+ic
                  n=0
                  do inu = 1, naocc
                     do imu = 1, inu
                        n=n+1
                        do ipt = 1, npts
                           ga(n,ipert) = ga(n,ipert)
     &                       + wt(ipt)*xc_hpt(ipt,irara)*
     &                         amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                         drhodb(ipt,1,iprt)
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
c
c           The cost of the following bit of code is (FLOPs)
c           Np + Np*Nb + Np*Nb*3*N + 3*N*Nb*Nb*Np*2
c           = 6*a*b*b*N**4 + 3*a*b*N**3 + a*b*N**2 + a*N
c
c           do ipt = 1, npts
c              wrr = wt(ipt)*xc_hpt(ipt,irara)
c              do inu = 1, nao
c                 wrrn(ipt,inu) = wrr*bfn_val(ipt,inu)
c                 do ipert = 1, npert
c                    drn(ipt,inu,ipert) = bfn_val(ipt,inu)*
c    &                                    drhodb(ipt,1,ipert)
c                 enddo
c              enddo
c           enddo
c           do ipert = 1, npert
c              do inu = 1, nao
c                 do imu = 1, nao
c                    do ipt = 1, npts
c                       ga(imu,inu,ipert) = ga(imu,inu,ipert)
c    &                    + wrrn(ipt,inu)*drn(ipt,imu,ipert)
c                    enddo
c                 enddo
c              enddo
c           enddo
c           endif
         else ! rkstyp_sw
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do ic=1,3
c                 ipert=lhs_atm_list(iatm)+ic
c                 iprt=3*(iatm-1)+ic
                  n=0
                  do inu = 1, naocc
                     do imu = 1, inu
                        n=n+1
                        do ipt = 1, npts
                           ga(n,ipert) = ga(n,ipert)
     &                       + wt(ipt)*xc_hpt(ipt,irara)*
     &                         amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                         drhodb(ipt,1,iprt)
     &                       + wt(ipt)*xc_hpt(ipt,irarb)*
     &                         amo_val(ipt,imu)*amo_val(ipt,inu)*
     &                         drhodb(ipt,2,iprt)
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
            do iprt=1,n_active_chf_prt
               ipert=active_chf_pert(iprt)
c           do iatm=1,n_lhs_atm
c              do ic=1,3
c                 ipert=lhs_atm_list(iatm)+ic
c                 iprt=3*(iatm-1)+ic
                  n=0
                  do inu = 1, nbocc
                     do imu = 1, nbocc
                        n=n+1
                        do ipt = 1, npts
                           gb(n,ipert) = gb(n,ipert)
     &                       + wt(ipt)*xc_hpt(ipt,irbrb)*
     &                         bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                         drhodb(ipt,2,iprt)
     &                       + wt(ipt)*xc_hpt(ipt,irarb)*
     &                         bmo_val(ipt,imu)*bmo_val(ipt,inu)*
     &                         drhodb(ipt,1,iprt)
                        enddo
                     enddo
                  enddo
c              enddo
            enddo
         endif ! rkstyp_sw
      endif ! gradcorr_sw
c
      return
      end
