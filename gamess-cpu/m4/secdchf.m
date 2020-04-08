c
c gdf:  this version all perturbations treated simultaneously
c
      subroutine chfdrv2(eps,lstop,skipp)
      implicit REAL  (a-h,o-z)
      logical lstop,skipp
      dimension skipp(*), eps(*)
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/global)
INCLUDE(common/vcore)
INCLUDE(common/gmempara)
      logical pg_create_inf, opg_root
      dimension iatms(3*maxat)
      character *9 fnm
      character *7 snm
      data fnm/"secdchf.m"/
      data snm/"chfdrv2"/

      npstar = 0
      npfin = np
      nov = nocca*nvirta
      nvv = nvirta*nvirta
      maxc = 50
      mynode = ipg_nodeid()

c  find number of perturbations
      ip = 0
      npi = 0
      do i = npstar + 1 , npfin
        ip = ip + 1
        if (skipp(i)) then
          if(opg_root()) write(iwr,*)
     +       'Skip CPHF for perturbation: ',ip
        else
          npi = npi + 1
          if(opg_root()) write(iwr,*)
     +       'CPHF for perturbation: ',ip 
          iatms(npi) = (i-1)/3+1
        end if
      end do
c
c  allocate distributed and replicated workspaces 
c
      if (pg_create_inf(0,nov,maxc*npi,'u',nov,1,g_u,fnm,snm,
     &                  IGMEM_DEBUG)) then
         call pg_distribution(g_u, mynode, ilo, ihi, jlo, jhi)
         write(iwr,918)  g_u, mynode, ilo, ihi, jlo, jhi
        il_u = ilo
        ih_u = ihi
        jl_u = jlo
        jh_u = jhi
        call pg_zero(g_u)
      else
        print *,'**GA error** failed to create GA ',g_u 
      end if

c gdf:   GA for solution vectors
c      if (pg_create(0,nov,maxc*npi,'solnv',nov,1,g_s)) then
c         call pg_distribution(g_s, mynode, ilo, ihi, jlo, jhi)
c         write(iwr,918)  g_s, mynode, ilo, ihi, jlo, jhi
c        il_s = ilo
c        ih_s = ihi
c        jl_s = jlo
c        jh_s = jhi
c        call pg_zero(g_s)
c      else
c        print *,'**GA error** failed to create GA ',g_s 
c      end if

      i00 = igmem_alloc_inf(ncoorb,fnm,snm,"e",IGMEM_DEBUG)
      i10 = igmem_alloc_inf(nov*npi,fnm,snm,"u",IGMEM_NORMAL)
      i20 = igmem_alloc_inf(nov*npi,fnm,snm,"unxt",IGMEM_NORMAL)
      i30 = igmem_alloc_inf(nov*npi,fnm,snm,"prhs",IGMEM_NORMAL)
      i40 = igmem_alloc_inf(maxc*npi,fnm,snm,"b",IGMEM_DEBUG)
      i50 = igmem_alloc_inf(maxc*npi,fnm,snm,"cc",IGMEM_DEBUG)
      i60 = igmem_alloc_inf(maxc*npi,fnm,snm,"uu",IGMEM_DEBUG)
      i70 = igmem_alloc_inf(maxc*maxc*npi,fnm,snm,"uau",IGMEM_DEBUG)
      lb = max( nvv , max( 4*maxc*npi , nov*npi ))
      i80 = igmem_alloc_inf(lb,fnm,snm,"buf",IGMEM_NORMAL)
      i90 = igmem_alloc_inf( max( nvv , nov ),fnm,snm,"buf2",
     &                       IGMEM_DEBUG)
      iz = igmem_alloc_inf(nov*npi,fnm,snm,"rhs",IGMEM_NORMAL)

c  fetch perturbations
      ip = 0
      npi = 0
      do i = npstar + 1 , npfin
         ip = ip + 1
         if (skipp(i)) then
c  skip perturbation
         else
            npi = npi + 1
            call fetch(Q(iz+(npi-1)*nov),nov,'rhs',ip)
         end if
      end do

c  in-core solver
            call secdchf_ga2(Q(i00),Q(iz),
     &           Q(i10),Q(i20),Q(i30),
     &           Q(i40),Q(i50),
     &           Q(i60),Q(i70),Q(i80),
     &           Q(i90),maxc,nov,eps,npi,iatms,Q(1))

c  stash  perturbations
      ip = 0
      npi = 0
      do i = npstar + 1 , npfin
         ip = ip + 1
         if (skipp(i)) then
            call vclr(Q(i80),1,nov)
            call stash(Q(i80),nov,'soln',ip)
         else
            npi = npi + 1
            call stash(Q(iz+(npi-1)*nov),nov,'soln',ip)
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

      call pg_destroy_inf(g_u,'u',fnm,snm)

      return
 6010 format (//1x,'insufficient store for chf equations'//1x,
     +        'store available ',i8/1x,'required - at least ',i8,
     +        ' and preferably ',i8)
918   format('secdhcf: distribution of GA [',i6,'] to node ',1i4,':'
     &    ,4x,'ilo =',1i6,3x,'ihi =',1i6,3x,'jlo =',1i6,3x,'jhi =',1i6)
      end
c
c gdf:  GA version with loop over perturbations inside
c  nb.  this version currently has nov*np memory requirements
c    --   aim to reduce this!
c
      subroutine secdchf_ga2(
     &     e,       ! ncoorb           i00
     &     rhs,     ! nov*npi          iz
     &     u,       ! nov*npi          i10
     &     unxt,    ! nov*npi          i20
     &     prhs,    ! nov*npi          i30
     &     b,       ! maxc*npi         i40
     &     cc,      ! maxc*npi         i50
     &     uu,      ! maxc*npi         i60
     &     uau,     ! maxc*maxc*npi    i70  
     &     buf,     ! nov*npi          i80
     &     buf2,    ! nvv              i90
     &     maxc,    !                  #iters
     &     nov,
     &     eps,            
     &     npi,
     &     iatms,   ! npi              atom number for each perturbation
     &     q)

      implicit REAL  (a-h,o-z)
INCLUDE(common/timez)
INCLUDE(common/common)
INCLUDE(common/cndx41)
INCLUDE(common/prnprn)
INCLUDE(common/iofile)
INCLUDE(common/timeperiods)
      integer a, npi
      common/small/alpha(50,50),aa(50,50),wk1(50),wk2(50)
      dimension e(*), rhs(nov,npi)
      dimension u(nov,npi),unxt(nov,npi),prhs(nov,npi)
      dimension b(maxc,npi),cc(maxc,npi),uu(maxc,npi),uau(maxc,maxc,npi)
      dimension buf(*), buf2(*)
      dimension eps(nov)
      dimension iatms(npi)
      dimension q(*)
      character*10 charwall
INCLUDE(common/global)
      logical opg_root 
      data smal/1.0d-13/,tich/1.0d-24/
      data four/4.0d0/,zero/0.0d0/

      call start_time_period(TP_MP2CHF)

      no = nocca
      nv = nvirta
      n = ncoorb
      no1 = no + 1
      nvv = nv*nv
      novnpi = nov*npi
c
c  iterative solution of simultaneous equations to give
c  coupled hartree-fock first order wavefunction
c
      nnodes = ipg_nnodes()
      mynode = ipg_nodeid()
CMR   write(iwr,*) 'CMR: secdchf_ga2 nnodes=',nnodes,' mynode=',mynode
      if (maxc.gt.50) maxc = 50
      uconv = 10.0d0**(-iconvv)
      iaa = 50
      ifail = 0
      if (oprn(12)) write (iwr,6010)
      if (oprn(13)) write (iwr,6020)

      call vclr(uau,1,maxc*maxc*npi)
      call vclr(b,1,maxc*npi)
      call vclr(uu,1,maxc*npi)

      write(iwr,6000) cpulft(1) ,charwall()

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  get zeroth order estimate
CMR   write(iwr,*) 'CMR: secdchf_ga2 npi=',npi,' nov=',nov
      do j = 1 , npi
        do i = 1 , nov
          v = -rhs(i,j)*eps(i)
          if (dabs(v).le.tich) v = zero
          b(1,j) = b(1,j) + v*v
          u(i,j) = v
        end do
      end do

c  copy U into initial (local) locns of global array
      do j = 1 , npi
        i = (j-1)*maxc + 1 
        if (i.ge.jl_u.and.i.le.jh_u)
     &     call pg_put(g_u, 1, nov, i, i, u(1,j), nov)
      end do

      if (oprn(6)) write(iwr,6100)
c  start of iterative solution of chf equations
c  50 iterations are allowed ---  usually less than 10 are necessary

      do 1 itr = 1 , maxc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        call vclr(unxt,1,novnpi)

c  parallel hessian*trial-vector product
        call mxmau_ga_secd2(u,unxt,buf,buf2,no,nv,nov,npi)

c  gsum product vector
c gdf:  loop over perturbations outside ?
c gdf:   (i.e. recompute orbital-hessian blocks but save on replicated w/s)
c       do i = 1 , npi
        call pg_dgop(101+itr,unxt,novnpi,'+')
c       end do
_IF(ccpdft)
c
c       Add the DFT contributions
c
        call au_dft(q,u,unxt,npi,iatms)
_ENDIF

c        do j = 1 , npi
c          call dbg_wrt_buf(unxt(1,j),nov,'ux pg ',j)
c        end do

c  scale unxt by difference of eigenvalues
        do j = 1 , npi 
          do i = 1 , nov
            unxt(i,j) = unxt(i,j)*eps(i)
          end do
        end do

c  norm of initial trial vectors
        do j = 1 , npi
          uu(itr,j) = ddot(nov,u(1,j),1,u(1,j),1)
c       write(6,*)'me ',ipg_nodeid(),' pertbn ',j,'uu: ',uu(itr,j)
        end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  norm of previous (local) trial vectors with current
        do j = 1 , npi
          do i = 1 , itr - 1
            ij = (j-1)*maxc + i
            if (ij.ge.jl_u.and.ij.le.jh_u) then
              call pg_get(g_u, 1, nov, ij, ij, buf, nov)
              if (i.lt.itr-1)
c    &            uau(itr,i,j) = ddot(nov,u(1,j),1,buf2,1)
c               uau(i,itr,j) = ddot(nov,buf2,1,unxt(1,j),1)
     &            uau(itr,i,j) = ddot(nov,u(1,j),1,buf,1)
                uau(i,itr,j) = ddot(nov,buf,1,unxt(1,j),1)
            end if
          end do
        end do

        if (opg_root()) then
c  one node does diagonal elements of uau
          do j = 1 , npi
            uau(itr,itr,j) = ddot(nov,u(1,j),1,unxt(1,j),1)
          end do
c  current u hasnt been saved to g_u yet 
          if (itr.gt.1) then 
            do j = 1 , npi
             uau(itr,itr - 1,j) = uu(itr,j)
            end do
          end if
        end if

c  pack message buffer  
        ij = 0 
        do j = 1 , npi
          do i = 1 , itr-1
            ij = ij + 1
            buf(ij) = uau(itr,i,j)
            ij = ij + 1
            buf(ij) = uau(i,itr,j)
          end do
          ij = ij + 1
          buf(ij) = uau(itr,itr,j)
        end do

        l999 = npi*(itr+itr-1)
c        call dbg_wrt_buf(buf,l999,'uau bg',itr)

c  global sum
        call pg_dgop(202+itr,buf,l999,'+')

c        call dbg_wrt_buf(buf,l999,'uau ag',itr)

c  unpack message buffer
       ij = 0 
       do j = 1 , npi
         do i = 1 , itr-1
           ij = ij + 1
           uau(itr,i,j) = buf(ij)
           ij = ij + 1
           uau(i,itr,j) = buf(ij)
         end do
         ij = ij + 1
         uau(itr,itr,j) = buf(ij)
       end do

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        do k = 1 , npi
          do i = 1 , itr
            do j = 1 , itr
              alpha(j,i) = uau(j,i,k)
            end do
          end do
          do i = 1 , itr
            alpha(i,i) = alpha(i,i) + uu(i,k)
          end do
c  nag routine to solve a small set of simultaneous equations
          call f04atf(alpha,iaa,b(1,k),itr,cc(1,k),aa,iaa,wk1,wk2,ifail)
c
c          call dbg_wrt_buf(cc(1,k),itr,'cc re ',k)
        end do

c  form new solution vector
        call dcopy(novnpi,rhs,1,prhs,1)
        call vclr(rhs,1,novnpi)

        do k = 1 , npi
          do i = 1 , itr
            ccjn = cc(i,k)
            if (dabs(ccjn).gt.tich) then
              ij = (k-1)*maxc + i
              if (ij.ge.jl_u.and.ij.le.jh_u) then
                call pg_get(g_u, 1, nov, ij, ij, buf, nov)
                call daxpy(nov,ccjn,buf,1,rhs(1,k),1)
              end if
            end if
          end do
        end do


        call pg_dgop(303+itr,rhs,novnpi,'+')

c        do j = 1 , npi
c          call dbg_wrt_buf(rhs(1,j),nov,'solnv ',j)
c        end do

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c  check for convergence of solution vectors
c  by rms difference of current and previous
        if (itr.ne.1) then
          gmax = zero
          do j = 1 , npi
            call vsub(rhs(1,j),1,prhs(1,j),1,prhs(1,j),1,nov)
            gnorm = ddot(nov,prhs(1,j),1,prhs(1,j),1)/dfloat(nov)
            gnorm = dsqrt(gnorm)
            gmax = dmax1(gmax,gnorm)
          end do
          if (oprn(6)) write (iwr,6050) itr , gmax
          if (gmax.le.uconv) then
            write (iwr,6090)
            its = itr
            go to 3
          end if
        end if

c  update next expansion vector
        call vclr(prhs,1,novnpi)
        do j = 1 , npi
          do i = 1 , itr
            ij = (j-1)*maxc + i
            if (ij.ge.jl_u.and.ij.le.jh_u) then
              fac = -uau(i,itr,j)/uu(i,j)
c             call pg_get(g_u, 1, nov, ij, ij, buf2, nov)
              call pg_get(g_u, 1, nov, ij, ij, buf, nov)
              call daxpy(nov,fac,buf,1,prhs(1,j),1)
            end if
          end do
        end do

        call pg_dgop(404+itr,prhs,novnpi,'+')

c        do j = 1 , npi
c          call dbg_wrt_buf(prhs(1,j),nov,'new v ',j)
c        end do

        call vadd(prhs,1,unxt,1,unxt,1,novnpi)

        do j = 1 , npi
          do i = 1 , nov
            if (dabs(unxt(i,j)).le.tich) unxt(i,j) = zero
          end do
        end do

c  check for convergence of next expansion vector
c  - rms of elements
        gmax = zero
        do j = 1 , npi
          gnorm = ddot(nov,unxt(1,j),1,unxt(1,j),1)/dfloat(nov)
          gnorm = dsqrt(gnorm)
          gmax = dmax1(gmax,gnorm)
        end do
        if (oprn(6)) write (iwr,6060) gmax
        if (gmax.le.smal) then
          write (iwr,*)'converged - new expansion vector negligible '
          its = itr
          go to 3
        end if

        nxtr = itr + 1

c  save current expansion vector
        do j = 1 , npi
          i = (j-1)*maxc + nxtr 
          if (i.ge.jl_u.and.i.le.jh_u) 
     &      call pg_put(g_u, 1, nov, i, i, unxt(1,j), nov)
        end do

        call dcopy(novnpi,unxt,1,u,1)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

1     continue
c  end of loop

      write (iwr,*) ' no full convergence after ',maxc,' iterations '
      write (iwr,*) ' this will require changes to the program ! '
3     call timit(3)
      write (iwr,6070) itr , cpulft(1) ,charwall()

      if (oprn(13)) then 
        do j = 1 , npi
          write (iwr,6080) j,(rhs(i,j),i=1,nov)
        end do
      end if

      call end_time_period(TP_MP2CHF)

      return
 99   format(1x,8f10.5)
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
      subroutine mxmau_ga_secd2(a1,a2,buf1,buf2,no,nv,nov,npi)
      implicit REAL  (a-h,o-z)
      dimension a1(nov,npi), a2(nov,npi), buf1(*), buf2(nv,*)
INCLUDE(common/global)
      integer i,j,k,l,kl,no,nv,nvv,ij,o1,o2,v1,v2,npi,nvt

      nvt = nv*(nv+1)/2
      nvv = nv*nv

_IF(ccpdft)
      hf_wght = CD_HF_exchange_weight()
_ENDIF

      ij = 0
      do o1=1,no
         do o2=1,o1
            ij = ij + 1
            if (ij.ge.jl_vvoo.and.ij.le.jh_vvoo) then

c  orbital hessian generated on-the-fly
              call pg_get(g_vovo,1,nvv,ij,ij,buf1,nvv)
              do k=1,nv
                do l=1,nv
_IF(ccpdft)
                  buf2(l,k)=4.0d0*buf1((l-1)*nv+k)
     &                     -hf_wght*buf1((k-1)*nv+l)
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
       
c  loop over perturbations 
            do ip = 1 , npi
              do v1 = 1,nv
                do v2 = 1,nv
                  ix1 = (v1-1)*no + o1
                  ix2 = (v2-1)*no + o2
                  a2(ix2,ip) = a2(ix2,ip) + buf2(v2,v1) * a1(ix1,ip)
                  if( o1 .ne. o2 )then
                   a2(ix1,ip) = a2(ix1,ip) + buf2(v2,v1) * a1(ix2,ip)
                  endif
                enddo
               enddo
            end do

          endif	
        enddo
      enddo
      return
      end

c
      subroutine dbg_wrt_buf(buf,len,string,tag)
c
c  parallel print (debug) routine - all lines of buf 
c  labelled with processor-id on output
c
      implicit none
      integer i,ist,ind,len,mynode,ipg_nodeid,tag
      REAL buf(*)
      character*6 string
      mynode = ipg_nodeid()
      ist = 1
      ind = 10
      if(ind.gt.len) ind = len
1     write(6,2) mynode, string, tag, (buf(i),i=ist,ind)
      ist = ist + 10
      ind = ind + 10
      if(ind.gt.len) ind = len
      if(ist.gt.len) return
      goto 1
2     format(1i4,2x,1a6,2x,1i3,2x,10f10.6)
      end

      subroutine ver_secdchf(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/secdchf.m,v $
     +     "/
      data revision /"$Revision: 6176 $"/
      data date /"$Date: 2010-08-10 16:49:47 +0200 (Tue, 10 Aug 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
