c
c  parallel test routines
c
      subroutine pg_test
      implicit none
      integer ibuff, ibuff2, igmem_alloc
      integer itmp1, itmp2, itmp3
INCLUDE(common/vcore)
INCLUDE(common/gadiis)
INCLUDE(common/parcntl)
      integer size, mult, i, j, power, ii, ix, ierr
      logical opg_root

      if(iptest.eq.0)return
	
      ibuff = igmem_alloc(1000000)
      ibuff2 = igmem_alloc(1000000)
c
c  broadcast
c
      do power=3,6
         size = 10**power
         mult = 1000000/size
         call start_test_timers
         do i=1,mult
            call pg_brdcst(101,Q(ibuff),size,0)
         enddo
         call end_test_timers('broadcast',size,mult)
      enddo

c synchronisation

      size=0
      mult = 10
      call start_test_timers
      do i=1,mult
         call pg_synch(101)
      enddo
      call end_test_timers('sync',size,mult)

c global sum

 101  continue


      do power=3,6

         size = 10**power
         mult = 1000000/size
         call start_test_timers
         do i=1,mult
            call pg_dgop(101,Q(ibuff),size,'+')
         enddo
         call end_test_timers('dgop',8*size,mult)

      enddo

c parallel diag

      itmp1 = igmem_alloc(1000)
      itmp2 = igmem_alloc(1000)
      itmp3 = igmem_alloc(1000000)

      do i=1,10
         size = 100*i
         mult=1

         ix = 0
         do ii =1,size
            do j=1, ii-1
               Q(ibuff + ix) = 0.0
               ix = ix + 1
            enddo
            Q(ibuff + ix) = 1.0
            ix = ix + 1
         enddo

_IF(ga)
         call declare_diis_storage(size,.false.)
         call init_diis_storage(size,Q(ibuff),
     &        Q(ibuff))
_ENDIF

_IF(diag_parallel)
         if(size.le.400)then
         call makeob(size, Q(ibuff),Q(ibuff2),
     &        Q(itmp1),Q(itmp2))
         call start_test_timers
         call gms_pdiag(size,Q(ibuff),
     &        Q(ibuff2),Q(itmp1),ierr)
         if(ierr.ne.0)call caserr('diag error')
         call pg_synch(100)
         call end_test_timers('diag',8*size,mult)
         endif
_ENDIF
_IF(ga)
c save vectors for parallel test
         call load_ga_from_square(ih_vec,Q(ibuff),size)

c GA mult2 
         call start_test_timers
         call mult2_ga(ih_ov,ih_vec, ih_scr, size)
         call pg_synch(100)
         call end_test_timers('p mult2',8*size,mult)
c GA mult2 replication 
         call start_test_timers
         call load_triangle_from_ga(Q(ibuff),ih_ov,size)
         call end_test_timers('ga to tr',8*size,mult)

         call start_test_timers
         call load_square_from_ga(Q(ibuff),ih_ov,size)
         call end_test_timers('ga to sq gop',8*size,mult)

         call start_test_timers
         call pg_get(ih_ov,1,size,1,size,
     &        Q(ibuff),size)
         call end_test_timers('ga to sq get',8*size,mult)

c serial mult2
         call start_test_timers
         call mult2(Q(ibuff2),Q(itmp3),
     &        Q(ibuff),size,size,size)
         call end_test_timers('s mult2',8*size,mult)

         call destroy_diis_storage (.false.)

_ENDIF
      enddo

      call gmem_free(itmp3)
      call gmem_free(itmp2)
      call gmem_free(itmp1)


c diis solver

c orfog
      call gmem_free(ibuff2)
      call gmem_free(ibuff)

      if(opg_root())write(6,*)'loads at end of tests'
      call pg_pproc

      return
      end


      subroutine start_test_timers
      REAL w0,c0,s0,u0,buff(3)
      common/test_timers/w0,c0,s0,u0

      call gms_cputime(buff)
      u0=buff(1)
      s0=buff(2)
      c0=buff(3)
      call walltime(w0)

      return
      end


      subroutine end_test_timers(text,size,mult)
      implicit none
INCLUDE(common/iofile)
      REAL w0,c0,s0,u0,buff(3)
      common/test_timers/w0,c0,s0,u0
      REAL w,c,s,u, fac, fac2
      integer size, mult
      character text*(*)
      logical opg_root
      logical test_timer_banner
      save  test_timer_banner
      data test_timer_banner/.false./

      if(.not.opg_root())return


      if(.not. test_timer_banner)then
         test_timer_banner = .true.
         write(iwr,98)
 98      format(34x,
     &        ' ---------total  times-------- ',
     &        ' ---------times  per call----- ',
     &        ' Mbyte/s ')
         write(iwr,99)
 99      format(1x,'      operation  ',' bytes','    mult ',
     &        '    cpu ','   user ','    sys ','   wall ',
     &        '    cpu ','   user ','    sys ','   wall ',
     &        '   wall ')
      endif


      call gms_cputime(buff)
      u=buff(1)-u0
      s=buff(2)-s0
      c=buff(3)-c0
      call walltime(w)
      w=w-w0
      fac = 1.0d0 / dble(mult)
      fac2 = 1.0d-6 * (size * dble(mult))

      if(size.gt.0)then
         write(iwr,100)text,size,mult,c,u,s,w,
     &        c*fac, u*fac, s*fac, w*fac,
     &        fac2/w
      else
         write(iwr,100)text,size,mult,c,u,s,w,
     &        c*fac, u*fac, s*fac, w*fac
      endif

 100  format(1x,a15,2i8,9f8.3)
 101  format(1x,a15,2i8,8f8.3,'  n/a ')

      return
      end
c
c
      subroutine makeob(nbfn, orbs, work, tmp1, tmp2)
      implicit double precision (a-h, o-z)

      dimension orbs(nbfn,nbfn), work(nbfn,nbfn)
      dimension tmp1(nbfn), tmp2(nbfn)
c
c     generate set of orthonormal vectors by orthoging twice
c     a set of random vectors ... don't do this at home!
c     ... should really diagonalize the overlap to get sym adaption
c
c      if (nodeid().eq.0) then
         call srand48(12345)
         do 10 i = 1, nbfn
            do 20 j = 1, nbfn
               work(j,i) = 0.0d0
               if(i.eq.j)work(i,j)=1.0d0
               orbs(j,i) = drand48(0)
 20         continue
 10      continue
         call orthv2(nbfn, orbs, work, tmp1, tmp2)
         call orthv2(nbfn, orbs, work, tmp1, tmp2)
c      endif
c      call brdcst(99+MSGDBL, orbs, mdtob(nbfn*nbfn), 0)
      end

      subroutine orthv2(n, v, s, work1, work2)

      implicit none
      integer n
      REAL v(n,n), s(n,n), work1(n), work2(n)
c
      integer i
      REAL  a
_IF(hp700)
      REAL `vec_$ddot'
_ELSEIF(cray,t3d)
      REAL sdot
_ELSE
      REAL ddot
_ENDIF
      character*1 xn, xt
      data xn,xt/'n','t'/
c     
c     orthonormalize vectors (v) in place over metric (s)
c     
c     v = matrix of column vectors
c     s = full square metric matrix
c
c     note ... is not suitable for high precision or if the
c              input vectors are close to linear dependence
c
c          ... assumes that the metric is symmetric
c     
c     normalize vector one
c
      call dgemv(xn,n,n,1.0d0,s,n
     &     ,v(1,1),1,0.0d0,work1,1)
      a = ddot(n, work1, 1, v(1,1), 1)
      call dscal(n, 1.0d0/dsqrt(a), v(1,1), 1)
c
      do 10 i = 2, n
c
c     orthog vector i to vectors j=1,...,i-1
c
         call dgemv(xn,n,n,1.0d0,s,n
     &        ,v(1,i),1,0.0d0,work1,1)
         call dgemv(xt,n,i-1,1.0d0,v,n
     &        ,work1,1,0.0d0,work2,1)
         call dgemv(xn,n,i-1,-1.0d0,v,n
     &        ,work2,1,1.0d0,v(1,i),1)
c
c     normalize vector i
c
         call dgemv(xn,n,n,1.0d0,s,n
     &        ,v(1,i),1,0.0d0,work1,1)
         a = ddot(n, work1, 1, v(1,i), 1)
         call dscal(n, 1.0d0/dsqrt(a), v(1,i), 1)
 10   continue
c
      end

      subroutine ver_ptest(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/ptest.m,v $
     +     "/
      data revision /"$Revision: 6131 $"/
      data date /"$Date: 2010-05-17 13:17:36 +0200 (Mon, 17 May 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
