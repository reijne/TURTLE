c
c  Adaptation of Intel EISCUBE diagonaliser for GAMESS-UK
c
c
      subroutine gms_pdiag(n, a, v, d, ierr)
      implicit none
INCLUDE(common/timeperiods)
      REAL a(*), v(*), d(*) 
      integer myid,p,n,m,job,ec_numcols,iprt,idiag,ierr
      integer sqmat,subdg,sigma,work1,work2
      integer cindx(1010) ,i
      integer ipg_nodeid, ipg_nnodes, igmem_alloc

INCLUDE(common/vcore)
      call start_time_period(TP_PDIAG)
      if (n .le. 1) return
      
      p = ipg_nnodes()
      myid = ipg_nodeid()
      m = ec_numcols(myid, n, p) 

c allocate workspace 

      sqmat = igmem_alloc(n*m)
      call vclr(Q(sqmat),1,n*m)

      subdg = igmem_alloc(n)
      call vclr(Q(subdg),1,n)

      sigma = igmem_alloc(n)
      call vclr(Q(sigma),1,n)

      if(p.gt.1) then
        i = max(2*(n+1),(m+1)*(m+1))
        work1 = igmem_alloc(i)
        work2 = igmem_alloc((m+1)*(m+1))
        call vclr(Q(work1),1,i)
        call vclr(Q(work2),1,(m+1)*(m+1))
      else if(p.eq.1) then 
        work1 = igmem_alloc(2*(n+1))
        call vclr(Q(work1),1,2*(n+1))
c dummy alloc        
        work2 = igmem_alloc(1)
      else 
        call caserr('internal error in pdiag')
      endif

c flags for analysis and output 
      job = 2
      idiag = 0 
      iprt = 0
c
      if (iprt.ne.0 .and. myid.eq.0) then
        write(*,'(1x,''intel ssd parallel eiscube demo,  05/06/92'')')
        write(*,'(1x,''symmetric matrix eigenvalues.'')')
        write(*,*)
        write(*,'(1x,''input:'')')
        write(*,'(1x,''n = rows of matrix '')') 
        write(*,*)
        write(*,'(1x,''output:'')')
        write(*,'(1x,''secs = node execution time in seconds'')')
        write(*,'(1x,''accuracy on trace, residual, orthogonality'')')
        write(*,*)
        write(*,9)
    9 format(3x,'task',8x,'p',4x,'n',10x,'secs',
     >       6x,'error   ',4x,'ok?')
      endif

c gdf:  construct column strip for this node  
      call ec_strip(n, a, Q(sqmat))

c      write(*,*)'   *** ec_driver ***'

      call ec_driver(p,n,m, d,v,v, 
     &     Q(sqmat),Q(subdg),Q(sigma),
     &     Q(work1),Q(work2),
     &     job,idiag,iprt,cindx) 

      call gmem_free(work2)
      call gmem_free(work1)
      call gmem_free(sigma)
      call gmem_free(subdg)
      call gmem_free(sqmat)



c gdf:  global sum for vectors on each node 
c gdf:      call ec_sumit(n, v, Q(sqmat))

c gdf:  fill matrix by strip-swapping    
      call ec_fill(n, v, cindx) 

      call end_time_period(TP_PDIAG)
      ierr = 0

      return
      end

      subroutine ec_driver(p,n,m, d,z,v,
     &     a,e,sigma,
     &     w1,w2, 
     &     job,idiag,iprt,ind)
      implicit none 
      integer p,n,m,job,idiag,myid,itcnt,ind(n),k,iprt
      integer i,j,l,root,info,mm,ierr,zero,idim
      REAL a(n,m),d(n),e(n),sigma(n),z(m,n),v(n,m),w1(*),w2(*)
      REAL en,ec_epslon,eps1,lb,ub,mesg(3),secs,err,err2,lwork(3)
      REAL tol1,tol2,times(4),dclock,click,vtime
      REAL trace,norma,t,resid,ortho
      equivalence(mesg(1),secs)
      equivalence(mesg(2),err)
      equivalence(mesg(3),err2)
      character*10 ok
      integer ipg_nodeid
      data zero/0/,root/0/

      myid = ipg_nodeid()


      idim=0
      en = n
      tol1 = ec_epslon(en)
      tol2 = 1000.0d0*tol1
      if(iprt.ne.0) then
	call pg_synch(100)
        click = dclock()
      endif
      trace = 0.0d0
      j = myid + 1
      do 15 l = 1, m
        trace = trace + a(j,l)
        j = j + p
   15 continue
      call vclr(d,1,n)
      call vclr(e,1,n)

      if (iprt.ne.0) then
         call pg_dgop(1001,trace,1,'+')
         secs = dclock() - click
         if (myid .eq. root) then
            write(*,720) p,n,secs
 720        format(/,1x,'generate: ',2i5,3x,f11.5)
         endif
         call pg_synch(100)
         click = dclock()
      endif
c------------------------------------------------------------------
c
c     tridiagonal factorization of a, accumulate transformations in z.
c
      call ec_psytre (a,n,n,m,p,myid,d,e,z,w1)

c------------------------------------------------------------------

      if (iprt.ne.0) then 
         secs = dclock() - click

         call pg_dgop(1002,secs,1,'+')
         secs = secs/p
c  check trace and compute f-norm of tridiagonal.
         norma = d(1)**2
         t = d(1)
         do 20 i = 2, n
            t = t + d(i)
            norma = norma + d(i)**2 + 2*e(i)**2
 20      continue
         norma = dsqrt(norma)
         err = dabs(trace-t)/norma
         if (myid .eq. root) then
            if (err .lt. tol1) then
               ok = 'ok        '
            elseif (err .lt. tol2) then
               ok = 'suspicious'
            else
               ok = 'trouble       !  '
            endif
            write(*,730) p,n,secs,err,ok
 730        format(1x,'tridiag:  ',2i5,3x,f11.5,1pd16.4,3x,a10)
            times(1) = secs
         endif
         call pg_synch(1003)
         click = dclock()
      endif 

      eps1 = 0.0d0
      do 30 i = 1, n
         sigma(i) = 0.0d0
         w1(i) = e(i)**2
   30 continue


c------------------------------------------------------------------
c
c     find shifts by bisection
c     mm = index of first eigenvalue for this node
c
      mm = myid*(n/p) + min(myid,mod(n,p)) + 1
      if (m.gt.0) then 
         call ec_tridib(n,eps1,d,e,w1,lb,ub,mm,m,sigma(mm),
     >        ind,ierr,w1,w1(n+2))
      endif 
      call pg_dgop(1004,sigma,n,'+')

c------------------------------------------------------------------


      if (iprt.ne.0) then
         if (ierr .ne. 0) write(*,*) ' nonzero ierr from tridib'
         secs = dclock() - click
         call pg_dgop(1005,secs,1,'+')
         secs = secs/p
         if (myid .eq. root) then
            t = 0.0d0
            do 35 i = 1, n
               t = t + sigma(i)
 35         continue
            err = dabs(trace-t)/norma
            if (err .lt. tol1) then
               ok = 'ok        '
            elseif (err .lt. tol2) then
               ok = 'suspicious'
            else
               ok = 'trouble    !  '
            endif
            write(*,40) p,n,secs,err,ok
 40         format(1x,'values:   ',2i5,3x,f11.5,1pd16.4,3x,a10)
            times(2) = secs
         endif
         call pg_synch(1006)
         click = dclock()
      endif

c------------------------------------------------------------------
c
c     eigenvectors by perfect shift qr with accumulation.
c
      call ec_psytqr (d,e,sigma,n,z,m,w1,job,info)


c------------------------------------------------------------------
c  d has the eigenvalues - z has m eigenvectors

      if (iprt.ne.0) then
         if (info .ne. 0) write(*,*) ' nonzero ierr from psytqr'
         secs = dclock() - click
         call pg_dgop(1006,secs,1,'+')
         secs = secs/p
         if (myid .eq. root) then
            write(*,50) p,n,secs
 50         format(1x,'vectors:  ',2i5,3x,f11.5)
            times(3) = secs
            itcnt = 0 
            do 57 k = 1, n
               itcnt = itcnt + k*e(k)
 57         continue
         endif
         call pg_synch(1007)
         click = dclock()
      endif

c------------------------------------------------------------------
c  convert row distribution z to column distribution v.
c  z and v are actually same storage, just different leading 
c  dimensions.
c  copy z->a

      call dcopy(m*n,z,1,a,1)

      if (p.gt.1) call ec_flip (myid,n,m,p,a,v,w1,w2)

c------------------------------------------------------------------

      if(iprt.ne.0) then
         secs = dclock() - click
         call pg_dgop(1007,secs,1,'+')
         secs = secs/p
         if (myid .eq. root) then
            write(*,60) p,n,secs
 60         format(1x,'flip:     ',2i5,3x,f11.5,/)
            times(4) = secs
         endif
         if (idiag.ne.0) then 
            click = dclock()
c  check residual and orthogonality
            call ec_psyxxx (a,n,n,m,p,myid,d,v,resid,ortho,e,w1)
            secs = dclock() - click
            err = resid/norma
            err2 = ortho
            call pg_dgop(1008,mesg,3,'+')
            secs = secs/p
            err = err/p
            err2 = err2/p
            if (myid .eq. root) then
c  total eigenvalue/eigenvector time
               vtime = times(1) + times(2) + times(3) + times(4)
               write(*,65) p,n,vtime
 65            format(1x,'eigentime:',2i5,3x,f11.5,/)
c  accuracy
               if (err .lt. tol1) then
                  ok = 'ok        '
               elseif (err .lt. tol2) then
                  ok = 'fair      '
               else
                  ok = 'poor      '
               endif
               write(*,70) p,n,secs,err,ok
 70            format(1x,'residual: ',2i5,3x,f11.5,1pd16.4,3x,a10)
               if (err2 .lt. tol1) then
                  ok = 'ok        '
               elseif (err2 .lt. tol2) then
                  ok = 'fair      '
               else
                  ok = 'poor      '
               endif
               write(*,71) p,n,err2,ok
 71            format(1x,'ortho:    ',2i5,3x,11x,1pd16.4,3x,a10)
               write(*,72) p,n,job,times,itcnt,err,err2
 72            format(i4,i6,i4,4f10.3,i10,1p2d16.4)
            endif
         endif
      endif 

      return
      end
c
c   ec_strip
c
c   this routine converts from complete triangle to column-wrapped strip 
c
      subroutine ec_strip(n, a, w)
      implicit none
      integer myid,n,m,p,i,j,ec_numcols,ic,k,jic
      integer ipg_nnodes, ipg_nodeid
      REAL a(*), w(*) 

c     parameters:
c     n - order of the matrix a
c     a - input triangle 
c     w - output strip

      p = ipg_nnodes()
      myid = ipg_nodeid()
      m = ec_numcols(myid, n, p)

      jic = 0 
      i = myid + 1 
      do ic = 1 , m
         do j = 1 , n
            if (i .gt. j) then  
               k = i*(i-1)/2 + j  
            else 
               k = j*(j-1)/2 + i 
            endif 
            jic = jic + 1
            w(jic) = a(k)
         enddo
         i = i + p 
      enddo
      return 
      end 
c     
c     ec_sumit
c     
c  this routine gathers the whole matrix on each node, 
c  using global sum, comm.s n-squared, workspace (b) also n*n. 
c
      subroutine ec_sumit(n, a, b)
      implicit none
      integer i,j,n,m,p,myid,ec_numcols
      integer ipg_nnodes, ipg_nodeid
      REAL a(n,n),b(n,n)

c     parameters:
c     n - order of the matrix a 
c     a - input matrix distributed by columns, in wrap fashion, so node 0
c     contains cols 1, p+1, 2*p+1, etc, as its 1st, 2nd, 3rd, etc, cols. 
c     the a array space is assumed for now to be of order n on input.
c     b - workspace
c
      p = ipg_nnodes()
      myid = ipg_nodeid()
      m = ec_numcols(myid, n, p) 
c
c  zero rest of matrix 
      do i = m+1 , n 
         do j = 1 , n 
            a(j,i) = 0.0d0 
         enddo 
      enddo 
c
c distribute first m cols through matrix  
      j = myid + 1 + (m-1)*p 
      do i = m , 1 , -1 
         call dswap(n,a(1,i),1,a(1,j),1) 
         j = j - p 
      enddo 
c
c  now do global sum
      call pg_dgop(1009,a,n*n,'+')
      return 
      end
c
c  returns the number of cols wrapped onto node
c
      integer function ec_numcols(id, n, p) 
      integer id,n,p,m 
      m = n/p
      if (id .lt. mod(n,p)) m = m + 1
      if (id .ge. p) m = 0
      ec_numcols = m 
      end 
c     
c     ec_fill
c     
      subroutine ec_fill(n, a, cindx)
      implicit none
      integer myid,n,m,p,i,j,k,myswap,icnt,ibf,ec_numcols
      integer mswp,cindx(n)
      REAL a(n,n)
      integer ipg_nnodes, ipg_nodeid
      integer lenmes, nodefrom
     
c     parameters:
c     myid - node number
c     n - order of the matrix a 
c     m - number of columns held by this node 
c     p - number of processors
c     a - input matrix distributed by columns, in wrap fashion, so node 0
c     contains cols 1, p+1, 2*p+1, etc, as its 1st, 2nd, 3rd, etc, cols. 
c     the a array space is assumed for now to be of order n on input.
c
c  this routine (based on ec_flip) gathers the whole matrix on each node, 
c  with the correct column ordering, starting with the matrix column-
c  wrapped in a round-robin fashion over the nodes (see above). the 
c  scheme is more involved than that in ec_sumit (which uses a simple 
c  global sum) but does not require additional workspace for a receive 
c  buffer and the message lengths scale as m*n. a workspace vector of n 
c  integers is required for book-keeping about the column orderings. 
c  the send buffer is therefore always the first m columns of a. 
c  the receive buffer is the next m' columns of a, in the first step, 
c  then the m' columns after that in the second swap-phase, etc, (m' 
c  is the number of columns held on the node you are swapping with). 

      p = ipg_nnodes()
      myid = ipg_nodeid()
      m = ec_numcols(myid, n, p) 

c  (m+1)th column of a is start of the first receive buffer 
      ibf = m + 1  
c  index the columns on this node 
      cindx(1) = myid + 1
      do i = 2 , m 
         cindx(i) = cindx(i-1) + p 
      enddo 
c  generate first power of 2 bigger than p 
c   (this is part of the method for shuffling the nodes into  
c    some sequence suitable for swapping the data in steps)  
      icnt=1
      do while (icnt.lt.p)
         icnt=icnt*2
      enddo
c  loop over range, missing out first, 0 (whose xor with myid gives myid), 
c    and last, (2^n)>p (which always gives an xor bigger than p).  
      do  k = 1, icnt-1
c  find node to swap with this step
         myswap = IXOR32(myid,k)
c  skip nodes bigger than p 
         if (myswap .le. p-1) then
c  mswp is the number of columns on the other node 
	    mswp = ec_numcols(myswap, n, p) 
c  index the columns to be received
            cindx(ibf) = myswap + 1 
	    do i = ibf + 1 , ibf + mswp - 1 
               cindx(i) = cindx(i-1) + p 
            enddo 
c  send/receive columns 
            if (myid .le. myswap) then
               call pg_snd(1010,a(1,1),8*m*n,myswap,1)
               call pg_rcv(1011,a(1,ibf),8*mswp*n,lenmes,-1,
     &              nodefrom,1)
            else
               call pg_rcv(1010,a(1,ibf),8*mswp*n,lenmes,myswap,
     &              nodefrom,1)
               call pg_snd(1011,a(1,1),8*m*n,myswap,1)
            endif
c  point to start of next receive buffer 
            ibf = ibf + mswp 
         endif
      enddo 
c  now re-shuffle columns into correct order 
      do i = 1, n-1  
         if (cindx(i) .ne. i) then  
            do j = i+1, n 
               if (cindx(j) .eq. i) then 
                  call dswap(n,a(1,i),1,a(1,j),1)   
                  cindx(j) = cindx(i)   
                  cindx(i) = i  
               endif 
            enddo 
         endif
      enddo 
      return
      end
c
c  ec_psytre
c
      subroutine ec_psytre (a,lda,n,m,p,id,d,e,z,work)
      implicit none
      integer lda,n,m,id,p
      REAL a(lda,m),d(n),e(n),work(*),z(m,n)
      REAL alpha,beta,gamma,dnrm2,atemp
_IF(hp700)
      REAL `vec_$ddot', `ddot'
_ELSE
      REAL ddot
_ENDIF
      integer i,j,k,l,ld,r,dpsize,type,zero
      integer jjj, ik, idim
      data dpsize/8/, type/1024/, zero /0/
      integer rank, size, ipg_nodeid, ipg_nnodes
      integer lenmes, nodefrom

      rank = ipg_nodeid()
      size = ipg_nnodes()

      ld=id+1
      do 20 i=1,m
      do 10 j=1,n
10    z(i,j)=0.0d0
      z(i,ld)=1.0d0
20    ld=ld+p

      l=1
      ld=id+1
      d(1)=0.0d0
      e(1)=0.0d0
      if (id.eq.0) d(1) = a(1,1)
      do 80 k = 1, n-1
         r = mod(k-1,p)
         if (id .eq. r) then
            atemp = a(ld,l)
            alpha = dnrm2(n-k,a(k+1,l),1)
            if (a(k+1,l) .lt. 0.0d0) alpha = -alpha
            if (alpha .ne. 0.0d0) then
               call dscal(n-k,1.0d0/alpha
     &                   ,a(k+1,l),1)
               a(k+1,l) = a(k+1,l) + 1.0d0
            endif
            do jjj = 0,size-1
               if (jjj .ne. rank) then
                  call pg_snd(1012,a(k+1,l),8*(n-k),jjj,1)
               endif
            enddo
            if (k .ne. 1) then
               ik = k-1
               if (beta .ne. 0.0d0) then
                  do 61 i = 1, m
                     gamma = ddot(n-ik,d(ik+1),1,z(i,ik+1),m)/beta
                     call daxpy(n-ik,gamma,d(ik+1),1,z(i,ik+1),m)
 61               continue
               endif
               if (p .ne. 1) e(ik+1) = 0.0d0
               d(ik+1) = atemp
            endif
            call dcopy(n-k,a(k+1,l),1,d(k+1),1)
            l = l + 1
            ld = ld + p
         else
            call pg_rcv(1012, d(k+1), 8*(n-k), lenmes, r,
     &           nodefrom, 1)
         endif

         beta = -d(k+1)

         if (beta .ne. 0.0d0) then
            call dcopy(n-k,0.0d0,0,e(k+1),1)
            i = ld
            do 30 j = l, m
               e(i) = e(i) + ddot(n-i+1,a(i,j),1,d(i),1)
               call daxpy(n-i,d(i),a(i+1,j),1,e(i+1),1)
               i = i + p
   30       continue
c
c ***********
            call gdcomb(n-k,e(k+1),zero,idim)
c       write(*,'(1i3,12f10.5)')n-k,(e(jjj),jjj=k+1,n)

            call dscal(n-k,1.0d0/beta,e(k+1),1)
            gamma = ddot(n-k,d(k+1),1,e(k+1),1)/(2.0d0*beta)
            call daxpy(n-k,gamma,d(k+1),1,e(k+1),1)

            i = ld
            do 40 j = l, m
               call daxpy(n-i+1,d(i),e(i),1,a(i,j),1)
               call daxpy(n-i+1,e(i),d(i),1,a(i,j),1)
               i = i + p
   40       continue

          if ((id .ne. mod(k,p)) .or. (k .eq. n-1)) then
            do 60 i = 1, m
               gamma = ddot(n-k,d(k+1),1,z(i,k+1),m)/beta
               call daxpy(n-k,gamma,d(k+1),1,z(i,k+1),m)
   60       continue
          endif

         endif

         if ((id .ne. mod(k,p)) .or. (k .eq. n-1)) then
           d(k+1) = 0.0d0
           e(k+1) = 0.0d0
           if (id .eq. mod(k,p)) d(k+1) = a(ld,l)
         endif

         if (id .eq. r) e(k+1) = -alpha
   80 continue

      call gdcomb (n, d, zero,idim)
      call gdcomb (n, e, zero,idim)

      return
      end
c
c tridib.f
c
      subroutine ec_tridib(n,eps1,d,e,e2,lb,ub,m11,m,w,ind,ierr,rv4,rv5) 
c 
      implicit none
      integer i,j,k,l,m,n,p,q,r,s,ii,m1,m2,m11,m22,tag,ierr,isturm 
      REAL d(n),e(n),e2(n),w(m),rv4(n),rv5(n) 
      REAL u,v,lb,t1,t2,ub,xu,x0,x1,eps1,tst1,tst2,ec_epslon 
      integer ind(m) 
c 
c     this subroutine is a translation of the algol procedure bisect, 
c     num. math. 9, 386-393(1967) by barth, martin, and wilkinson. 
c     handbook for auto. comp., vol.ii-linear algebra, 249-256(1971). 
c 
c     this subroutine finds those eigenvalues of a tridiagonal 
c     symmetric matrix between specified boundary indices, 
c     using bisection. 
c 
c     on input 
c 
c        n is the order of the matrix. 
c 
c        eps1 is an absolute error tolerance for the computed 
c          eigenvalues.  if the input eps1 is non-positive, 
c          it is reset for each submatrix to a default value, 
c          namely, minus the product of the relative machine 
c          precision and the 1-norm of the submatrix. 
c 
c        d contains the diagonal elements of the input matrix. 
c 
c        e contains the subdiagonal elements of the input matrix 
c          in its last n-1 positions.  e(1) is arbitrary. 
c 
c        e2 contains the squares of the corresponding elements of e. 
c          e2(1) is arbitrary. 
c 
c        m11 specifies the lower boundary index for the desired 
c          eigenvalues. 
c 
c        m specifies the number of eigenvalues desired.  the upper 
c          boundary index m22 is then obtained as m22=m11+m-1. 
c 
c     on output 
c 
c        eps1 is unaltered unless it has been reset to its 
c          (last) default value. 
c 
c        d and e are unaltered. 
c 
c        elements of e2, corresponding to elements of e regarded 
c          as negligible, have been replaced by zero causing the 
c          matrix to split into a direct sum of submatrices. 
c          e2(1) is also set to zero. 
c 
c        lb and ub define an interval containing exactly the desired 
c          eigenvalues. 
c 
c        w contains, in its first m positions, the eigenvalues 
c          between indices m11 and m22 in ascending order. 
c 
c        ind contains in its first m positions the submatrix indices 
c          associated with the corresponding eigenvalues in w -- 
c          1 for eigenvalues belonging to the first submatrix from 
c          the top, 2 for those belonging to the second submatrix, etc.. 
c 
c        ierr is set to 
c          zero       for normal return, 
c          3*n+1      if multiple eigenvalues at index m11 make 
c                     unique selection impossible, 
c          3*n+2      if multiple eigenvalues at index m22 make 
c                     unique selection impossible. 
c 
c        rv4 and rv5 are temporary storage arrays. 
c 
c     note that subroutine tql1, imtql1, or tqlrat is generally faster 
c     than tridib, if more than n/4 eigenvalues are to be found. 
c 
c     questions and comments should be directed to burton s. garbow, 
c     mathematics and computer science div, argonne national laboratory 
c 
c     this version dated april 1983. 
c 
c     ------------------------------------------------------------------ 
c 
      ierr = 0 
      tag = 0 
      xu = d(1) 
      x0 = d(1) 
      u = 0.0d0 
c     .......... look for small sub-diagonal entries and determine an 
c                interval containing all the eigenvalues .......... 
      do 40 i = 1, n 
         x1 = u 
         u = 0.0d0 
         if (i .ne. n) u = dabs(e(i+1)) 
         xu = dmin1(d(i)-(x1+u),xu) 
         x0 = dmax1(d(i)+(x1+u),x0) 
         if (i .eq. 1) go to 20 
         tst1 = dabs(d(i)) + dabs(d(i-1)) 
         tst2 = tst1 + dabs(e(i)) 
         if (tst2 .gt. tst1) go to 40 
   20    e2(i) = 0.0d0 
   40 continue 
c 
      x1 = n 
      x1 = x1 * ec_epslon(dmax1(dabs(xu),dabs(x0))) 
      xu = xu - x1 
      t1 = xu 
      x0 = x0 + x1 
      t2 = x0 
c     .......... determine an interval containing exactly 
c                the desired eigenvalues .......... 
      p = 1 
      q = n 
      m1 = m11 - 1 
      if (m1 .eq. 0) go to 75 
      isturm = 1 
   50 v = x1 
      x1 = xu + (x0 - xu) * 0.5d0 
      if (x1 .eq. v) go to 980 
      go to 320 
   60 if (s - m1) 65, 73, 70 
   65 xu = x1 
      go to 50 
   70 x0 = x1 
      go to 50 
   73 xu = x1 
      t1 = x1 
   75 m22 = m1 + m 
      if (m22 .eq. n) go to 90 
      x0 = t2 
      isturm = 2 
      go to 50 
   80 if (s - m22) 65, 85, 70 
   85 t2 = x1 
   90 q = 0 
      r = 0 
c     .......... establish and process next submatrix, refining 
c                interval by the gerschgorin bounds .......... 
  100 if (r .eq. m) go to 1001 
      tag = tag + 1 
      p = q + 1 
      xu = d(p) 
      x0 = d(p) 
      u = 0.0d0 
c 
      do 120 q = p, n 
         x1 = u 
         u = 0.0d0 
         v = 0.0d0 
         if (q .eq. n) go to 110 
         u = dabs(e(q+1)) 
         v = e2(q+1) 
  110    xu = dmin1(d(q)-(x1+u),xu) 
         x0 = dmax1(d(q)+(x1+u),x0) 
         if (v .eq. 0.0d0) go to 140 
  120 continue 
c 
  140 x1 = ec_epslon(dmax1(dabs(xu),dabs(x0))) 
      if (eps1 .le. 0.0d0) eps1 = -x1 
      if (p .ne. q) go to 180 
c     .......... check for isolated root within interval .......... 
      if (t1 .gt. d(p) .or. d(p) .ge. t2) go to 940 
      m1 = p 
      m2 = p 
      rv5(p) = d(p) 
      go to 900 
  180 x1 = x1 * (q - p + 1) 
      lb = dmax1(t1,xu-x1) 
      ub = dmin1(t2,x0+x1) 
      x1 = lb 
      isturm = 3 
      go to 320 
  200 m1 = s + 1 
      x1 = ub 
      isturm = 4 
      go to 320 
  220 m2 = s 
      if (m1 .gt. m2) go to 940 
c     .......... find roots by bisection .......... 
      x0 = ub 
      isturm = 5 
c 
      do 240 i = m1, m2 
         rv5(i) = ub 
         rv4(i) = lb 
  240 continue 
c     .......... loop for k-th eigenvalue 
c                for k=m2 step -1 until m1 do -- 
c                (-do- not used to legalize -computed go to-) .......... 
      k = m2 
  250    xu = lb 
c     .......... for i=k step -1 until m1 do -- .......... 
         do 260 ii = m1, k 
            i = m1 + k - ii 
            if (xu .ge. rv4(i)) go to 260 
            xu = rv4(i) 
            go to 280 
  260    continue 
c 
  280    if (x0 .gt. rv5(k)) x0 = rv5(k) 
c     .......... next bisection step .......... 
  300    x1 = (xu + x0) * 0.5d0 
         if ((x0 - xu) .le. dabs(eps1)) go to 420 
         tst1 = 2.0d0 * (dabs(xu) + dabs(x0)) 
         tst2 = tst1 + (x0 - xu) 
         if (tst2 .eq. tst1) go to 420 
c     .......... in-line procedure for sturm sequence .......... 
  320    s = p - 1 
         u = 1.0d0 
c 
         do 340 i = p, q 
            if (u .ne. 0.0d0) go to 325 
            v = dabs(e(i)) / ec_epslon(1.0d0) 
            if (e2(i) .eq. 0.0d0) v = 0.0d0 
            go to 330 
  325       v = e2(i) / u 
  330       u = d(i) - x1 - v 
            if (u .lt. 0.0d0) s = s + 1 
  340    continue 
c 
         go to (60,80,200,220,360), isturm 
c     .......... refine intervals .......... 
  360    if (s .ge. k) go to 400 
         xu = x1 
         if (s .ge. m1) go to 380 
         rv4(m1) = x1 
         go to 300 
  380    rv4(s+1) = x1 
         if (rv5(s) .gt. x1) rv5(s) = x1 
         go to 300 
  400    x0 = x1 
         go to 300 
c     .......... k-th eigenvalue found .......... 
  420    rv5(k) = x1 
      k = k - 1 
      if (k .ge. m1) go to 250 
c     .......... order eigenvalues tagged with their 
c                submatrix associations .......... 
  900 s = r 
      r = r + m2 - m1 + 1 
      j = 1 
      k = m1 
c 
      do 920 l = 1, r 
         if (j .gt. s) go to 910 
         if (k .gt. m2) go to 940 
         if (rv5(k) .ge. w(l)) go to 915 
c 
         do 905 ii = j, s 
            i = l + s - ii 
            w(i+1) = w(i) 
            ind(i+1) = ind(i) 
  905    continue 
c 
  910    w(l) = rv5(k) 
         ind(l) = tag 
         k = k + 1 
         go to 920 
  915    j = j + 1 
  920 continue 
c 
  940 if (q .lt. n) go to 100 
      go to 1001 
c     .......... set error -- interval cannot be found containing 
c                exactly the desired eigenvalues .......... 
  980 ierr = 3 * n + isturm 
 1001 lb = t1 
      ub = t2 
      return 
      end 
c
c psytqr
c
      subroutine ec_psytqr (d,e,sigma,n,x,mm,work,job,info)
      implicit none
      integer mm,n,job,info
      REAL d(n),e(n),sigma(n),x(mm,n),work(*)
c
c     psytqr, parallel symmetric tridiagonal qr.
c     psytqr computes the eigenvalues and, optionally, the eigenvectors
c     of a real symmetric tridiagonal matrix.
c
c     on entry
c
c        d       real*8 (n)
c                the diagonal of a symmetric tridiagonal matrix.
c
c        e       real*8 (n)
c                the offdiagonal of the matrix, e(2) through e(n).
c
c        sigma   real*8 (n)
c                eigenvalue estimates to be used as initial shifts.
c                probably obtained from parallel bisection.
c
c        n       integer
c                the order of the matrix.
c
c        x       real*8 (mm,n)
c                usually, the output of psytre.
c                if the tridiagonal matrix is the primary data,
c                x should be mm rows of the identity matrix.
c
c        mm      integer
c                the number of rows of x .
c
c        work    real*8 (2*n)
c                work is a scratch array, used only if job > 0,
c                in which case its dimension must be at least 2*n.
c
c        job     integer
c                = 0  eigenvalues only.
c                = 1  eigenvalues and eigenvectors, usual case.
c                = 2  eigenvalues and eigenvectors, special case, see below.
c
c     on return
c
c        d       the eigenvalues, in increasing order.
c
c        x       some rows of the corresponding eigenvectors, if requested.
c
c        info    integer
c                = 0, usual case, results are accurate.
c                = nonzero, results may be inaccurate, see below.
c
c     the two special cases, job = 2 and info nonzero, involve tradeoffs
c     between speed and accuracy that apply only to so-called graded
c     matrices whose elements vary over several orders of magnitude in
c     a uniform way.  graded matrices should be permuted so that their
c     large elements are in the upper left corner.  job = 1 initiates a
c     fast version of an algorithm which leads to vectors whose residuals
c     are small relative to the norm of the input matrix.  job = 2
c     initiates a slower version of the algorithm which leads to eigenvectors
c     whose individual components may be more accurate.  the execution times
c     of the two versions may differ by as much as 50 percent.
c
c     in rare situations associated with underflow, the output value of info
c     may be a nonzero value, m, indicating that the iteration for the m-th
c     eigenvalue did not converge satisfactorily.  this is a warning that the
c     computed eigenvalues d(1) through d(m) may be inaccurate.  better
c     results might be obtained by multiplying the matrix by a huge scale
c     factor, say 1.0e30, computing the eigenvalues again, and then dividing
c     the computed eigenvalues by the scale factor.
c
c     eispack 8x. this version dated 01/28/87 .
c     cleve moler, intel scientific computers 
c
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     internal variables
c
      REAL f,c,s,r,t,shift,tst1,tst2
      integer i,j,k,l,m,iter,maxit,phase,its(2)
c
c        f          the element being "chased" by the implicit qr iteration.
c        c,s        cosine and sine of givens transformation.
c        r,t        temporary variables in shift calculation and qr step.
c        shift      the qr shift.
c        tst1,tst2  used in convergence test.
c        i          index, most likely candidate for vectorization.
c        j,k        other indices.
c        l          start of submatrix.
c        m          end of submatrix and index of eigenvalue being found.
c        iter       iteration counter.
c        maxit      maximum number of iterations.
c        phase      there are two phases.  in phase 1, the implicit tridiagonal
c                   qr algorithm with a shift from the lower 2 by 2 is used to
c                   compute a single eigenvalue, without accumulation of
c                   transformations.  if job = 0, only phase 1 is used.
c                   in phase 2, the qr algorithm with a shift obtained from
c                   phase 1 (a "perfect" shift) is used with accumulation of
c                   transformations to obtain eigenvectors.  if job = 1,
c                   only one step of phase 2 is used for each eigenvalue.
c                   this is fast, and produces accurate eigenvectors for most
c                   matrices.  if job = 2, the phase 2 steps are repeated until
c                   a stringent local convergence test is satisfied.  this takes
c                   more time, but produces more accurate eigenvectors for
c                   graded matrices.
c
      maxit = 30*n
      iter = 0
      info = 0
      do 80 m = n, 2, -1
         its(1) = 0
         its(2) = 0
         phase = 1
c
c        save a copy of the current tridiagonal matrix for phase 2
c
         if (job .gt. 0) then
            call dcopy(m,d,1,work,1)
            call dcopy(m,e,1,work(n+1),1)
         endif
   10    continue
c
c           find l so that e(l) is negligible.
c           if the iteration count is too large, neglect e(m).
c           see the comments below about underflow.
c
            do 20 l = m, 2, -1
               tst1 = dabs(d(l-1)) + dabs(d(l))
               tst2 = tst1 + dabs(e(l))
               if (iter .gt. maxit) then
                  tst2 = tst1
                  if (info .eq. 0) info = m
               endif
               if (tst2 .eq. tst1) go to 30
   20       continue
            l = 1
c
c           if e(m) is negligible, then accept d(m) as an eigenvalue.
c
   30       if (l .eq. m) then
               if (job .eq. 0 .or. phase .eq. 2) go to 70
               phase = 2
               shift = d(m)
               call dcopy(m,work,1,d,1)
               call dcopy(m,work(n+1),1,e,1)
               go to 10
            endif
            iter = iter + 1
            its(phase) = its(phase) + 1
c
c           phase 1 shift calculation.
c
            if (phase .eq. 1) then
               t = (d(m-1) - d(m))/(2.0d0*e(m))
               r = sqrt(1.0d0+t*t)
               if (t .lt. 0.0d0) r = -r
               shift = d(m) - e(m)/(t + r)
               if (its(1) .eq. 1) then
                  j = m
                  do 40 i = 1, m
                     if (dabs(sigma(i)-shift) .lt.
     >                  dabs(sigma(j)-shift)) j = i
   40             continue
                  shift = sigma(j)
                  sigma(j) = sigma(m)
                  sigma(m) = shift
               endif
            endif
c
c           implicit qr iteration, chase nonzero  f  down matrix.
c
            e(l) = d(l) - shift
            f = e(l+1)
            do 60 j = l, m-1
               r = sqrt(e(j)*e(j)+f*f)
               if (r .ne. 0.0d0) then
                  c = e(j)/r
                  s = f/r
               endif
               t = s*(d(j+1) - d(j)) + 2.0d0*c*e(j+1)
               e(j) = r
               e(j+1) = e(j+1) - c*t
               t = s*t
               d(j) = d(j) + t
               d(j+1) = d(j+1) - t
               if (j .lt. m-1) then
                  f = s*e(j+2)
c                 underflow in the computation of f jeopardizes convergence,
c                 but it is too late to do anything now, except rescale the
c                 input matrix and start all over.
                  e(j+2) = -c*e(j+2)
               endif
c
c              accumulate transformations during phase 2.
c
               if (phase .eq. 2) then
                  call ec_update(x(1,j),mm,s,c)
c                  do 50 i = 1, mm
c                     t = x(i,j)
c                     x(i,j) = c*t + s*x(i,j+1)
c                     x(i,j+1) = s*t - c*x(i,j+1)
c   50             continue
               endif
   60       continue
            e(l) = 0.0d0
         if (phase .eq. 1 .or. job .eq. 2) go to 10
   70    work(m) = its(1)
         work(n+m) = its(2)
   80 continue
      work(1) = 0
      work(n+1) = 0
c gdf:      if (.true.) return
c
c     sort the eigenvalues, and possibly the eigenvectors.
c
      do 130 j = 1, n-1
         k = j
         do 110 i = j+1, n
            if (d(i) .lt. d(k)) k = i
  110    continue
         if (k .ne. j) then
            call dswap(1,d(j),1,d(k),1)
            call dswap(mm,x(1,j),1,x(1,k),1)
         endif
  130 continue
      return
      end
c     
c   ec_flip
c     
      subroutine ec_flip (myid,n,m,p,ar,ac,w1,w2)
      implicit none
      integer myid,n,m,p
      REAL ar(m,n),ac(n,m),w1(*),w2(*)
c     
c     convert matrix distributed by rows to matrix distributed by columns.
c     both row matrix and column matrix are row mapped. this routine
c     does not assume an even distribution of rows/columns among the
c     nodes. this version uses a large working space w to reduce the
c     number of messages and increase efficiency
c     
c     parameters:
c     
c     myid - node number
c     n - length of rows (initially) and columns (finally). this is the
c     long dimension, which must be the same on all the nodes. it is
c     also the first dimension of the matrix ac.
c     m - length of columns (initially) and rows (finally). this is the
c     short dimension and may vary between nodes when n is not divisible
c     by p. it is also the first dimension of ar.
c     p - number of processors
c     ar - input matrix distributed by rows, in wrap fashion, so node 0
c     contains rows 1, p+1, 2*p+1, etc. ar is unchanged on output.
c     ac - output matrix distributed by columns in wrap fashion so node 0
c     contains columns 1, p+1, 2*p+1, etc.
c     w - work array. its length should be at least (m+1)*(m+1)*dpsize
c     
c     richard chamberlain   4 may 92
c     
      integer ik,j,k
      integer myswap,inum,mswlen
      integer icnt, lenmes, nodefrom

c      integer tflip, pid, dpsize
c      data tflip/8005/, pid/0/, dpsize/8/

c     
c     loop on other nodes to exchange data with

      icnt=1
      do while (icnt.lt.p)
         icnt=icnt*2
      enddo
      do 30 k = 1, icnt-1
c      do 30 k = 1, p-1
c find node to swap with this step
c         myswap = xor(myid,k)
c         call xor(myid,k,myswap)
         myswap = IXOR32(myid,k)
         if (myswap.le.p-1)then

c     inum counts the length of the communication to node myswap
            inum = 0
c     copy data from ar to buffer w1
            do 10 j = myswap+1,n,p
               call dcopy(m,ar(1,j),1,w1(inum+1),1)
               inum = inum + m
 10         continue
c     send and receive data
c
c NB: intel version was asynchronous here
c    call csend(tflip+k,w,inum*dpsize,myswap,pid)
c    call crecv(tflip+k,w,(m+1)*m*dpsize)
c
            if(myid.le.myswap)then
               call pg_snd(1013,w1,8*inum,myswap,1)
               call pg_rcv(1014,w2,8*(m+1)*m,lenmes,myswap,
     &              nodefrom, 1)
            else
               call pg_rcv(1013,w2,8*(m+1)*m,lenmes,myswap,
     &              nodefrom, 1)
               call pg_snd(1014,w1,8*inum,myswap,1)
            endif

c     inum is counter of buffer as we copy from buffer to ac
            inum = 0
c     ik is counter of which column to put data
            ik = 1

c     mswlen is length of columns originally in the swapping node -
c     this will be m-1, m or m+1, but we work it out here
            mswlen =  n/p
            if (myswap .lt. (n - mswlen*p)) mswlen = mswlen + 1
c     copy data into correct place in ac
            do 20 j = myid+1,n,p
               call dcopy(mswlen,w2(inum+1),1,ac(myswap+1,ik),p)
               inum = inum + mswlen
               ik = ik + 1
 20         continue
         endif
 30   continue

c     finally move data, which stays on this node, from ar to ac
      do 40 j = 1,m
         call dcopy(m,ar(j,myid+1),m*p,ac(myid+1+(j-1)*p,1),n)
 40   continue
      return
      end
c
c psyxxx
c
      subroutine ec_psyxxx (a,lda,n,m,p,myid,d,v,resid,ortho,e,w)
      implicit none
      integer lda,n,m,p,myid
      integer zero,idim
      data zero /0/
      REAL a(lda,m),d(n),v(n,m),resid,ortho,e(n),w(n)
c
c     resid = max over j of norm(a*v(:,j) - d(j)*v(:,j), 1)
c     ortho = max over j of norm(v'*v(:,j) - i(:,j), 1)
c
c     this version dated 01/24/87.
c     cleve moler, intel scientific computers.
c
      REAL t,dasum
_IF(hp700)
      REAL `vec_$ddot'
_ELSE
      REAL ddot
_ENDIF
      integer i,j,k,l,r

      resid = 0.0d0
      ortho = 0.0d0
c
      l = 1
      do 30 k = 1, n
         r = mod(k-1,p)
c
c        distribute e = v(:,k)
         
         if (myid .eq. r) call dcopy(n,v(1,l),1,e,1)
         call pg_brdcst(1016,e,8*n,r)
c
c        form w = a*e - d(k)*e in node r
c        form w(n+1) = norm(v'*e - i(:,k)) distributed
c
         do 10 j = 1, m
            w(j) = ddot(n,v(1,j),1,e,1)
   10    continue
         if (myid .eq. r) w(l) = w(l) - 1.0d0
         t = dasum(m,w,1)
         call dcopy(n,0.0d0,0,w,1)
         i = myid+1
         do 20 j = 1, m
c
c   sum e(i)*a(all,j) into w 
c
            call daxpy(n,e(i),a(1,j),1,w,1)
            i = i + p
   20    continue
         w(n+1) = t

         call gdcomb (n+1,w,zero,idim)

         if (myid .eq. r) then
            call daxpy(n,-d(k),e,1,w,1)
            t = dasum(n,w,1)
            resid = dmax1(resid,t)
            t = w(n+1)
            ortho = dmax1(ortho,t)
            l = l+1
         endif
   30 continue
      return
      end
c
c ec_epslon.f
c
      REAL function ec_epslon (x) 
      REAL x 
c 
c     estimate unit roundoff in quantities of size x. 
c 
      REAL a,b,c,eps 
c 
c     this program should function properly on all systems 
c     satisfying the following two assumptions, 
c        1.  the base used in representing floating point 
c            numbers is not a power of three. 
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in floating point variables 
c            that are stored in memory. 
c     the statement number 10 and the go to 10 are intended to 
c     force optimizing compilers to generate code satisfying 
c     assumption 2. 
c     under these assumptions, it should be true that, 
c            a  is not exactly equal to four-thirds, 
c            b  has a zero for its last bit or digit, 
c            c  is not exactly equal to one, 
c            eps  measures the separation of 1.0 from 
c                 the next larger floating point number. 
c     the developers of eispack would appreciate being informed 
c     about any systems where these assumptions do not hold. 
c 
c     this version dated 4/6/83. 
c 
      a = 4.0d0/3.0d0 
   10 b = a - 1.0d0 
      c = b + b + b 
      eps = dabs(c-1.0d0) 
      if (eps .eq. 0.0d0) go to 10 
      ec_epslon = eps*dabs(x) 
      return 
      end 
      
c
c  ec_update.f
c
      subroutine ec_update(x,mm,s,c)
c givens rotation on two part-columns of x. this is used by
c psytqr - see commented out 50 loop in psytqr.f. this file
c has been separated out so that vectorisation can be applied
c here but not to the rest of psytqr.
      implicit none
      integer mm, i
      REAL s,c,t,x(mm*2)
      do 10 i = 1,mm
         t = x(i)
         x(i) = c*t + s*x(i+mm)
         x(i+mm) = s*t - c*x(i+mm)
 10   continue
      return
      end
c
c gdcomb implemented as dgop
c
      subroutine gdcomb(count,sendbuf,zero,idim) 
      implicit none 
      integer count,zero,idim 
      REAL sendbuf(*),recvbuf(*) 
      call pg_dgop(2000,sendbuf,count,'+')
      end

      subroutine ver_pdiag(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/pdiag.m,v $
     +     "/
      data revision /"$Revision: 6131 $"/
      data date /"$Date: 2010-05-17 13:17:36 +0200 (Mon, 17 May 2010) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
