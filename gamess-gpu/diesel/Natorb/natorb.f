      program natorb
      implicit real*8 (a-h,o-y)
      implicit character *8 (z)
      parameter (maxorb=2048)
      common/junk/v,v2,v3,p,occ,eigv,isym,imap
      dimension v(maxorb*maxorb)
      dimension v2(maxorb*maxorb)
      dimension v3(maxorb*maxorb)
      dimension p(maxorb*(maxorb+1)/2)
      dimension occ(maxorb),eigv(maxorb)
      dimension isym(maxorb)
      dimension irrep(8)
      dimension imap(8,maxorb)
      dimension si(2*maxorb),is(maxorb+1)
      common/junkc/zcom(19),ztitle(10)
c
      ind3(i,j)=max0(i,j)*(max0(i,j)-1)/2+min0(i,j)
      ind4(i,j,k)=(i-1)*k+j
c -> i = MO, j = AO (j snel lopend)
c
      read(5,*) num
      read(5,*) ncore,ndiscard
      nactive=num-ncore-ndiscard
      do i=1,num
         read(5,*)jj,isym(i),d1,d2,occ(i)
      enddo
      istep=num/6
      if (mod(num,6).ne.0) istep=istep+1
      do i=1,num
         istart=1
         do k=1,istep
            iend=istart+5
            iend=min(iend,num)
            read(5,'(6F11.6)')(v(ind4(i,j,num)),j=istart,iend)
            istart=istart+6
         enddo
      enddo
      print *,'hf orbitals'
      call prev(v,occ,num,num,num)
      do i=1,8
         irrep(i)=0
      enddo
      do i=ncore+1,ncore+nactive
         j=isym(i)
         irrep(isym(i))=irrep(isym(i))+1
         imap(j,irrep(isym(i)))=i
      enddo
      nirrep=0
      do i=1,8
         if (irrep(i).ne.0) nirrep=nirrep+1
      enddo
      print *,'nirrep',nirrep,(irrep(i),i=1,nirrep)
      do i=1,num*(num+1)/2
         p(i)=0.0d0
      enddo
      do i=1,ncore
         p(ind3(i,i))=2.0d0
      enddo
      nlines=0
      do i=1,nirrep
         nlines=nlines+irrep(i)*irrep(i)
      enddo
      print *,'remco nlines',nlines
      do i=1,nlines
          read(5,*)isyma,isymb,ia,ib,dd
          jj=imap(isyma+1,ia)
          kk=imap(isymb+1,ib)
          p(ind3(jj,kk))=dd
      enddo
      sum=0.0d0
      do i=1,num
         sum=sum+p(ind3(i,i))
      enddo
      print *,'remco density matrix',sum
      call prtri(p,num)
      call jacodiag(p,num,v2,eigv,si,is)
      print *,'NOs on MO basis'
      call prev(v2,eigv,num,num,num)
      call matmult(v,v2,v3,num)
      print *,'NOs on AO basis'
      call prev(v3,eigv,num,num,num)
      print *,'-----'
      do i=1,num
         istart=1
         do k=1,istep
            iend=istart+5
            iend=min(iend,num)
            write(6,'(6F11.6)')(v3(ind4(i,j,num)),j=istart,iend)
            istart=istart+6
         enddo
      enddo
      istart=1
      iend=istart+5
      iend=min(iend,num)
      write(6,'(A6,2x,I4,6F11.6)')'occup ',num,(eigv(j),j=istart,iend)
      istart=istart+6
      do k=2,istep
         iend=istart+5
         iend=min(iend,num)
         write(6,'(6F11.6)')(eigv(j),j=istart,iend)
         istart=istart+6
      enddo
c     call putq(zcomm,ztit,eigv,eigv,num,num,num,0,1,
c    &v3,mpos,iblkq)
      end
   
      subroutine matmult(x,c,cnew,ndim)
      implicit real*8 (a-h,o-z)
      dimension x(ndim,ndim),c(ndim,ndim)
      dimension cnew(ndim,ndim)
      do i=1,ndim
         do j=1,ndim
            cnew(i,j)=0.0d0
            do k=1,ndim
               cnew(i,j)=cnew(i,j)+x(i,k)*c(k,j)
            enddo
         enddo
      enddo
      return
      end

      subroutine prev(v,e,m,n,ndim)
c
c     ----- print out e and v-matrices
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(ndim,*),e(*)
      iwr=6 
      max = 10
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. m) imax = m
      write (iwr,9008) 
      write (iwr,8068) (e(i),i = imin,imax)
      write (iwr,9008) 
      write (iwr,8028) (i,i = imin,imax)
      write (iwr,9008)
      do 150 j = 1,n
  150 write (iwr,8048) j,(v(j,i),i = imin,imax)
 140  if (imax .lt. m) go to 100
      return
 9008 format(/)
 8028 format(17x,10(3x,i3,3x))
 8048 format(i5,12x,10f9.4)
 8068 format(17x,10f9.4)
 9028 format(17x,12(6x,i3,6x))
 9048 format(i5,2x,10x,7f15.10)
 9068 format(17x,7f15.10)
      end

      subroutine prtri(d,n)
c     
c     ----- print out a triangular matrix -----
c     
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension d(*),dd(12)
      iwr=6 
      max = 12
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. n) imax = n
      write (iwr,9008) 
      write (iwr,8028) (i,i = imin,imax)
      do 160 j = 1,n
      k = 0
      do 140 i = imin,imax
      k = k+1
      m = max0(i,j)*(max0(i,j)-1)/2 + min0(i,j)
  140 dd(k) = d(m)
      write (iwr,8048) j,(dd(i),i = 1,k)
  160 continue 
      if (imax .lt. n) go to 100
      return
 9008 format(/)
 9028 format(6x,7(6x,i3,6x))
 9048 format(i5,1x,7f15.10)
 8028 format(6x,12(3x,i3,3x))
 8048 format(i5,1x,12f9.4)
      end
      
      subroutine jacodiag(a,n,evecs,evals,si,is)
c
      implicit real*8   (a-h,o-z) , integer   (i-n)
c....
c.... diagonalises triangular matrix with dimension n eigenvectors
c.... are returned in evecs, eigenvalues in evals, 
c.... si is een scratch array 2*n
c
      dimension a(*),evecs(*),evals(*),si(*),is(*)
c
      small = 1.0d-8
      do i=1,n+1
         is(i)=i*(i-1)/2
      enddo
c
      call jacobi(a,n,is,evecs,n,evals,0,3,small,si)
      return
      end
      subroutine jacobi(h,nbasis,iky,v,nrow,e,init,iorder,small,y)
c
      implicit real*8   (a-h,o-z) , integer   (i-n)
c
c.....
c.....standard jacobi diagonaliser. h is diagonalised, v will contain
c.....the eigenvectors, e the eigenvalues. iky should contain i*(i-1)/2
c.....at i. y is a scratch array of length 2*nbasis.
c.....init   : = 1 ? => v is assumed to be initialised elsewhere
c.....  ,,   : other values will cause it to be a unit matrix at first
c.....iorder : = 0 ? => no ordering in the vectors, no e.val. at all
c.....  ,,   : = 1 ? => no ordering in the vectors, e.val's in e.
c.....  ,,   : = 2 ? => increasing eigenvalues/vectors
c.....  ,,   : = 3 ? => decreasing eigenvalues/vectors
c.....  ,,   : = 4 ? => lock mode (i.e. try to change the initial
c.....                             vectors the least possible)
c.....
      dimension h(*),v(nrow,nbasis),e(*),y(*),iky(*)     
      if (init.ne.1) then
         call zero(v,nrow*nbasis)
         do 10 i=1,nbasis
            v(i,i) = 1.0d0
10       continue
      end if
      if (nbasis.eq.1) then
         e(1)   = h(1)
         return
      end if
11    rlarge = 0.0d0
      do 30 i=2,nbasis
         ii = iky(i)
         do 20 j=1,i-1
           if (dabs(h(ii+j)).gt.rlarge) rlarge = dabs(h(ii+j))
20       continue
30    continue
      if (rlarge.gt.small) then
         reason = rlarge * .1d0
         do 60 i=1,nbasis - 1
            call fmove( h(iky(i)+1),y,i)
            call gather(nbasis-i,y(i+1),h(i+1),iky(i+1))
            hii = y(i)
            do 50 j=i+1,nbasis
               if (dabs(y(j)).ge.reason) then
                  vii  = hii * 0.5d0
                  vij  = y(j)
                  vjj  = h(iky(j)+j) * 0.5d0
                  diff = vii - vjj
                  root = dsqrt( diff*diff + vij*vij )
                  if (diff.lt.0.0d0) root = -root
                  cosr = (diff + root)/vij
                  sinr = dsqrt( 1.0 d0/ (1.0d0+cosr*cosr) )
                  cosr = cosr * sinr
                  diff = vii + vjj
                  call drotj(j-1,y,1,h(iky(j)+1),1,cosr,sinr)
                  hii  = diff + root
                  h( iky(j)+j ) = diff - root
                  y(j) = 0.0d0
                  do 40 k=j+1,nbasis
                     vij   = y(k)
                     kj    = iky(k) + j
                     y(k)  = vij   * cosr + h(kj) * sinr
                     h(kj) = h(kj) * cosr - vij   * sinr
40                continue
                  call drotj(nbasis,v(1,i),1,v(1,j),1,cosr,sinr)
               end if
50          continue
            call fmove(y,h(iky(i)+1),i-1)
            call scatter(nbasis-i,h(i+1),iky(i+1),y(i+1))
            h(iky(i)+i) = hii
60       continue
         goto 11
      end if
      if (iorder.eq.0) then
         return
      else if (iorder.eq.1) then
         call gather(nbasis,e,h,iky(2))
         return
      else if (iorder.eq.2) then
         call gather(nbasis,e,h,iky(2))
         call order(e,v,nbasis,nrow,y,y(nbasis+1), 1)
         call scatter(nbasis,h,iky(2),e)
      else if (iorder.eq.3) then
         call gather(nbasis,e,h,iky(2))
         call order(e,v,nbasis,nrow,y,y(nbasis+1),-1)
         call scatter(nbasis,h,iky(2),e)
      end if
      return
      end
      subroutine gather(n,r,a,map)
c
      implicit real*8 (a-h,o-z), integer (i-n)
      dimension r(n),a(*),map(n)
c
      do 10 loop=1,n
   10 r(loop) = a(map(loop))
c
      return
      end

      subroutine igather(n,r,a,map)
      implicit integer (a-z)
      dimension a(*),r(n),map(n)
c
      do 10 loop=1,n
   10 r(loop) = a(map(loop))
c
      return
      end

      subroutine scatter(n,r,index,a)
c
c...  cray scilib imitation/ but not as in my manual (jvl 1986)
c...  arguments as in cyber205 atmol-library
c
      implicit real*8 (a-h,o-z), integer (i-n)
      dimension r(*),index(n),a(n)
c
      do 10 i=1,n
10    r(index(i)) = a(i)
c
      return
      end
      subroutine order(e,v,n,nrow,y,imap,ipar)
c
      implicit real*8   (a-h,o-z) , integer   (i-n)
c
      dimension e(*),v(nrow,n),y(*),imap(*)
      do 10 i=1,n
         imap(i) = i
10    continue
      call bubble(e,n,ipar,imap)
      do 20 iold=1,n-1
         inew = imap(iold)
         imap(iold) = 0
         if (iold.eq.inew) goto 20
         call fmove(v(1,iold),y,n)
         call fmove(v(1,inew),v(1,iold),n)
         call fmove(y,v(1,inew),n)
         imap(locati(iold,imap,n)) = inew
20    continue
      return
      end
      subroutine bubble(a,n,ipar,imap)
c
      implicit real*8   (a-h,o-z) , integer   (i-n)
c
c.....
c.....simple bubble sort. imap keeps track of the permutations
c.....ipar enables one to use this routine for sorting increasingly
c.....(ipar > 0 => first element smallest) or the other way around
c.....
      dimension a(n),imap(n)
      logical ready
10    ready = .true.
      do 20 i=2,n
         if (a(i-1)*ipar.gt.a(i)*ipar) then
            r         = a(i-1)
            a(i-1)    = a(i  )
            a(i  )    = r
            ready     = .false.
            ii        = imap(i-1)
            imap(i-1) = imap(i  )
            imap(i  ) = ii
         end if
20    continue
      if (.not.ready) goto 10
      return
      end
      subroutine zero(v,n)
c
c...  zero vector
c
      implicit double precision (a-h,o-z), integer(i-n)
      dimension v(*)
c
      do 10 i=1,n
         v(i)=0.0d0
10    continue
c
      return
      end
      subroutine fmove(a,b,n)
c
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(n),b(n)
c
c...  b = a
c
      do 10 i=1,n
         b(i)=a(i)
10    continue
c
      return
      end
      function locati(i,ii,n)
      implicit real*8(a-h,o-z)
      dimension ii(*)
      do 10 j=1,n
         if (ii(j).eq.i) then
            locati = j
            return
         end if
10    continue
      locati = 0
      return
      end

      subroutine drotj(n,sx,incx,sy,incy,sc,ss)
      implicit real*8  (a-h,o-z)
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       sx  double precision vector with n elements
c     incx  storage spacing between elements of sx
c       sy  single precision vector with n elements
c     incy  storage spacing between elements of sy
c       sc  element of rotation matrix
c       ss  element of rotation matrix
c
c     --output--
c       sx  rotated vector sx (unchanged if n .le. 0)
c       sy  rotated vector sy (unchanged if n .le. 0)
c
c     multiply the 2 x 2 matrix  ( sc ss) times the 2 x n matrix (sx**t)
c                                (-ss sc)                        (sy**t)
c     where **t indicates transpose.  the elements of sx are in
c     sx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = (-incx)*n, and similarly for sy using ly and incy.
      dimension sx(*),sy(*)
      data zero,one/0.d0,1.d0/
      if(n .le. 0 .or. (ss .eq. zero .and. sc .eq. one)) go to 40
      if(.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
c
           nsteps=incx*n
           do 10 i=1,nsteps,incx
                w=sx(i)
                z=sy(i)
                sx(i)=sc*w+ss*z
                sy(i)=-ss*w+sc*z
   10           continue
           go to 40
c
   20 continue
           kx=1
           ky=1
c
           if(incx .lt. 0) kx=1-(n-1)*incx
           if(incy .lt. 0) ky=1-(n-1)*incy
c
           do 30 i=1,n
                w=sx(kx)
                z=sy(ky)
                sx(kx)=sc*w+ss*z
                sy(ky)=-ss*w+sc*z
                kx=kx+incx
                ky=ky+incy
   30           continue
   40 continue
c
      return
      end

