      subroutine matre2(value,supers,superh,superg,
     &                  wone,wtwo,nelec,nalfa,
     &                  scr,s,g,dtotal)
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
c..   this thing calculates non-orthogonal matrix-elements
c..   dumb but robust (for debugging)
c..   non-zero contributions are returned in
c..      wone :  1-st order cofactors
c..      h    :  corresponding one-electron integrals
c..      wtwo :  2-nd order cofactors
c..      g    :  corresponding integrals
c..
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
      common/detcop/ir2(melec),ic2(melec)
c.....
c.....for debugging : dets are given by symblo
c.....
c.....
c.....in this routine the matrix-element is calculated ***debug***
c.....
      dimension wone(nelec**2),
     &          s(nelec**2),
     &          wtwo (*),scr (*),
     &          g (*),
     &          supers(*),superh(*),superg(*)
c.....
      logical asing1,asing2,bsing1,bsing2
      value = 0.0d0
      dtotal = 0.0d0
      asing1 = .false.
      asing2 = .false.
      bsing1 = .false.
      bsing2 = .false.
      nbeta = nelec - nalfa
      call makemat(s,nalfa,ir2,ic2,supers)
      call bio2(s,nalfa,ir2,ic2,ira,deta,ipa)
      na = nalfa * nalfa
      call makemat(s(na+1),nbeta,ir2(nalfa+1),ic2(nalfa+1),supers)
      call bio2(s(na+1),nbeta,ir2(nalfa+1),ic2(nalfa+1),irb,detb,ipb)
      nsing = (nalfa - ira) + (nbeta - irb)
      if (nalfa.eq.0) deta=1
      if (nbeta.eq.0) detb=1
      if (nsing.gt.2) return
      if (nalfa-ira.eq.1) asing1 = .true.
      if (nalfa-ira.eq.2) asing2 = .true.
      if (nbeta-irb.eq.1) bsing1 = .true.
      if (nbeta-irb.eq.2) bsing2 = .true.
c.....
c.....four cases are distinguished : - no singularities
c.....                               - one singularity
c.....                               - two singularities in one block
c.....                               - two singularities in two blocks
c.....
      call cofac1(s      ,wone      ,nalfa,ira,deta)
      call cofac1(s(na+1),wone(na+1),nbeta,irb,detb)
      if (nsing.lt.2) then
c.....
c.....   here the one-electron contribution per block is calculated
c.....
         if (.not.bsing1) then
            call makemat(g,nalfa,ir2,ic2,superh)
            value = value + ddot(nalfa*nalfa,g,1,wone,1) * detb
         end if
         if (.not.asing1) then
            call makemat(g,nbeta,ir2(nalfa+1),ic2(nalfa+1),superh)
            value = value + ddot(nbeta*nbeta,g,1,wone(na+1),1) * deta
         end if
      end if
c.....
c.....   two-electron contributions involving one block at a time
c.....
      if (.not.bsing1.and..not.bsing2) then
         call cofac2(nalfa,wtwo,scr,wone,s,ira,deta)
         call maktwo(g,nalfa,superg,ir2,ic2)
         n2a = nalfa*(nalfa-1)/2
         n2a = n2a * n2a
         value = value + ddot(n2a,g,1,wtwo,1) * detb
      end if
      if (.not.asing1.and..not.asing2) then
         call cofac2(nbeta,wtwo,scr,wone(na+1),s(na+1),irb,detb)
         call maktwo(g,nbeta,superg,ir2(nalfa+1),ic2(nalfa+1))
         n2b = nbeta*(nbeta-1)/2
         n2b = n2b * n2b
         value = value + ddot(n2b,g,1,wtwo,1) * deta
      end if
      if (.not.asing2.and..not.bsing2) then
c.....
c.....   two electron contributions involving the two blocks at a time
c.....
         ib = nalfa + 1
         value = value + gmixed(wone,ir2,ic2,nalfa,
     &                          wone(na+1),ir2(ib),ic2(ib),nbeta,
     &                          superg)
      end if
      value = value * ipa * ipb
      if (bsing1) detb=0
      if (asing1) deta=0
      if (nsing.eq.0) then
         dtotal = deta * detb  * ipa * ipb
      else
         dtotal = 0.0d0
      end if
      return
      end
      subroutine bio2(s,ndim,ir,ic,irank,det,ipar)
c
      implicit real*8 (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine the square matrix s is decomposed in
c.....a left-under/diagonal/right-upper matrix ("biorthogonalised").
c.....ir and ic keep track of the permutations that are induced
c.....by the pivot search. singularity is allowed for.
c.....irank will contain the rank of s. det is the product
c.....of the non-zero diagonal elements times the parity change.
c.....only in non-singular cases det in fact is the determinant.
c.....
      parameter (small=1.d-10)
      dimension s(ndim,ndim),ir(ndim),ic(ndim)
      irank = 0
      det = 1.0d0
      ipar = 1
      call pivo2(s,ndim,1,piv,ir,ic,ipar)
      if (dabs(piv).lt.small) then
         det = 0.0d0
         return
      end if
      irank = irank + 1
      tol = small * piv
      do 100 i = 1,ndim-1
         det = det * s(i,i)
         do 30 j = i+1,ndim
c.....
c.....form i-th row of r
c.....
            s(i,j) = -s(i,j)/s(i,i)
c.....
c.....adapt s-matrix to r
c.....
            do 10 k = i+1,ndim
               s(k,j) = s(k,j) + s(k,i)*s(i,j)
10          continue
c.....
c.....adapt r-matrix to itsself
c.....
            do 20 l = 1,i-1
               s(l,j) = s(l,j) + s(l,i)*s(i,j)
20          continue
30       continue
         call pivo2(s,ndim,i+1,piv,ir,ic,ipar)
c.....
c.....form (i+1)-th row of l
c.....
         do 40 m=1,i
            s(i+1,m) = -s(i+1,m)/s(m,m)
40       continue
c.....
c.....adapt l to itsself
c.....
         sum = s(i+1,1)
         do 60 n=1,i
            do 50 nn = n+1,i
               sum = sum + s(i+1,nn)*s(nn,n)
50          continue
            s(i+1,n) = sum
            sum = s(i+1,n+1)
60       continue
         if(dabs(piv).lt.dabs(tol)) then
            if(irank.eq.ndim-2) then
c.....
c.....      the last row of l still must be made !!!!
c.....
               do 70 m=1,irank
                  s(ndim,m) = -s(ndim,m)/s(m,m)
70             continue
               do 90 m=1,irank
                  sum = s(ndim,m)
                  do 80 n=m+1,ndim-1
                     sum = sum + s(n,m) * s(ndim,n)
80                continue
                  s(ndim,m) = sum
90             continue
            end if
            return
         end if
         irank = irank + 1
100   continue
      det = det * s(ndim,ndim)
      return
      end
      subroutine pivo2(a,ndim,n,piv,ir,ic,ipar)
c
      implicit real*8 (a-h,o-z) , integer   (i-n)
c
c.....
c.....this searches for the (absolute) largest element in the submatrix
c.....a(n..ndim,n..ndim) and permutes the rows and columns so that this
c.....element occurs in the (n,n)-th position. ir and ic are arrays that
c.....keep track of the changes. the sign of ipar is changed accordingly
c.....
      dimension a(ndim,ndim),ir(ndim),ic(ndim)
      piv = 0.d0
      irt = n
      ict = n
c.....
c.....search pivot
c.....
      do 20 j = n,ndim
         do 10 i = n,ndim
            if (dabs(piv).lt.dabs(a(i,j))) then
               piv = a(i,j)
               irt = i
               ict = j
            end if
10       continue
20    continue
c.....
      if (irt.ne.n) then
c.....
c.....permute rows
c.....
         do 30 k = 1,ndim
            aa = a(irt,k)
            a(irt,k) = a(n,k)
            a(n,k) = aa
30       continue
         iii = ir(irt)
         ir(irt) = ir(n)
         ir(n) = iii
         ipar = -ipar
      end if
c.....
      if (ict.ne.n) then
c.....
c.....permute columns
c.....
         do 40 l = 1,ndim
            aa = a(l,ict)
            a(l,ict) = a(l,n)
            a(l,n) = aa
40       continue
         iii = ic(ict)
         ic(ict) = ic(n)
         ic(n) = iii
         ipar = -ipar
      end if
      return
      end
      subroutine cofac1(s,w,ndim,irank,det)
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
c.....
c.....in this routine the first order cofactors are calculated.
c.....see: f. prosser and s. hagstrom, int.j.quantum chem. 2, 89 (1968).
c.....s is the biorthogonalised overlap matrix in question. irank is its
c.....rank. det is the product of the non-zero diagonal elements of s. w
c.....is the first order adjugate of s transposed (!).
c.....
      dimension s(ndim,ndim),w(ndim,ndim)
c.....
      if (ndim.lt.1) return
      if (ndim.eq.1) then
c.....
c.....   the first order cofactor of a one-dimensional matrix is 1.
c.....
         w(1,1)=1.0d0
         return
      end if
c.....
      if (irank.lt.ndim-1) then
c.....
c.....   nullity.ge.2
c.....   no non-zero first order cofactors, therefore no one-electron
c.....   contribution to the matrix element.
c.....
         return
      end if
c.....
      if (irank.lt.ndim) then
c.....
c.....   nullity.eq.1
c.....   so that in principle all cofactors are non-zero.
c.....
         do 20 i=1,ndim-1
            w(ndim,i)=s(i,ndim)*det
            w(i,ndim)=s(ndim,i)*det
            do 10 j=1,ndim-1
               w(j,i)=w(ndim,i)*s(ndim,j)
10          continue
20       continue
         w(ndim,ndim)=det
         return
      end if
c.....
c.....nullity.eq.0
c.....the rank of the s-matrix is equal to its dimension,
c.....therefore a straight forward calculation.
c.....
      call vclr(w,1,ndim*ndim)
c.....
c.....calculate (adjd * l) transposed.
c.....
      do 30 i=1,ndim
         w(i,i)=det/s(i,i)
30    continue
      do 50 i=2,ndim
         do 40 j=1,i-1
            w(j,i)=w(i,i)*s(i,j)
40       continue
50    continue
c.....
c.....calculate (adjd * l) transposed * r transposed.
c.....
      do 80 j=1,ndim
         do 70 i=1,ndim
            do 60 k=max(i,j+1),ndim
               w(i,j)=w(i,j)+w(i,k)*s(j,k)
60          continue
70       continue
80    continue
c.....
      return
      end
      function gmixed(wik,irik,icik,nik,
     &                wjl,irjl,icjl,njl,
     &                superg)
c.....
c.....
c.....
      implicit real*8 (a-h,o-z)
      dimension irik(nik),icik(nik),irjl(njl),icjl(njl),
     &          superg(*),wik(nik*nik),wjl(njl*njl)
c
      value = 0.0d0
c.....
c.....find the alpha/beta two-electron integrals
c.....
      mik = 0
      do 100 k=1,nik
         do 90 i=1,nik
            mik = mik + 1
            mjl = 0
            do 80 l=1,njl
               do 70 j=1,njl
                  mjl = mjl + 1
                  ikjl = intpos(irik(i),icik(k),irjl(j),icjl(l))
                  value = value + wik(mik)*wjl(mjl) * superg(ikjl)
70             continue
80          continue
90       continue
100   continue
      gmixed = value
      return
      end
      subroutine makemat(a,ndim,ir,ic,super)
c
      implicit real*8 (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine a is constructed from super using the labels from
c.....ir(ow) and ic(olumn). super is assumed to be a triangular matrix.
c.....
      dimension a(ndim,ndim),ir(ndim),ic(ndim),super(*)
      ind(i,j) = max(i,j) * (max(i,j)-1)/2 + min(i,j)
      it = 0
      do 20 j=1,ndim
         do 10 i=1,ndim
            ipos = ind( ir(i),ic(j) )
            a(i,j) = super(ipos)
10       continue
20    continue
      return
      end
      subroutine maktwo(g,ndim,super,ir,ic)
c
      implicit real*8 (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine the two electron integrals are gathered
c.....
      dimension g(*),super(*),ir(ndim),ic(ndim)
      it = 0
      do 40 i=2,ndim
         do 30 j=1,i-1
            do 20 k=2,ndim
               do 10 l=1,k-1
                  it = it + 1
                  ikjl = intpos(ir(i),ic(k),ir(j),ic(l))
                  iljk = intpos(ir(i),ic(l),ir(j),ic(k))
                  gc = super(ikjl)
                  ge = super(iljk)
                  g(it) = gc - ge
10             continue
20          continue
30       continue
40    continue
      return
      end
      subroutine cofac2(ndim,w,x,y,s,irank,det)
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
c.....
c.....in this routine the second order cofactors are calculated.
c.....see: f. prosser and s. hagstrom, int.j.quantum chem. 2, 89 (1968),
c.....and: j. verbeek (thesis, 1990). s is the biorthogonalised overlap
c.....matrix in question. irank is its rank. det is the product of the
c.....non-zero diagonal elements of s. y is the first order adjugate of
c.....s transposed (!). w is the second order adjugate of s. i and j re-
c.....fer to the left hand slater determinant (rows), k and l refer to
c.....the right hand slater determinant (columns). w has the elements of
c.....the second order adjugate in column order. x is a scratch array.
c.....i < j and k < l
c.....
      dimension w((ndim*(ndim-1)/2)**2),x(ndim,ndim),y(ndim,ndim),
     &          s(ndim,ndim)
c.....
      if (ndim.lt.2) return
      if (ndim.eq.2) then
c.....
c.....   the second order cofactor of a two-dimensional matrix is 1.
c.....
         w(1)=1.0d0
         return
      end if
c.....
      if (irank.lt.ndim-2) then
c.....
c.....   nullity.ge.3
c.....   no non-zero second order cofactors, therefore no two-electron
c.....   contribution to the matrix element.
c.....   already handled by matre2: if (nsing.gt.2) return.
c.....
         return
      end if
c.....
      if (irank.lt.ndim-1) then
c.....
c.....   nullity.eq.2
c.....
         k=ndim-1
         l=ndim
c.....
c.....   calculate last column of second order compound matrix of r,
c.....   element ij,kl is in x(i,j).
c.....   calculate last row    of second order compound matrix of l,
c.....   element kl,ij is in x(j,i).
c.....
c.....   i < j and j < k
         do 20 j=2,k-1
            do 10 i=1,j-1
               x(i,j)=(s(i,k)*s(j,l)-s(i,l)*s(j,k))*det
               x(j,i)= s(k,i)*s(l,j)-s(l,i)*s(k,j)
10          continue
20       continue
c.....
         do 30 i=1,k-1
c.....      i < k and j = k
            x(i,k)=(s(i,k)*s(k,l)-s(i,l))*det
            x(k,i)= s(k,i)*s(l,k)-s(l,i)
c.....      i < k and j = l
            x(i,l)=s(i,k)*det
            x(l,i)=s(k,i)
30       continue
c.....
c.....   i = k and j = l
         x(k,l)=1.0d0*det
         x(l,k)=1.0d0
c.....
         it=0
         do 70 l=2,ndim
            do 60 k=1,l-1
               do 50 j=2,ndim
                  do 40 i=1,j-1
                     it=it+1
                     w(it)=x(i,j)*x(l,k)
40                continue
50             continue
60          continue
70       continue
      return
      end if
c.....
      if (irank.eq.ndim-1) then
c.....
c.....   nullity.eq.1
c.....   use factorised cofactor algorithm.
c.....
         call vclr(x,1,ndim*ndim)
c.....
         do 100 i=1,ndim-1
            x(i,i)=1.0d0/s(i,i)
100      continue
         do 120 i=2,ndim-1
            do 110 j=1,i-1
               x(i,j)=x(i,i)*s(i,j)
110         continue
120      continue
c.....
         do 150 j=1,ndim-1
            do 140 i=1,ndim-1
               do 130 k=max(i+1,j),ndim-1
                  x(i,j)=x(i,j)+x(k,j)*s(i,k)
130            continue
140         continue
150      continue
c.....
         it=0
         do 190 l=2,ndim
            do 180 k=1,l-1
               do 170 j=2,ndim
                  do 160 i=1,j-1
                      it=it+1
                      w(it)=x(i,k)*y(l,j)+x(j,l)*y(k,i)-
     &                      x(i,l)*y(k,j)-x(j,k)*y(l,i)
160               continue
170            continue
180         continue
190      continue
         return
      end if
c.....
c.....nullity.eq.0
c.....use jacobi ratio theorem.
c.....
      it=0
      do 230 l=2,ndim
         do 220 k=1,l-1
            do 210 j=2,ndim
               do 200 i=1,j-1
                  it=it+1
                  w(it)=(y(k,i)*y(l,j)-y(l,i)*y(k,j))/det
200            continue
210         continue
220      continue
230   continue
c.....
      return
      end
      subroutine check32(seed,array,n,check,test)
c
c...  calculate and return/test a 32(or 64 for i8) -bits checksum
c
      implicit integer (a-z)
      dimension array(n)
      character*(*) test
c
      cc = seed
      do i=1,n
         cc = ieor(cc,array(i))
      end do
c
      if (test.eq.'test') then
         if (cc.ne.check) call caserr('checksum error')
      else if (test.eq.'calc') then
         check = cc
      else
         call caserr('wrong call to checksum')
      end if
c
      return
      end
      subroutine remco(a)
      character*(*) a
      logical oroot
c      if (.not.oroot()) return
      print *,ipg_nodeid(),': remco : ',a
      call flush(6)
      end
      subroutine remcoi(a,i)
      character*(*) a
      logical oroot
c      if (.not.oroot()) return
      print *,ipg_nodeid(),': remco : ',a,' with ',i
      call flush(6)
      end
      subroutine remcof(a,f)
      real*8 f
      character*(*) a
      print *,ipg_nodeid(),': remco : ',a,' with ',f
      call flush(6)
      end
      subroutine writetxi(fn,q,nelem,mode)
c
      integer q,d
      integer nelem,i
      character*12 fn
      character*1  mode
      dimension q(nelem)
c
      print *,'Writing ',fn
      if (mode.eq.'r') then
         open(unit=31,file=fn,form='formatted')
         rewind(31)
      elseif (mode.eq.'a') then
         open(unit=31,file=fn,form='formatted')
         rewind(31)
10       read(unit=31,fmt='(f50.35)',err=20) d
         goto 10
20       continue
      else
         print *,'Write modus should be a or r'
         return
      endif
      write(31,'(i10)') nelem
      do i=1,nelem
         write(31,'(i50)') q(i)
      enddo
      close(31)
      return
      end
      subroutine readtxi(fn,q)
c
      integer q
      integer nelem,i
      character*12 fn
      dimension q(*)
c
      print *,'Reading ',fn
      open(unit=31,file=fn,form='formatted')
      rewind(31)
      read(31,'(i10)') nelem
      do i=1,nelem
         read(31,'(i50)') q(i)
      enddo
      close(31)
      return
      end
      subroutine writetxt(fn,q,nelem,mode)
c
      real*8 q,d
      integer nelem,i
      character*12 fn
      character*1 mode
      character*75 line
      dimension q(nelem)
c
      print *,'Writing ',fn
      if (mode.eq.'r') then
         open(unit=31,file=fn,form='formatted')
         rewind(31)
      elseif (mode.eq.'a') then
         open(unit=31,file=fn,form='formatted')
         rewind(31)
10       read(unit=31,fmt='(a75)',end=20) line
         goto 10
20       continue
      else
         print *,'Write modus should be a or r'
         return
      endif
      write(31,'(i10)') nelem
      do i=1,nelem
         write(31,'(Z16,f50.35)') q(i),q(i)
      enddo
      close(31)
      return
      end
      subroutine readtxt(fn,q)
c
      real*8 q,d
      integer nelem,i
      character*12 fn
      dimension q(*)
c
      print *,'Reading ',fn
      open(unit=31,file=fn,form='formatted')
      rewind(31)
      read(31,'(i10)') nelem
      do i=1,nelem
         read(31,'(Z16,f50.35)') q(i),d
      enddo
      close(31)
      return
      end 
