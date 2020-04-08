c***
c***   this code for machines with BLAS2 (dgemmv) and 
c***   BLAS3 (dgemm)
c***
      subroutine mxm(a,nar,b,nac,c,nbc)
      implicit REAL  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(nar,nac),b(nac,nbc),c(nar,nbc)
_IF(convex,hp700)
c
c ..  this sequence appears to circumvent the dgemm problems
c ..  in the convex dgemm (at least on cxa)
c
      call vclr(c,1,nar*nbc)
      call mxmb(a,1,nar,b,1,nac,c,1,nar,nar,nac,nbc)
_ELSE
      data xn /'n'/
      call dgemm(xn , xn , nar , nbc , nac 
     &    , 1.0d0 , a , nar , b , nac , 0.0d0 
     &    , c , nar )
_ENDIF
      return
      end
_IF1(i)@process directive('*vdir:')
      subroutine mxmaa(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
      implicit REAL  (a-h,p-w),integer  (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(*),b(*)
      logical done
_IFN1(x)      data xn, xt /'n','t'/
cjvl   clear result if nlink 0 (e.g. for mcscf)
      if (nlink.le.0) then
         ir=1
         do 1001 i=1,nrow
         irr = ir
            do 2002 j=1,ncol
            r(irr) = 0.0d0
 2002    irr = irr + mcolr
 1001 ir = ir + mrowr
      end if
cjvl
      if(ncol.le.0. or. nrow.le.0 .or. nlink.le.0) return
      done = .false.
c
c        if any of the main dimensions is 1 then we have not got a
c        true matrix - matrix multiply. skip to the fortran code.
c
      index = 0
      if ( ncol .gt. 1 )  index = index + 4
      if ( nrow .gt. 1 )  index = index + 2
      if ( nlink .gt. 1 )  index = index + 1
c
      if ( index .lt. 7 )  then
         if ( index .ge. 5 )  then
c
c        index = 6 -- ncol > 1  nrow > 1  nlink = 1
c        index = 5 -- ncol > 1  nrow = 1  nlink > 1
c
            go to 1000
         elseif ( index .eq. 4 )  then
c
c        index = 4 -- ncol > 1  nrow = 1  nlink = 1
c
            ir = 1
            ia = 1
            do 10  i = 1,ncol
               r(ir) = b(1) * a(ia)
               ir = ir + mcolr
               ia = ia + mcola
 10         continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 3 )  then
c
c        index = 3 -- ncol = 1  nrow > 1  nlink > 1
c
            ir = 1
            do 20  j = 1,nrow
               r(ir) = 0.0d0
               ir = ir + mrowr
20          continue
c
            ir = 1
            ib = 1
            do 30  j = 1,nrow
               ibb = ib
               ia = 1
               do 40  k = 1,nlink
                  r(ir) = r(ir) + b(ibb) * a(ia)
                  ibb = ibb + mcolb
                  ia = ia + mrowa
40             continue
               ir = ir + mrowr
               ib = ib + mrowb
30          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 2 )  then
c
c        index = 2 -- ncol = 1  nrow > 1  nlink = 1
c
            ir = 1
            ib = 1
            do 50  i = 1,nrow
               r(ir) = a(1) * b(ib)
               ib = ib + mrowb
               ir = ir + mrowr
50          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 1 )  then
c
c        index = 1 -- ncol = 1  nrow = 1  nlink > 1
c
            ia = 1
            ib = 1
            r(1) = 0.0d0
            do 60  i = 1,nlink
               r(1) = r(1) + a(ia) * b(ib)
               ia = ia + mrowa
               ib = ib + mcolb
60          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 0 )  then
c
c        index = 0 -- ncol = 1  nrow = 1  nlink = 1
c
            r(1) = a(1) * b(1)
            done = .true.
            go to 1000
c
         endif
      endif
c
_IFN(convex)
      index = 1
      if ( mcola .gt. 1 )  index = index + 4
      if ( mcolb .gt. 1 )  index = index + 2
      if ( mcolr .gt. 1 )  index = index +1
      go to (100,200,300,400,500,600,700,800) , index
c
c  r = a * b
c
  100 continue
      call dgemm(xn, xn, ncol, nrow, nlink, 1.0d0
     &     , a, mrowa , b, mrowb, 0.0d0, r, mrowr )
      done = .true.
      go to 1000
c
c  r(t) = a * b
c
  200 continue
      if ( mrowr .eq. 1 )  then
         call dgemm(xt, xt, nrow, ncol, nlink, 1.0d0
     &        , b, mrowb , a, mrowa, 0.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a * b(t)
c
  300 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 )  then
         call dgemm(xn, xt, ncol, nrow, nlink, 1.0d0
     &        , a, mrowa , b, mcolb, 0.0d0, r, mrowr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r(t) = a * b(t)
c
  400 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xt, nrow, ncol, nlink, 1.0d0
     &        , b, mcolb, a, mrowa, 0.0d0, r, mcolr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r = a(t) * b
c
  500 continue
      if ( mrowa .eq. 1 )  then
         call dgemm(xt, xn, ncol, nrow, nlink, 1.0d0
     &        , a, mcola, b, mrowb, 0.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b
c
  600 continue
      if ( mrowa .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xt, xn, nrow, ncol, nlink, 1.0d0
     &        , b, mrowb, a, mcola, 0.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a(t) * b(t)
c
  700 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 )  then
         call dgemm(xt, xt, ncol, nrow, nlink, 1.0d0
     &        , a, mcola, b, mcolb,  0.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b(t)
c
  800 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xn, nrow, ncol, nlink, 1.0d0
     &        , b, mcolb, a, mcola, 0.0d0, r, mcolr )
         done = .true.
      endif
_ENDIF
c
 1000 continue
c
c        check to see if we have already done the calculation in any of
c        the special cases.
c
      if ( done )  return
c
      ir=1
      ib=1
      do 1 j=1,nrow
      irr = ir
      do 44  loop=1,ncol
         r(irr) = 0.0d0
         irr = irr + mcolr
44    continue
      ibb=ib
      ia=1
      icount=1
      do 22 k=1,nlink
      fac=b(ibb)
      if(fac)4,2,4
    4 goto (11,12,13,14),icount
   11 iaa1=ia
      fac1=fac
      icount=2
      goto 2
   12 iaa2=ia
      fac2=fac
      icount=3
      goto 2
   13 iaa3=ia
      fac3=fac
      icount=4
      goto 2
   14 irr=ir
       iaa=ia
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
      do 3 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)+fac3*a(iaa3)+fac*a(iaa)
      irr=irr+mcolr
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
    3  iaa=iaa+mcola
        icount=1
    2   ibb=ibb+mcolb
   22     ia=ia+mrowa
       irr=ir
       goto (88,31,32,33),icount
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   31   do 41 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)
       iaa1=iaa1+mcola
   41  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   32  do 42 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
   42  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   33  do 43 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)+fac3*a(iaa3)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
   43  irr=irr+mcolr
   88 ir=ir+mrowr
    1 ib=ib+mrowb
      return
      end
c
c
c

_IF1(i)@process directive('*vdir:')
      subroutine mxmb(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
      implicit REAL  (a-h,p-w),integer  (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(*),b(*)
      logical done
      data xn, xt /'n','t'/
      if(ncol.le.0. or. nrow.le.0 .or. nlink.le.0) return
      done = .false.
c
c        if any of the main dimensions is 1 then we have not got a
c        true matrix - matrix multiply. skip to the fortran code.
c
      index = 0
      if ( ncol .gt. 1 )  index = index + 4
      if ( nrow .gt. 1 )  index = index + 2
      if ( nlink .gt. 1 )  index = index + 1
c
      if ( index .lt. 7 )  then
         if ( index .ge. 5 )  then
c
c        index = 6 -- ncol > 1  nrow > 1  nlink = 1
c        index = 5 -- ncol > 1  nrow = 1  nlink > 1
c
            go to 1000
         elseif ( index .eq. 4 )  then
c
c        index = 4 -- ncol > 1  nrow = 1  nlink = 1
c
            ir = 1
            ia = 1
            do 10  i = 1,ncol
               r(ir) = r(ir) + b(1) * a(ia)
               ir = ir + mcolr
               ia = ia + mcola
10          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 3 )  then
c
c        index = 3 -- ncol = 1  nrow > 1  nlink > 1
c
            ir = 1
            ib = 1
            do 20  j = 1,nrow
               ibb = ib
               ia = 1
               do 30  k = 1,nlink
                  r(ir) = r(ir) + b(ibb) * a(ia)
                  ibb = ibb + mcolb
                  ia = ia + mrowa
30             continue
               ir = ir + mrowr
               ib = ib + mrowb
20          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 2 )  then
c
c        index = 2 -- ncol = 1  nrow > 1  nlink = 1
c
            ir = 1
            ib = 1
            do 40  i = 1,nrow
               r(ir) = r(ir) + a(1) * b(ib)
               ib = ib + mrowb
               ir = ir + mrowr
40          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 1 )  then
c
c        index = 1 -- ncol = 1  nrow = 1  nlink > 1
c
            ia = 1
            ib = 1
            do 50  i = 1,nlink
               r(1) = r(1) + a(ia) * b(ib)
               ia = ia + mrowa
               ib = ib + mcolb
50          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 0 )  then
c
c        index = 0 -- ncol = 1  nrow = 1  nlink = 1
c
            r(1) = r(1) + a(1) * b(1)
            done = .true.
            go to 1000
c
         endif
      endif
c
      index = 1
      if ( mcola .gt. 1 )  index = index + 4
      if ( mcolb .gt. 1 )  index = index + 2
      if ( mcolr .gt. 1 )  index = index +1
      go to (100,200,300,400,500,600,700,800) , index
c
c  r = a * b
c
  100 continue
      call dgemm(xn, xn, ncol, nrow, nlink, 1.0d0
     &     , a, mrowa, b, mrowb, 1.0d0, r, mrowr )
      done = .true.
      go to 1000
c
c  r(t) = a * b
c
  200 continue
      if ( mrowr .eq. 1 )  then
         call dgemm(xt, xt, nrow, ncol, nlink, 1.0d0
     &        , b, mrowb, a, mrowa, 1.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a * b(t)
c
  300 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 )  then
         call dgemm(xn, xt, ncol, nrow, nlink, 1.0d0
     &        , a, mrowa, b, mcolb, 1.0d0, r, mrowr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r(t) = a * b(t)
c
  400 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xt, nrow, ncol, nlink, 1.0d0
     &        , b, mcolb, a, mrowa, 1.0d0, r, mcolr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r = a(t) * b
c
  500 continue
      if ( mrowa .eq. 1 )  then
         call dgemm(xt, xn, ncol, nrow, nlink, 1.0d0
     &        , a, mcola, b, mrowb, 1.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b
c
  600 continue
      if ( mrowa .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xt, xn, nrow, ncol, nlink, 1.0d0
     &        , b, mrowb, a, mcola, 1.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a(t) * b(t)
c
  700 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 )  then
         call dgemm(xt, xt, ncol, nrow, nlink, 1.0d0
     &        , a, mcola, b, mcolb,  1.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b(t)
c
  800 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xn, nrow, ncol, nlink, 1.0d0
     &        , b, mcolb, a, mcola, 1.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
 1000 continue
c
c        check to see if we have already done the calculation in any of
c        the special cases.
c
      if ( done )  return
c
      ir=1
      ib=1
      do 1 j=1,nrow
      ibb=ib
      ia=1
      icount=1
      do 22 k=1,nlink
      fac=b(ibb)
      if(fac)4,2,4
    4 goto (11,12,13,14),icount
   11 iaa1=ia
      fac1=fac
      icount=2
      goto 2
   12 iaa2=ia
      fac2=fac
      icount=3
      goto 2
   13 iaa3=ia
      fac3=fac
      icount=4
      goto 2
   14 irr=ir
       iaa=ia
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
      do 3 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)+fac3*a(iaa3)+fac*a(iaa)
      irr=irr+mcolr
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
    3  iaa=iaa+mcola
        icount=1
    2   ibb=ibb+mcolb
   22     ia=ia+mrowa
       irr=ir
       goto (88,31,32,33),icount
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   31   do 41 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)
       iaa1=iaa1+mcola
   41  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   32  do 42 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
   42  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   33  do 43 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)+fac3*a(iaa3)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
   43  irr=irr+mcolr
   88 ir=ir+mrowr
    1 ib=ib+mrowb
      return
      end
_IF1(i)@process directive('*vdir:')
      subroutine mxmbn(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
      implicit REAL  (a-h,p-w),integer  (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(*),b(*)
      logical done
      data xn, xt /'n','t'/
      if(ncol.le.0. or. nrow.le.0 .or. nlink.le.0) return
      done = .false.
c
c        if any of the main dimensions is 1 then we have not got a
c        true matrix - matrix multiply. skip to the fortran code.
c
      index = 0
      if ( ncol .gt. 1 )  index = index + 4
      if ( nrow .gt. 1 )  index = index + 2
      if ( nlink .gt. 1 )  index = index + 1
c
      if ( index .lt. 7 )  then
         if ( index .ge. 5 )  then
c
c        index = 6 -- ncol > 1  nrow > 1  nlink = 1
c        index = 5 -- ncol > 1  nrow = 1  nlink > 1
c
            go to 1000
         elseif ( index .eq. 4 )  then
c
c        index = 4 -- ncol > 1  nrow = 1  nlink = 1
c
            ir = 1
            ia = 1
            do 10  i = 1,ncol
               r(ir) = r(ir) - b(1) * a(ia)
               ir = ir + mcolr
               ia = ia + mcola
10          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 3 )  then
c
c        index = 3 -- ncol = 1  nrow > 1  nlink > 1
c
            ir = 1
            ib = 1
            do 20  j = 1,nrow
               ibb = ib
               ia = 1
               do 30  k = 1,nlink
                  r(ir) = r(ir) - b(ibb) * a(ia)
                  ibb = ibb + mcolb
                  ia = ia + mrowa
30             continue
               ir = ir + mrowr
               ib = ib + mrowb
20          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 2 )  then
c
c        index = 2 -- ncol = 1  nrow > 1  nlink = 1
c
            ir = 1
            ib = 1
            do 40  i = 1,nrow
               r(ir) = r(ir) - a(1) * b(ib)
               ib = ib + mrowb
               ir = ir + mrowr
40          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 1 )  then
c
c        index = 1 -- ncol = 1  nrow = 1  nlink > 1
c
            ia = 1
            ib = 1
            do 50  i = 1,nlink
               r(1) = r(1) - a(ia) * b(ib)
               ia = ia + mrowa
               ib = ib + mcolb
50          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 0 )  then
c
c        index = 0 -- ncol = 1  nrow = 1  nlink = 1
c
            r(1) = r(1) - a(1) * b(1)
            done = .true.
            go to 1000
c
         endif
      endif
c
      index = 1
      if ( mcola .gt. 1 )  index = index + 4
      if ( mcolb .gt. 1 )  index = index + 2
      if ( mcolr .gt. 1 )  index = index +1
      go to (100,200,300,400,500,600,700,800) , index
c
c  r = a * b
c
  100 continue
      call dgemm(xn, xn, ncol, nrow, nlink, -1.0d0
     &     , a, mrowa, b, mrowb, 1.0d0, r, mrowr )
      done = .true.
      go to 1000
c
c  r(t) = a * b
c
  200 continue
      if ( mrowr .eq. 1 )  then
         call dgemm(xt, xt, nrow, ncol, nlink, -1.0d0
     &        , b, mrowb, a, mrowa, 1.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a * b(t)
c
  300 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 )  then
         call dgemm(xn, xt, ncol, nrow, nlink, -1.0d0
     &        , a, mrowa, b, mcolb, 1.0d0, r, mrowr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r(t) = a * b(t)
c
  400 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xt, nrow, ncol, nlink, -1.0d0
     &        , b, mcolb, a, mrowa, 1.0d0, r, mcolr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r = a(t) * b
c
  500 continue
      if ( mrowa .eq. 1 )  then
         call dgemm(xt, xn, ncol, nrow, nlink, -1.0d0
     &        , a, mcola, b, mrowb, 1.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b
c
  600 continue
      if ( mrowa .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xt, xn, nrow, ncol, nlink, -1.0d0
     &        , b, mrowb, a, mcola, 1.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a(t) * b(t)
c
  700 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 )  then
         call dgemm(xt, xt, ncol, nrow, nlink, -1.0d0
     &        , a, mcola, b, mcolb,  1.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b(t)
c
  800 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xn, nrow, ncol, nlink, -1.0d0
     &        , b, mcolb, a, mcola, 1.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
 1000 continue
c
c        check to see if we have already done the calculation in any of
c        the special cases.
c
      if ( done )  return
c
      ir=1
      ib=1
      do 1 j=1,nrow
      ibb=ib
      ia=1
      icount=1
      do 22 k=1,nlink
      fac=b(ibb)
      if(fac)4,2,4
    4 goto (11,12,13,14),icount
   11 iaa1=ia
      fac1=fac
      icount=2
      goto 2
   12 iaa2=ia
      fac2=fac
      icount=3
      goto 2
   13 iaa3=ia
      fac3=fac
      icount=4
      goto 2
   14 irr=ir
       iaa=ia
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
      do 3 loop=1,ncol
       r(irr)=r(irr)-fac1*a(iaa1)-fac2*a(iaa2)-fac3*a(iaa3)-fac*a(iaa)
      irr=irr+mcolr
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
    3  iaa=iaa+mcola
        icount=1
    2   ibb=ibb+mcolb
   22     ia=ia+mrowa
       irr=ir
       goto (88,31,32,33),icount
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   31   do 41 loop=1,ncol
       r(irr)=r(irr)-fac1*a(iaa1)
       iaa1=iaa1+mcola
   41  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   32  do 42 loop=1,ncol
       r(irr)=r(irr)-fac1*a(iaa1)-fac2*a(iaa2)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
   42  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   33  do 43 loop=1,ncol
       r(irr)=r(irr)-fac1*a(iaa1)-fac2*a(iaa2)-fac3*a(iaa3)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
   43  irr=irr+mcolr
   88 ir=ir+mrowr
    1 ib=ib+mrowb
      return
      end
_IF1(i)@process directive('*vdir:')
      subroutine mxmd(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
      implicit REAL  (a-h,p-w),integer  (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(*),b(*)
      logical done
      data xn, xt /'n','t'/
      if(ncol.le.0. or. nrow.le.0 .or. nlink.le.0) return
      done = .false.
c
c        if any of the main dimensions is 1 then we have not got a
c        true matrix - matrix multiply. skip to the fortran code.
c
      index = 0
      if ( ncol .gt. 1 )  index = index + 4
      if ( nrow .gt. 1 )  index = index + 2
      if ( nlink .gt. 1 )  index = index + 1
c
      if ( index .lt. 7 )  then
         if ( index .ge. 5 )  then
c
c        index = 6 -- ncol > 1  nrow > 1  nlink = 1
c        index = 5 -- ncol > 1  nrow = 1  nlink > 1
c
            go to 1000
         elseif ( index .eq. 4 )  then
c
c        index = 4 -- ncol > 1  nrow = 1  nlink = 1
c
            ir = 1
            ia = 1
            do 10  i = 1,ncol
               r(ir) = b(1) * a(ia)
               ir = ir + mcolr
               ia = ia + mcola
10          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 3 )  then
c
c        index = 3 -- ncol = 1  nrow > 1  nlink > 1
c
            ir = 1
            do 20  j = 1,nrow
               r(ir) = 0.0d0
               ir = ir + mrowr
20          continue
c
            ir = 1
            ib = 1
            do 30  j = 1,nrow
               ibb = ib
               ia = 1
               do 40  k = 1,nlink
                  r(ir) = r(ir) + b(ibb) * a(ia)
                  ibb = ibb + mcolb
                  ia = ia + mrowa
40             continue
               ir = ir + mrowr
               ib = ib + mrowb
30          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 2 )  then
c
c        index = 2 -- ncol = 1  nrow > 1  nlink = 1
c
            ir = 1
            ib = 1
            do 50  i = 1,nrow
               r(ir) = a(1) * b(ib)
               ib = ib + mrowb
               ir = ir + mrowr
50          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 1 )  then
c
c        index = 1 -- ncol = 1  nrow = 1  nlink > 1
c
            ia = 1
            ib = 1
            r(1) = 0.0d0
            do 60  i = 1,nlink
               r(1) = r(1) + a(ia) * b(ib)
               ia = ia + mrowa
               ib = ib + mcolb
60          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 0 )  then
c
c        index = 0 -- ncol = 1  nrow = 1  nlink = 1
c
            r(1) = a(1) * b(1)
            done = .true.
            go to 1000
c
         endif
      endif
c
      index = 1
      if ( mcola .gt. 1 )  index = index + 4
      if ( mcolb .gt. 1 )  index = index + 2
      if ( mcolr .gt. 1 )  index = index +1
      go to (100,200,300,400,500,600,700,800) , index
c
c  r = a * b
c
  100 continue
      call dgemm(xn, xn, ncol, nrow, nlink, 1.0d0
     &     , a, mrowa, b, mrowb, 0.0d0, r, mrowr )
      done = .true.
      go to 1000
c
c  r(t) = a * b
c
  200 continue
      if ( mrowr .eq. 1 )  then
         call dgemm(xt, xt, nrow, ncol, nlink, 1.0d0
     &        , b, mrowb, a, mrowa, 0.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a * b(t)
c
  300 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 )  then
         call dgemm(xn, xt, ncol, nrow, nlink, 1.0d0
     &        , a, mrowa, b, mcolb, 0.0d0, r, mrowr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r(t) = a * b(t)
c
  400 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xt, nrow, ncol, nlink, 1.0d0
     &        , b, mcolb, a, mrowa, 0.0d0, r, mcolr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r = a(t) * b
c
  500 continue
      if ( mrowa .eq. 1 )  then
         call dgemm(xt, xn, ncol, nrow, nlink, 1.0d0
     &        , a, mcola, b, mrowb, 0.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b
c
  600 continue
      if ( mrowa .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xt, xn, nrow, ncol, nlink, 1.0d0
     &        , b, mrowb, a, mcola, 0.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a(t) * b(t)
c
  700 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 )  then
         call dgemm(xt, xt, ncol, nrow, nlink, 1.0d0
     &        , a, mcola,  b, mcolb,  0.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b(t)
c
  800 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xn, nrow, ncol , nlink, 1.0d0
     &        , b, mcolb, a, mcola, 0.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
 1000 continue
c
c        check to see if we have already done the calculation in any of
c        the special cases.
c
      if ( done )  return
c
      do 70  i = 1,nrow * ncol
         r(i) = 0.0d0
70    continue
      ir=1
      ib=1
      do 1 j=1,nrow
      ibb=ib
      ia=1
      icount=1
      do 22 k=1,nlink
      fac=b(ibb)
      if(fac)4,2,4
    4 goto (11,12,13,14),icount
   11 iaa1=ia
      fac1=fac
      icount=2
      goto 2
   12 iaa2=ia
      fac2=fac
      icount=3
      goto 2
   13 iaa3=ia
      fac3=fac
      icount=4
      goto 2
   14 irr=ir
       iaa=ia
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
      do 3 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)+fac3*a(iaa3)+fac*a(iaa)
      irr=irr+mcolr
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
    3  iaa=iaa+mcola
        icount=1
    2   ibb=ibb+mcolb
   22     ia=ia+mrowa
       irr=ir
       goto (88,31,32,33),icount
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   31   do 41 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)
       iaa1=iaa1+mcola
   41  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   32  do 42 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
   42  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   33  do 43 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)+fac3*a(iaa3)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
   43  irr=irr+mcolr
   88 ir=ir+mrowr
    1 ib=ib+mrowb
      return
      end
_IF1(i)@process directive('*vdir:')
      subroutine mxmg(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
      implicit REAL  (a-h,p-w),integer  (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(*),b(*)
      logical done
      data xn, xt /'n','t'/
      if(ncol.le.0. or. nrow.le.0 .or. nlink.le.0) return
      done = .false.
c
c        if any of the main dimensions is 1 then we have not got a
c        true matrix - matrix multiply. skip to the fortran code.
c
      index = 0
      if ( ncol .gt. 1 )  index = index + 4
      if ( nrow .gt. 1 )  index = index + 2
      if ( nlink .gt. 1 )  index = index + 1
c
      if ( index .lt. 7 )  then
         if ( index .ge. 5 )  then
c
c        index = 6 -- ncol > 1  nrow > 1  nlink = 1
c        index = 5 -- ncol > 1  nrow = 1  nlink > 1
c
            go to 1000
         elseif ( index .eq. 4 )  then
c
c        index = 4 -- ncol > 1  nrow = 1  nlink = 1
c
            ir = 1
            ia = 1
            do 10  i = 1,ncol
               r(ir) = b(1) * a(ia)
               ir = ir + mcolr
               ia = ia + mcola
10          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 3 )  then
c
c        index = 3 -- ncol = 1  nrow > 1  nlink > 1
c
            ir = 1
            do 20  j = 1,nrow
               r(ir) = 0.0d0
               ir = ir + mrowr
20          continue
c
            ir = 1
            ib = 1
            do 30  j = 1,nrow
               ibb = ib
               ia = 1
               do 40  k = 1,nlink
                  r(ir) = r(ir) + b(ibb) * a(ia)
                  ibb = ibb + mcolb
                  ia = ia + mrowa
40             continue
               ir = ir + mrowr
               ib = ib + mrowb
30          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 2 )  then
c
c        index = 2 -- ncol = 1  nrow > 1  nlink = 1
c
            ir = 1
            ib = 1
            do 50  i = 1,nrow
               r(ir) = a(1) * b(ib)
               ib = ib + mrowb
               ir = ir + mrowr
50          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 1 )  then
c
c        index = 1 -- ncol = 1  nrow = 1  nlink > 1
c
            ia = 1
            ib = 1
            r(1) = 0.0d0
            do 60  i = 1,nlink
               r(1) = r(1) + a(ia) * b(ib)
               ia = ia + mrowa
               ib = ib + mcolb
60          continue
            done = .true.
            go to 1000
c
         elseif ( index .eq. 0 )  then
c
c        index = 0 -- ncol = 1  nrow = 1  nlink = 1
c
            r(1) = a(1) * b(1)
            done = .true.
            go to 1000
c
         endif
      endif
c
      index = 1
      if ( mcola .gt. 1 )  index = index + 4
      if ( mcolb .gt. 1 )  index = index + 2
      if ( mcolr .gt. 1 )  index = index +1
      go to (100,200,300,400,500,600,700,800) , index
c
c  r = a * b
c
  100 continue
      call dgemm(xn, xn, ncol, nrow, nlink, 1.0d0
     &     , a, mrowa, b, mrowb, 0.0d0, r, mrowr )
      done = .true.
      go to 1000
c
c  r(t) = a * b
c
  200 continue
      if ( mrowr .eq. 1 )  then
         call dgemm(xt, xt, nrow, ncol, nlink, 1.0d0
     &        , b, mrowb, a, mrowa, 0.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a * b(t)
c
  300 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 )  then
         call dgemm(xn, xt, ncol, nrow, nlink, 1.0d0
     &        , a, mrowa, b, mcolb, 0.0d0, r, mrowr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r(t) = a * b(t)
c
  400 continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if ( mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xt, nrow, ncol, nlink, 1.0d0
     &        , b, mcolb, a, mrowa, 0.0d0, r, mcolr )
         done = .true.
      endif
_ENDIF
      go to 1000
c
c  r = a(t) * b
c
  500 continue
      if ( mrowa .eq. 1 )  then
         call dgemm(xt, xn, ncol, nrow, nlink, 1.0d0
     &        , a, mcola, b, mrowb, 0.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b
c
  600 continue
      if ( mrowa .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xt, xn, nrow, ncol, nlink, 1.0d0
     &        , b, mrowb, a, mcola, 0.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
c  r = a(t) * b(t)
c
  700 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 )  then
         call dgemm(xt, xt, ncol, nrow, nlink, 1.0d0
     &        , a, mcola, b, mcolb,  0.0d0, r, mrowr )
         done = .true.
      endif
      go to 1000
c
c  r(t) = a(t) * b(t)
c
  800 continue
      if ( mrowa .eq. 1 .and. mrowb .eq. 1 .and. mrowr .eq. 1 )  then
         call dgemm(xn, xn, nrow, ncol, nlink, 1.0d0
     &        , b, mcolb, a, mcola, 0.0d0, r, mcolr )
         done = .true.
      endif
      go to 1000
c
 1000 continue
c
c        check to see if we have already done the calculation in any of
c        the special cases.
c
      if ( done )  return
c
      ir=1
      ib=1
      do 1 j=1,nrow
      irr = ir
      do 80  loop=1,ncol
         r(irr) = 0.0d0
         irr = irr + mcolr
80    continue
      ibb=ib
      ia=1
      icount=1
      do 22 k=1,nlink
      fac=b(ibb)
      if(fac)4,2,4
    4 goto (11,12,13,14),icount
   11 iaa1=ia
      fac1=fac
      icount=2
      goto 2
   12 iaa2=ia
      fac2=fac
      icount=3
      goto 2
   13 iaa3=ia
      fac3=fac
      icount=4
      goto 2
   14 irr=ir
       iaa=ia
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
      do 3 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)+fac3*a(iaa3)+fac*a(iaa)
      irr=irr+mcolr
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
    3  iaa=iaa+mcola
        icount=1
    2   ibb=ibb+mcolb
   22     ia=ia+mrowa
       irr=ir
       goto (88,31,32,33),icount
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   31   do 41 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)
       iaa1=iaa1+mcola
   41  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   32  do 42 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
   42  irr=irr+mcolr
       goto 88
_IF1(x)c$dir no_recurrence
_IF1(ath)cvd$  nodepchk
_IF1(i)c*vdir: ignore recrdeps
   33  do 43 loop=1,ncol
       r(irr)=r(irr)+fac1*a(iaa1)+fac2*a(iaa2)+fac3*a(iaa3)
       iaa1=iaa1+mcola
       iaa2=iaa2+mcola
       iaa3=iaa3+mcola
   43  irr=irr+mcolr
   88 ir=ir+mrowr
    1 ib=ib+mrowb
      return
      end
cccccccccc_IFN(hp700)
      subroutine mxmm(a,mrowa,b,r,mrowr,ncol,nrow)
      implicit REAL  (a-h,p-w),integer    (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension a(mrowa,*),b(*),r(mrowr,*)
      data xn /'n'/
      m=1
_IF(hp700)
      do 1 i=1,nrow
         do 2 j=1,i
            top=b(m)
            if(top)4,2,4
 4          call vec_$dmult_add(r(1,i),a(1,j),ncol,top,r(1,i))
 2          m=m+1
_ELSE
      lda = mrowa
      if ( nrow .eq. 1 )  lda = ncol
      do 1 i=1,nrow
         call dgemv(xn , ncol , i , 1.0d0 , a , lda , b(m) , 1
     *                    , 1.0d0 , r(1,i) , 1 )
         m = m + i
_ENDIF
    1 continue
      return
      end
_IF(hp700,linux)
c
c  as libvec.a BLAS on hp misses dgemv, these mappings are 
c  used (temporary..)
c
c  problem with mult2 using -lblas on pentium II
      subroutine mxmtt(a,aa,mcol,mrow,ne,nf,b)
      implicit REAL  (a-h,o-z),integer (i-n)
      dimension aa(*),a(*),b(*)
      m=0
      n=0
      do 1 loop=1,ne
      call mxmb(a,mcol,mrow,aa(n+1),mrow,1,b(m+1),1,1,loop,nf,1)
      m=m+loop
1     n=n+mcol
       return
      end
      subroutine mxmt(q,mcol,mrow,f,r,na,nbas)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),f(*),r(*)
      ll=1
      mm=1
      do 1 loop=1,nbas
      call mxmb(q,mcol,mrow,f(mm),1,1,r(ll),1,1,na,loop,1)
      mm=mm+loop
1     ll=ll+na
      return
      end
_ELSE
      subroutine mxmt(q,mcol,mrow,f,r,na,nbas)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension q(*),f(*),r(*)
      data xn, xt /'n','t'/
      ll=1
      mm=1
      if ( mcol .eq. 1 )  then
         lda = mrow
         if ( nbas .eq. 1 )  lda = na
         do 10  loop=1,nbas
            call dgemv(xn, na, loop, 1.0d0, q, lda, f(mm), 1
     *                      , 1.0d0, r(ll), 1 )
            mm = mm + loop
            ll = ll + na
10       continue
      elseif ( mrow .eq. 1 )  then
         lda = mcol
         if ( na .eq. 1 )  lda = nbas
         do 20  loop = 1,nbas
            call dgemv(xt, loop, na, 1.0d0, q, lda, f(mm), 1
     *                      , 1.0d0, r(ll), 1 )
            mm = mm + loop
            ll = ll + na
20       continue
      else
         do 30   loop=1,nbas
            call mxmb(q,mcol,mrow,f(mm),1,1,r(ll),1,1,na,loop,1)
            mm=mm+loop
            ll=ll+na
30       continue
      endif
      return
      end
      subroutine mxmtt(a,aa,mcol,mrow,ne,nf,b)
      implicit REAL  (a-h,o-z),integer (i-n)
      character *1 xn, xt
      dimension aa(*),a(*),b(*)
      data xn, xt /'n','t'/
      m=1
      n=1
      if ( mcol .eq. 1 )  then
         lda = mrow
         if ( nf .eq. 1 )  lda = ne
         do 10  loop = 1,ne
            call dgemv(xn, loop, nf, 1.0d0, a, lda, aa(n), mrow
     *                      , 1.0d0, b(m), 1 )
            m = m + loop
            n = n + mcol
10       continue
      elseif ( mrow .eq. 1 )  then
         lda = mcol
         do 20  loop = 1,ne
            if ( loop .eq. 1 )  lda = nf
            call dgemv(xt, nf, loop, 1.0d0, a, lda, aa(n), mrow
     *                      , 1.0d0, b(m), 1 )
            m = m + loop
            n = n + mcol
20       continue
      else
         do 30  loop=1,ne
            call mxmb(a,mcol,mrow,aa(n),mrow,1,b(m),1,1,loop,nf,1)
            m=m+loop
            n=n+mcol
30       continue
      endif
       return
      end
_ENDIF
      subroutine mxms (a,mcola,mrowa,b,mcolb,mrowb,
     >     r,mcolr,mrowr, ncol,nlink,nrow)
      implicit REAL  (a-h,o-z)
      dimension a(*),b(*),r(*)
      call mxmaa (a,mcola,mrowa,b,mcolb,mrowb,
     >      r,mcolr,mrowr, ncol,nlink,nrow )
      return
      end

      subroutine mxma(a,mcola,mrowa,b,mcolb,mrowb,
     *r,mcolr,mrowr,ncol,nlink,nrow)
      implicit REAL  (a-h,p-w),integer  (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension r(*),a(*),b(*)
      logical done
      data xn,xt / 'n','t' /
      if (ncol.le.0 .or. nrow.le.0 .or. nlink.le.0) return

      done = .false.
c
c        if any of the main dimensions is 1 then we have not got a
c        true matrix - matrix multiply. skip to the fortran code.
c
      index = 0
      if (ncol.gt.1) index = index + 4
      if (nrow.gt.1) index = index + 2
      if (nlink.gt.1) index = index + 1
c
      if (index.lt.7) then
c
c        index = 6 -- ncol > 1  nrow > 1  nlink = 1
c        index = 5 -- ncol > 1  nrow = 1  nlink > 1
c
         if (index.ge.5) go to 160
         if (index.eq.4) then
c
c        index = 4 -- ncol > 1  nrow = 1  nlink = 1
c
            ir = 1
            ia = 1
            do 20 i = 1 , ncol
               r(ir) = b(1)*a(ia)
               ir = ir + mcolr
               ia = ia + mcola
 20         continue
            done = .true.
            go to 160
c
         else if (index.eq.3) then
c
c        index = 3 -- ncol = 1  nrow > 1  nlink > 1
c
            ir = 1
            do 30 j = 1 , nrow
               r(ir) = 0.0d0
               ir = ir + mrowr
 30         continue
c
            ir = 1
            ib = 1
            do 50 j = 1 , nrow
               ibb = ib
               ia = 1
               do 40 k = 1 , nlink
                  r(ir) = r(ir) + b(ibb)*a(ia)
                  ibb = ibb + mcolb
                  ia = ia + mrowa
 40            continue
               ir = ir + mrowr
               ib = ib + mrowb
 50         continue
            done = .true.
            go to 160
c
         else if (index.eq.2) then
c
c        index = 2 -- ncol = 1  nrow > 1  nlink = 1
c
            ir = 1
            ib = 1
            do 60 i = 1 , nrow
               r(ir) = a(1)*b(ib)
               ib = ib + mrowb
               ir = ir + mrowr
 60         continue
            done = .true.
            go to 160
c
         else if (index.eq.1) then
c
c        index = 1 -- ncol = 1  nrow = 1  nlink > 1
c
            ia = 1
            ib = 1
            r(1) = 0.0d0
            do 70 i = 1 , nlink
               r(1) = r(1) + a(ia)*b(ib)
               ia = ia + mrowa
               ib = ib + mcolb
 70         continue
            done = .true.
            go to 160
c
         else if (index.eq.0) then
c
c        index = 0 -- ncol = 1  nrow = 1  nlink = 1
c
            r(1) = a(1)*b(1)
            done = .true.
            go to 160
c
         end if
      end if
c
      index = 1
      if (mcola.gt.1) index = index + 4
      if (mcolb.gt.1) index = index + 2
      if (mcolr.gt.1) index = index + 1
      go to (80,90,100,110,120,130,140,150) , index
c
c  r = a * b
c
 80   call dgemm(xn,xn,ncol,nrow,nlink,1.0d0
     +  ,a,mrowa,b,mrowb,0.0d0,r,mrowr)
      done = .true.
      go to 160
c
c  r(t) = a * b
c
 90   if (mrowr.eq.1) then
         call dgemm(xt,xt,nrow,ncol,nlink,1.0d0
     +   ,b,mrowb,a,mrowa,0.0d0,r,mcolr)
         done = .true.
      end if
      go to 160
c
c  r = a * b(t)
c
 100  continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if (mrowb.eq.1) then
         call dgemm(xn,xt,ncol,nrow,nlink,1.0d0
     +   ,a,mrowa,b,mcolb,0.0d0,r,mrowr)
         done = .true.
      end if
_ENDIF
      go to 160
c
c  r(t) = a * b(t)
c
 110  continue
_IF(mips2)
*     dgemm errors under mips2 on r4400
*     bypass dgmemm call
_ELSE
      if (mrowb.eq.1 .and. mrowr.eq.1) then
         call dgemm(xn,xt,nrow,ncol,nlink,1.0d0
     +  ,b,mcolb,a,mrowa,0.0d0,r,mcolr)
         done = .true.
      end if
_ENDIF
      go to 160
c
c  r = a(t) * b
c
 120  if (mrowa.eq.1) then
         call dgemm(xt,xn,ncol,nrow,nlink,1.0d0
     +   ,a,mcola,b,mrowb,0.0d0,r,mrowr)
         done = .true.
      end if
      go to 160
c
c  r(t) = a(t) * b
c
 130  if (mrowa.eq.1 .and. mrowr.eq.1) then
         call dgemm(xt,xn,nrow,ncol,nlink,1.0d0
     +  ,b,mrowb,a,mcola,0.0d0,r,mcolr)
         done = .true.
      end if
      go to 160
c
c  r = a(t) * b(t)
c
 140  if (mrowa.eq.1 .and. mrowb.eq.1) then
         call dgemm(xt,xt,ncol,nrow,nlink,1.0d0
     +  ,a,mcola,b,mcolb,0.0d0,r,mrowr)
         done = .true.
      end if
      go to 160
c
c  r(t) = a(t) * b(t)
c
 150  if (mrowa.eq.1 .and. mrowb.eq.1 .and. mrowr.eq.1) then
         call dgemm(xn,xn,nrow,ncol,nlink,1.0d0
     +  ,b,mcolb,a,mcola,0.0d0,r,mcolr)
         done = .true.
      end if
c
c
c        check to see if we have already done the calculation in any of
c        the special cases.
c
 160  if (done) return
c
      ir = 1
      ib = 1
      do 310 j = 1 , nrow
         irr = ir
         do 170 loop = 1 , ncol
            r(irr) = 0.0d0
            irr = irr + mcolr
 170     continue
         ibb = ib
         ia = 1
         icount = 1
         do 230 k = 1 , nlink
            fac = b(ibb)
            if (fac.ne.0) then
               go to (180,190,200,210) , icount
            else
               ibb = ibb + mcolb
               ia = ia + mrowa
               go to 230
            end if
 180        iaa1 = ia
            fac1 = fac
            icount = 2
            ibb = ibb + mcolb
            ia = ia + mrowa
            go to 230
 190        iaa2 = ia
            fac2 = fac
            icount = 3
            ibb = ibb + mcolb
            ia = ia + mrowa
            go to 230
 200        iaa3 = ia
            fac3 = fac
            icount = 4
            ibb = ibb + mcolb
            ia = ia + mrowa
            go to 230
 210        irr = ir
            iaa = ia
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
            do 220 loop = 1 , ncol
               r(irr) = r(irr) + fac1*a(iaa1) + fac2*a(iaa2)
     +                  + fac3*a(iaa3) + fac*a(iaa)
               irr = irr + mcolr
               iaa1 = iaa1 + mcola
               iaa2 = iaa2 + mcola
               iaa3 = iaa3 + mcola
               iaa = iaa + mcola
 220        continue
            icount = 1
            ibb = ibb + mcolb
            ia = ia + mrowa
 230     continue
         irr = ir
         go to (300,240,260,280) , icount
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 240     do 250 loop = 1 , ncol
            r(irr) = r(irr) + fac1*a(iaa1)
            iaa1 = iaa1 + mcola
            irr = irr + mcolr
 250     continue
         go to 300
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 260     do 270 loop = 1 , ncol
            r(irr) = r(irr) + fac1*a(iaa1) + fac2*a(iaa2)
            iaa1 = iaa1 + mcola
            iaa2 = iaa2 + mcola
            irr = irr + mcolr
 270     continue
         go to 300
_IF1(ct)cdir$ ivdep
_IF1(a)cvd$  nodepck
_IF1(x)c$dir no_recurrence
 280     do 290 loop = 1 , ncol
            r(irr) = r(irr) + fac1*a(iaa1) + fac2*a(iaa2) + fac3*a(iaa3)
            iaa1 = iaa1 + mcola
            iaa2 = iaa2 + mcola
            iaa3 = iaa3 + mcola
            irr = irr + mcolr
 290     continue
 300     ir = ir + mrowr
         ib = ib + mrowb
 310  continue
      return
      end
_IF(hp800)
***** this piece of stupidity is caused by the loader under HPUX 10.20
***** not picking up daxpy from the BLAS library ... it does OK
***** in the matrix benchmarks .. must be loader problem ...
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
C$Id: mxm_blas3.m,v 1.8 2006-01-13 17:54:55 jmht Exp $
      REAL dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
_ENDIF
