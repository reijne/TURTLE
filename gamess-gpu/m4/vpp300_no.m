      function pack2(nw,m)
      implicit real*8  (a-h,p-w),integer (i-n),logical  (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c these version for 32-bit integer machines
c
      dimension  i4(2)
      equivalence (i4(1),pack)
      i4(1)=nw
      i4(2)=m
      pack2 = pack
      return
      end
c======================================================================
      function dor(a,b)
c
c...   64-bit and function
c...   a and b are real*8 in outside world
c
      real*8  rr, dor
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
      data r/2*0/
c
      r(1) = ior(a(1),b(1))
      r(2) = ior(a(2),b(2))
      dor = rr
c
      return
      end
c======================================================================
      function dxor(a,b)
c
c...   64-bit xor function
c...   a and b are real*8 in outside world
c
      real*8  rr, dxor
      integer*4 a(2),b(2),r(2)
      equivalence (rr,r)
      data r/2*0/
c
      r(1) = ieor(a(1),b(1))
      r(2) = ieor(a(2),b(2))
      dxor = rr
c
      return
      end
c======================================================================
      function shift(x,n)
c
c...   real*8 shift (with wraparound)
c...   n <0 v >64 shifts aborts
c...   lshft or rshft by 32 does not empty word !!
c
      implicit real*8  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shift called wrongly')
       call caserr('invalid argument in shiftc')
      endif
c

      a = x
      if (n.eq.0.or.n.eq.64) then
         b = a
      else if (n.lt.32) then
         ib(1) = ior(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = ior(ishft(ia(2),n),ishft(ia(1),n-32))
      else if (n.eq.32) then
         ib(1) = ia(2)
         ib(2) = ia(1)
      else
         ib(1) = ior(ishft(ia(2),n-32),ishft(ia(1),n-64))
         ib(2) = ior(ishft(ia(1),n-32),ishft(ia(2),n-64))
      end if
c
      shift = b
c
       return
      end
c======================================================================
      function shiftl(x,n)
c
c...   real*8  left shift (without wraparound)
c...   n <0 v >64 shifts aborts
c...   shift by 64 bits must clear word
c...   lshft or rshft by 32 does not empty word !!
c
      implicit real*8  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftl called wrongly')
       call caserr('invalid argument in shiftl')
      endif
c
      a = x
      if (n.eq.0) then
         b = a
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else if (n.lt.32) then
         ib(1) = ior(ishft(ia(1),n),ishft(ia(2),n-32))
         ib(2) = ishft(ia(2),n)
      else if (n.eq.32) then
         ib(1) = ia(2)
         ib(2) = 0
      else
         ib(1) = ishft(ia(2),n-32)
         ib(2) = 0
      end if
c
      shiftl = b
c
      return
      end
c======================================================================
      function shiftr(x,n)
c
c...   real*8 shiftr  (right shift without propagation)
c...   n <0 v >64 shifts aborts
c...   n = 64 must clear word !!
c...   lshft or rshft by 32 does not empty word !!
c
      implicit real*8  (a-h,o-z), integer(i-n)
c
      dimension ia(2),ib(2)
      equivalence (a,ia),(b,ib)
c
      if (n.gt.64.or.n.lt.0) then
       write(6,1000)
1000   format(' shiftr called wrongly')
       call caserr('invalid argument in shiftr')
      endif
c
      a = x
c
      if (n.eq.0) then
         b = a
      else if (n.eq.32) then
         ib(2) = ia(1)
         ib(1) = 0
      else if (n.lt.32) then
         ib(2) = ior(ishft(ia(2),-n),ishft(ia(1),32-n))
         ib(1) = ishft(ia(1),-n)
      else if (n.eq.64) then
         ib(1) = 0
         ib(2) = 0
      else
         ib(2) = ishft(ia(1),32-n)
         ib(1) = 0
      end if
c
      shiftr = b
c
      return
      end
