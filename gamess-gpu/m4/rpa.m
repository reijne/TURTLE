      subroutine daxmyz(n,a,x,y,z)
      implicit REAL (a-h,o-z)
      dimension x(*),y(*),z(*)
      if(n.eq.0) return
      if(n.lt.0) call caserr('daxmyz called with vector length < 0')
      do 10 i=1,n
10    z(i)=a*x(i)-y(i)
      return
      end
      function dnrsq (n,x,incx)
      implicit REAL (a-h,p-z)
      dimension x(*)
      dnrsq = 0.0d0
      if (n .le. 0) return
      incxa = iabs(incx)
      ix = 1
      do 10 i = 1 , n
          dnrsq = dnrsq + x(ix) * x(ix)
          ix = ix + incxa
10    continue
      return
      end
      subroutine getdia (odirct,ipo,nbasis,diag,ind,q)
c
c---------------------------------------------------------------------
c   1. Gets the diagonal of the RPA matrix A (for conventional RPA)
c      or the diagonal consisting of the orbital energy differences
c      eps_a - eps_i (for direct RPA).
c   2. Determines the lowest entries of the diagonal.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit REAL (a-h,p-z), integer (i-n), logical (o)
INCLUDE(common/sizes)
      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,
     +            newbas,ncol,ivalue,ioccup,ispa,inx(maxorb),
     +            ipoint(maxorb),isit(maxorb),itogrp(maxat),
     +            ilabel(maxorb),iok
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
INCLUDE(common/infoa)
INCLUDE(common/iofile)
      common/symmet/isym,isize,neig,iperm(8,8),maxspa,intspa,ineed0,
     + ipadsymmet(10)
c
      dimension diag(*),ind(*),q(*),ipo(nbasis,*)
c
      if (odirct) then
         do 20 ii = 1 , na
            do 10 ia = nc1 , num
               if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
                 q(ipo(ia,ii)) = evalue(ia) - evalue(ii)
               end if
10          continue
20       continue
      else
         jblo = jstart(isym)
         np = npair(isym)
         ipair = 0
         do 40 ibox = 1 , nbox(isym)
            if (ibox .eq. nbox(isym)) np = isize - ipair
            npp = np * isize
            call rdedx_less(q,npp,jblo,num8)
            do 30 n = 1 , np
               diag(ipair+n) = q((n-1)*isize+ipair+n)
30          continue
            ipair = ipair + np
            jblo = jblo + lensec(2*npp)
40       continue
         call dcopy(isize,diag,1,q,1)
      end if
c
      osort = .false.
      if (odebug(19)) osort = .true.
      call quick (q,isize,nstart,ind,osort)
      if (odebug(19)) then
         write (iwr,50)
50       format (/'the lowest entries of the diagonal:'/)
          do 90 i = 1 , nstart
            do 80 ii = 1 , na
               do 70 ia = nc1 , num
               if (iperm(isymmo(ii),isymmo(ia)) .ne. isym) go to 70
               if (ipo(ia,ii) .eq. ind(i)) then
                  write (iwr,60) i,ii,ia,q(i)
60                format (2i6,' -->',i6,8x,f14.9)
                  go to 90
                  end if
70             continue
80          continue
90       continue
      end if
c
      return
      end
      function imaxf(n,x,incx)
      implicit integer (a-z)
      dimension x(*)
      if (n .gt. 0)then
        imaxf = x(1)
        incxa = iabs(incx)
        ix = 1 + incxa
        do 10 i = 2 , n
            imaxf = max(imaxf,x(ix) )
            ix = ix + incxa
10      continue
      else
        imaxf = 0
      endif
      return
      end
      subroutine inbuck (iad,iblo,irun,isym,ipos,jpos,itype,q)
c
c---------------------------------------------------------------------
c   Distributes an RPA integral into the correct bucket
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit REAL (a-h,p-z), integer (i-n), logical (o)
_IFN1(ck)      integer*2 itemp(4)
_IF1(ck)      integer itemp(4)
      logical btest
INCLUDE(common/sizes)
INCLUDE(common/atmblk)
      common/blkin /gin(511)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      parameter (nvcmax=9953)
      common/dipblk/ifill(nvcmax),index(nvcmax)
      dimension q(*),iad(980)
_IFN1(ck)      equivalence (pack,itemp(1))
_IFN1(ck)      data itemp /4*0/
      ibox = nbox0(isym) - 1 + (ipos-1)/npair(isym) + 1
      j = ifill(ibox) + 1
      mlower = (ibox-1)*511
      q(mlower+j) = gin(irun)
      kpos = mlower + num2ep + j/2
      if (btest(j,0)) then
         itemp(1) = 10*ipos + itype
         itemp(2) = jpos
_IF(cray,t3d)
         pack = shiftl(itemp(1),48).or.shiftl(itemp(2),32)
_ENDIF
        q(kpos+1) = pack
      else
        pack = q(kpos)
         itemp(3) = 10*ipos + itype
         itemp(4) = jpos
_IF(cray,t3d)
         pack = pack.or.(shiftl(itemp(3),16).or.itemp(4))
_ENDIF
        q(kpos) = pack
      end if
      if (j .ge. num2ep) then
        q(mlower+511) = dble(iad(ibox))
        call wrt3 (q(mlower+1),511,iblo,5)
        iad(ibox) = iblo*1000 + num2ep
        iblo = iblo + 1
        call vclr(q(mlower+1),1,511)
        ifill(ibox) = 0
      else
        ifill(ibox) = j
      end if
      return
      end
      function isum (n,x,incx)
      implicit integer (a-z)
      dimension x(*)
      isum = 0
      if (n .le. 0) return
      incxa = iabs(incx)
      ix = 1
      do 10 i = 1 , n
          isum = isum + x(ix)
          ix = ix + incxa
10    continue
      return
      end
      subroutine outeig (itable,iw,eigen,neig,niter,ntot,nevlo1,mode)
c
c---------------------------------------------------------------------
c   MODE = 1 :  Write eigenvalues of each iteration to fortran 
c               file ITABLE
c   MODE = 2 :  Read eigenvalues from ITABLE and write them to
c               standard output in the correct order
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit REAL (a-h,p-z)
      logical opg_root
      dimension eigen(*),buf(5)
      m = 1
      n = 5
      kk = 1
      go to (10,30), mode

c   write mode
c
10    nn = min(n,neig)
      if(opg_root()) then
       write (itable) niter,ntot,(eigen(i),i=m,nn)
      endif
      if (neig .le. n) go to 20
      m = m + 5
      n = n + 5
      go to 10
20    return
c
c   read mode
c
30    nn = min(n,neig)
      if(opg_root()) then
      write (iw,40) (nevlo1+i,i=m,nn)
40    format (15x,5i14)
      rewind itable
       do i = 1 , niter
          do j = 1 , kk-1
           read (itable)
          enddo
          read (itable) ibuf,jbuf,(buf(k),k=1,nn-m+1)
         write (iw,60) ibuf,jbuf,(buf(k),k=1,nn-m+1)
60        format (i5,i11,3x,5f14.9)
          do j = kk+1 , (neig-1)/5+1
           read (itable)
          enddo
       enddo
      endif
      if (neig .le. n) return
      if(opg_root()) write (iw,90)
90    format (/)
      m = m + 5
      n = n + 5
      kk = kk + 1
      go to 30
c
      end
      subroutine outosc
c
c---------------------------------------------------------------------
c   Output of total oscillator strengths
c   (c) Carsten Fuchs 1993
c---------------------------------------------------------------------
c
      implicit REAL (a-h,p-z), integer (i-n), logical (o)
INCLUDE(common/sizes)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,
     +              oskip,ototal,oall
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/rpacom)
c
      data c23/0.6666666666666666666666667d0/

      tottda = 2.0d0 * c23 * tottda
      totrpa = c23 * totrpa
      write (iwr,70)
      write (iwr,10)
      if (modus .le. 2) then
         write (iwr,20) 'TDA',fetda
         if (ototal) then
            if (oall) then
               write (iwr,30) tottda
            else
               write (iwr,40) tottda
            end if
            write (iwr,60) tottda - fetda
         end if
      end if
      if (modus .ge. 2) then
         write (iwr,20) 'RPA',ferpa
         if (ototal) then
            if (oall) then
               write (iwr,30) totrpa
            else
               write (iwr,40) totrpa
            end if
            write (iwr,60) totrpa - ferpa
         end if
         write (iwr,50) dble(ne)
      end if
      write (iwr,10)
      write (iwr,70)

      return
10    format (75('*'))
20    format ('*',73x,'*'/'*   ',a3,67x,'*'/'*   ---',
     +            67x,'*'/'*',73x,'*'/
     +        '*   Sum of oscillator strengths calculated:',
     +         15x,f9.4,7x,'*')
30    format ('*   Total oscillator strength expected:',19x,f9.4,7x,'*')
40    format ('*   Total oscillator strength', 
     +        ' in the irreps considered:   ',f9.4,7x,'*')
50    format ('*   Total oscillator strength',
     +        ' with a complete basis set:  ',f9.4,7x,'*'/'*',73x,'*')
60    format ('*   Difference (oscillator strength not exhausted):',
     +         7x,f9.4,7x,'*')
70    format (//122('=')//)
      end
      subroutine outpol (nfreq,freq,pola,ipola,buf)

c---------------------------------------------------------------------
c   Output of the polarizability tensors
c   (c) Carsten Fuchs 1991
c---------------------------------------------------------------------

      implicit REAL (a-h,p-z)
INCLUDE(common/iofile)
      dimension freq(nfreq),pola(6,*),ipola(3,3),buf(*)
      write (iwr,10)
10    format (/'Polarizability tensors'/22('=')/)
      m = 1
      n = 3
20    if (nfreq .lt. m) go to 100
      n = min(n,nfreq)
      write (iwr,30) ('frequency ',freq(i),i=m,n)
30    format (3(a10,f9.4,' a.u.',25x))
      write (iwr,40)
40    format (/)
      do 50 icoor1=1,3
50    write (iwr,60) ((pola(ipola(icoor1,icoor2),i),icoor2=1,3),i=m,n)
60    format (3(3f13.4,10x))
      write (iwr,70) ('averaged polarizability:',i=m,n)
70    format (/3(a24,25x))
      j = 0
      do 80 i=m,n
      j = j+1
80    buf(j) = (pola(1,i)+pola(4,i)+pola(6,i)) / 3.0d0
      write (iwr,90) (buf(i),i=1,j)
90    format (/3(f13.4,36x))
      m = m + 3
      n = n + 3
      go to 20
100   continue
      return
      end
      subroutine outres (iw,inorms,res,nevlo1,neig,niter,mode)
c
c----------------------------------------------------------------------
c MODE = 1: Write norms of residue vectors of each iteration to fortran 
c            file INORMS
c MODE = 2: Read norms of residue vectors from INORMS and write them to
c               standard output in the correct order, representing
c               converged roots as blanks
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
c
      implicit REAL (a-h,p-z)
      logical opg_root
      character*100 line
      character*1 digit(10)
      dimension res(*),buf(20)
      data digit/'0','1','2','3','4','5','6','7','8','9'/
c
c   standard: 5 roots per line, format f14.9
c
      data nroots,iform1,iform2 / 5,14,9 /
c
      m = 1
      n = nroots
      kk = 1
      go to (10,200), mode
c
c   write mode
c
10    nn = min(n,neig)
      if (opg_root()) write (inorms) niter,(res(i),i=m,nn)
      if (neig .le. n) then
         return
      end if
      m = m + nroots
      n = n + nroots
      go to 10
c
c   read mode
c
200   if (opg_root()) then
30    continue
      if (nroots .gt. 20) return
      nn = min(n,neig)
      write (iw,40) (nevlo1+i,i=m,nn)
c...next format should be changed to (15x,<nroots>i<iform1>)
c...if you want different values of nroots,iform1
40    format (15x,5i14)
      rewind inorms
      do 140 i = 1 , niter
         do 50 j = 1 , kk-1
            read (inorms)
50       continue
         read (inorms) idum,(buf(k),k=1,nn-m+1)
         do 110 k = 1 , nn-m+1
            temp = buf(k) + 0.5d0 * (10.0d0**(-iform2))
            if (temp .lt. 0.0d0) go to 90
            if (temp .ge. 10.0d0**(iform1-iform2-1)) go to 70
            t = temp
            ipos = k * iform1 - iform2
            line(ipos:ipos) = '.'
            do 60 ii = iform1-iform2-2 , -iform2 , -1
               power = 10.0d0**ii
               u = t / power
               iu = u
               ipos = k * iform1 - iform2 - 1 - ii + min(1,max(0,-ii))
               line(ipos:ipos) = digit(iu+1)
               if (iu.eq.0 .and. ii.ge.1) line(ipos:ipos) = ' '
               t = t - dble(iu) * power
60          continue
            go to 110
70          do 80 ipos = (k-1)*iform1+1 , k*iform1
               line(ipos:ipos) = '*'
80          continue
            go to 110
90          do 100 ipos = (k-1)*iform1+1 , k*iform1
               line(ipos:ipos) = ' '
100         continue
110      continue
         write (iw,120) i,(line(ipos:ipos),ipos=1,(nn-m+1)*iform1)
c...next format should be changed to (i5,14x,<nroots*iform1>a1)
c...if you want different values of nroots,iform1
120      format (i5,14x,70a1)
         do 130 j = kk+1 , (neig-1)/nroots+1
            read (inorms)
130      continue
140   continue
      if (neig .le. n) return
      write (iw,150)
150   format (/)
      m = m + nroots
      n = n + nroots
      kk = kk + 1
      go to 30
      endif

      end
      subroutine outtdm(ipo, nbasis, method,scfvec,vec,lenvec,
     &     jroot,rnorm,qdmo,qdao,q)

      implicit REAL (a-h,p-y), integer (i-n), logical (o)
      implicit character*8 (z)


INCLUDE(common/sizes)
      common/infoa /nat,ich,mul,nba,nx,ne,ncore,nee,czan(maxat),
     +              c(3,maxat),amass(maxat),iannn(2,maxat),
     +              ipseud(maxat),symz(maxat),lpseud
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     ?              lcore(8),lvirt(8),iofset(8,8),npair(8),
     ?              nbox(8),nbox0(8),jstart(8),iofsym(8)

cjmht      common/block /iky(maxorb+1),ilifq(maxorb+1)
INCLUDE(common/blocko)

      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     ?              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     ?              inozer,nozer(3),idone

      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv

      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)

      dimension q(*), qdmo(*), qdao(*)

      character*3 method
      character*8 sym, root
      character*80 tdm_file
      dimension scfvec(*), vec(*), ipo(nbasis,*)

      do 10 i=1,lenvec
10    vec(i) = vec(i)/rnorm
      write(sym,200) isym
      write(root,200) jroot
200   format(i1)
      l1=index(sym,' ')-1
      tdm_file = 'tdm'//'_'//method//'_'//sym(1:l1)//'_'//root
      ltri =  ((nba+1)*nba)/2
c     lsqr =  nba*nba
      idsq1 = 1
      do 40 i = 1 , maxorb
         k = i*(i-1)/2
        iky(i) = k
40    continue
      do 1 i=1,ltri
         qdmo(i) = 0.0d0
 1    continue
      do 100 ii = 1 , ncore
         do 90 ia = nc1 , nba
            k = iperm(isymmo(ii),isymmo(ia))
            if( k .eq. isym ) then
               qdmo(iky(ia) + ii) = vec(ipo(ii,ia))
            endif
 90      continue
 100  continue
c
      call demoao(qdmo,qdao,scfvec,q(idsq1),
     &     nba,nba,nba)
c
      call square (q(idsq1),qdao,nba,nba)
c
      iden=99
      open(unit=iden,form='unformatted',status='unknown',file=tdm_file)
      call prden(iden,q(idsq1),nba)
      close(unit=iden)
c
      return
      end
      subroutine prden(iunit,den,n)

      integer iunit
      REAL den(n,n)

      write(iunit) den

      end
      subroutine quick (q,n,k,index,osort)
c
c--------------------------------------------------------------------
c   Determines the k lowest entries of the vector q of dimension n
c   using a quicksort type algorithm
c   If osort=true, the lowest entries are ordered
c   (c) Carsten Fuchs 1991-1993
c--------------------------------------------------------------------
c
      implicit REAL (a-h,p-z), logical (o)
      dimension q(*),index(*)
      if (k .gt. n) go to 70
      do 10 i=1,n
         index(i) = i
10    continue
      ip = 1
      ir = n
20    t = q(ip)
      i = ip - 1
      j = ir + 1
30    j = j - 1
      if (q(j) .gt. t) go to 30
40    i = i + 1
      if (q(i) .lt. t) go to 40
      if (i .lt. j) then
        qhelp = q(i)
        q(i) = q(j)
        q(j) = qhelp
        ihelp = index(i)
        index(i) = index(j)
        index(j) = ihelp
        go to 30
      end if
      if (k .gt. j) then
        ip = j + 1
      else
        ir = j
      end if
      if (ip .ne. ir) go to 20
      if (osort) then
         do 60 i=1,k
            m = i
            qmin = q(i)
            do 50 j=i+1,k
               if (q(j) .lt. qmin) then
                 m = j
                  qmin = q(j)
               end if
50          continue
            t = q(i)
            q(i) = q(m)
            q(m) = t
            l = index(i)
            index(i) = index(m)
            index(m) = l
60       continue
      end if
      return
70    call caserr ('invalid parameters for routine quick')
      stop
      end
      subroutine redeig (maxr,ntot,sred1,ered1,bred,sred2,ered2)

c----------------------------------------------------------------------
c  Constructs the reduced matrices
c
c             |  Sigma~  Delta~  |                  |  A~  B~  |
c      S~  =  |                  |    and    E~  =  |          |
c             | -Delta~ -Sigma~  |                  |  B~  A~  |
c
c  (in arrays SRED2,ERED2) for the iterative response eigenvalue
c  algorithm from the matrices SRED1,ERED1 and the vector BRED,
c  which contain the entries of S~ and E~ in a packed form.
c  The upper triangle of SRED1 contains the upper triangle of Sigma~,
c  which will then be completed to a symmetric matrix, and
c  the strict lower triangle of SRED1 contains the strict lower triangle
c  of Delta~, which will be completed to an antisymmetric matrix.
c  The upper triangle of ERED1 contains the upper triangle of A~,
c  which will then be completed to a symmetric matrix, and
c  the strict lower triangle of ERED1 contains the strict lower triangle
c  of B~, which will also be completed to an symmetric matrix, together
c  with the diagonal of B~ supplied in the vector BRED.
c  (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
  
      implicit REAL (a-h,p-z)
      dimension sred1(maxr,maxr),ered1(maxr,maxr),bred(maxr)
      dimension sred2(2*maxr,2*maxr),ered2(2*maxr,2*maxr)

      do 200 i = 1 , ntot
         ii = i + ntot
         ered2(i,i)   = ered1(i,i)
         ered2(ii,ii) = ered1(i,i)
         sred2(i,i)   = sred1(i,i)
         sred2(ii,ii) = - sred1(i,i)
         ered2(i,ii)  = bred(i)
         ered2(ii,i)  = bred(i)
         sred2(i,ii)  = 0.0d0
         sred2(ii,i)  = 0.0d0
         do 100 j = i+1 , ntot
            jj = j + ntot
            t = ered1(i,j)
            ered2(i,j)   = t
            ered2(j,i)   = t
            ered2(ii,jj) = t
            ered2(jj,ii) = t
            t = sred1(i,j)
            sred2(i,j)   = t
            sred2(j,i)   = t
            sred2(ii,jj) = - t
            sred2(jj,ii) = - t
            t = ered1(j,i)
            ered2(i,jj)  = t
            ered2(jj,i)  = t
            ered2(j,ii)  = t
            ered2(ii,j)  = t
            t = sred1(j,i)
            sred2(j,ii)  = t
            sred2(ii,j)  = t
            sred2(i,jj)  = - t
            sred2(jj,i)  = - t
 100     continue
 200  continue
      return
      end
      subroutine redeqs (maxr,ntot,freq,sred,ered,bred,ared,fred2,
     +                   ipvt,ierr)
c
c----------------------------------------------------------------------
c Constructs and solves the reduced linear response equation system
c
c         (omega * S~ - E~) c~ = f~
c
c The matrix ARED containing  omega * S~ - E~  is constructed from
c matrices SRED,ERED and vector BRED as described in subroutine REDEIG.
c The right hand side is provided in FRED2, and FREQ is the 
c frequency omega.
c   (c) Carsten Fuchs 1991-1993
c-----------------------------------------------------------------------
c
      implicit REAL (a-h,p-z), logical (o)
      common/debug /odebug(100)
      dimension sred(maxr,maxr),ered(maxr,maxr),bred(maxr),
     +          ared(2*maxr,2*maxr)
      dimension fred2(2*maxr),ipvt(2*maxr)
c
      do 200 i = 1 , ntot
         ii = i + ntot
         ared(i,i)   =  freq * sred(i,i) - ered(i,i)
         ared(ii,ii) = -freq * sred(i,i) - ered(i,i)
         ared(i,ii)  = -bred(i)
         ared(ii,i)  = -bred(i)
         do 100 j = i+1 , ntot
            jj = j + ntot
            ared(i,j)   =  freq * sred(i,j) - ered(i,j)
            ared(j,i)   =  ared(i,j)
            ared(ii,jj) = -freq * sred(i,j) - ered(i,j)
            ared(jj,ii) =  ared(ii,jj)
            ared(i,jj)  = -freq * sred(j,i) - ered(j,i)
            ared(jj,i)  =  ared(i,jj)
            ared(j,ii)  =  freq * sred(j,i) - ered(j,i)
            ared(ii,j)  =  ared(j,ii)
100      continue
200   continue
c
      ntot2 = 2 * ntot
      if (odebug(90)) call outsqr (ared,ntot2,ntot2,ntot2,'matrix ared')
      call dsifa (ared,2*maxr,ntot2,ipvt,ierr)
      if (ierr .ne. 0) return
      call dsisl (ared,2*maxr,ntot2,ipvt,fred2)
c
      return
      end
      subroutine rpalgo (diag,b1,b2,sred1,ered1,bred,sred2,ered2,
     +  etare,etaim,sigma,xred,h1,h2,v1,v2,resarr,scr)
c
c----------------------------------------------------------------------
c   Performs the iterative solution of the generalized RPA eigenvalue
c   equation
c                 
c       (         |  S   0  |   |  A   B  | ) | y |   | 0 |
c       ( omega * |         | - |         | ) |   | = |   |
c       (         |  0  -S  |   |  B   A  | ) | z |   | 0 |
c
c   using the iterative algorithm by Olsen & Jorgensen 
c   (J. Comp. Phys. 74, 265 (1988))
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
                   
      implicit REAL (a-h,p-z), integer (i-n), logical (o)
INCLUDE(common/sizes)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      parameter (nvcmax=9953)
      common/dipblk/iconv(nvcmax),index(nvcmax)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr
      common/junk2 /ix(1360)
INCLUDE(common/rpacom)
INCLUDE(common/rpfile)
      common/symmet/isym,isize,neig,ipadsymmet(77)
      dimension diag(isize)
      dimension b1(isize,maxr),b2(isize,maxr)
      dimension sred1(maxr,maxr),ered1(maxr,maxr),bred(maxr)
      dimension sred2(maxr2,maxr2),ered2(maxr2,maxr2)
      dimension etare(maxr2),etaim(maxr2),sigma(maxr2),xred(maxr2,maxr2)
      dimension h1(isize),h2(isize)
      dimension v1(isize,nstart),v2(isize,nstart)
      dimension resarr(*),scr(*)
      data ieigen,iresid / 20,21 /

      ntot1 = 1
      ntot = nstart
      if (odebug(69)) call outvec (diag,isize,'diag in rpalgo')
10    niter = niter + 1
      if (niter .gt. maxit) go to 1000
c
c   calculate the vectors e * b(i) and fill up matrices sred1,ered1
c
      do 20 i=ntot1,ntot
      time = cpulft(1)
      call rpmult (.true.,jstart(isym),npair(isym),isize,num8,
     +   b1(1,i),b2(1,i),h1,h2,scr)
      tim(2) = tim(2) + cpulft(1) - time
      bred(i) = ddot(isize,b1(1,i),1,h2,1) 
     +        + ddot(isize,b2(1,i),1,h1,1)
      if (odebug(30)) then
        print*,'result vector no',i
        call outvec (h1,isize,'h1')
        call outvec (h2,isize,'h2')
      end if
      do 20 j=1,i
      ered1(j,i) = ddot(isize,b1(1,j),1,h1,1) 
     +    + ddot(isize,b2(1,j),1,h2,1)
      sred1(j,i) = ddot(isize,b1(1,i),1,b1(1,j),1) 
     +    - ddot(isize,b2(1,i),1,b2(1,j),1)
      if (j .eq. i) go to 20
      ered1(i,j) = ddot(isize,b1(1,j),1,h2,1) 
     +    + ddot(isize,b2(1,j),1,h1,1)
      sred1(i,j) = ddot(isize,b1(1,i),1,b2(1,j),1) 
     +    - ddot(isize,b2(1,i),1,b1(1,j),1)
20    continue
c
c   construct ered2,sred2 for library programs
c
      call redeig (maxr,ntot,sred1,ered1,bred,sred2,ered2)
      ntot2 = 2 * ntot
      if (odebug(17)) then
         print*,'sred2 & ered2, iteration step ',niter
         call outsqr (sred2,maxr2,ntot2,ntot2,'sred2')
         call outsqr (ered2,maxr2,ntot2,ntot2,'ered2')
      end if
c
c   solve reduced eigenvalue problem
c
      if (oreduc) go to 405
c
c   qz algorithm
c
      time = cpulft(1)
      call qzhes (maxr2,ntot2,ered2,sred2,.true.,xred)
      call qzit (maxr2,ntot2,ered2,sred2,0.0d0,.true.,xred,ierr)
      call qzval (maxr2,ntot2,ered2,sred2,etare,etaim,sigma,.true.,xred)
      if (ierr .ne. 0) go to 400
      call qzvec (maxr2,ntot2,ered2,sred2,etare,etaim,sigma,xred)
      tim(3) = tim(3) + cpulft(1) - time
      if (odebug(18)) then
         write (iwr,30) niter
30       format ('eigenvalues from QZ algorithm, iteration step ',i3/
     +           '-------------------------------------------------'/)
         call outvec (etare,ntot2,'real parts of eta')
         call outvec (etaim,ntot2,'imaginary parts of eta')
         call outvec (sigma,ntot2,'sigma')
      end if
c
c   collect real positive eigenvalues
c
      nrealp = 0
      do 40 i=1,ntot2
      t = etare(i)
      u = etaim(i)
      v = sigma(i)
      if (1.0d0+u .ne. 1.0d0 .or. 1.0d0+v .eq. 1.0d0) go to 40
      t = t / v
      if (t .lt. 0.0d0) go to 40
      nrealp = nrealp + 1
      etare(nrealp) = t
      call dcopy(ntot2,xred(1,i),1,xred(1,nrealp),1)
40    continue
      if (nrealp .eq. 0) go to 500
      if (nrealp .lt. nstart) then
         write (iwr,50) niter,nrealp
50       format ('warning: in iteration step ',i2,' there are only ',
     +   i2,' positive real eigenvalues.'/'trying to continue ...')
         call vclr(etare(nrealp+1),1,nstart-nrealp)
      end if
      if (odebug(18)) write (iwr,51) niter,nrealp
51    format ('iteration step ',i3,': ',i3,' positive real eigenvalues')
      nev2 = min(nrealp,nstart)
c
c   sort eigenvalues
c
      do 60 i=1,nev2
      k = i
      if (nrealp .gt. 1) then
         emin = etare(i)
         do 55 j=i+1,nrealp
         if (etare(j) .lt. emin) then
            k = j
            emin = etare(j)
         end if
55    continue
      end if
      t = etare(i)
      etare(i) = etare(k)
      etare(k) = t
      call dswap(ntot2,xred(1,i),1,xred(1,k),1)
60    continue
c      nev3 = nev2 - nevlo1
      go to 65
c
c   reduction to ordinary eigenvalue problem
c
405   time = cpulft(1)
      call reduc (maxr2,ntot2,sred2,ered2,etaim,ierr)
      if (ierr .ne. 0) go to 600
      call tred2(maxr2,ntot2,sred2,etare,sigma,xred)
      call tql2(maxr2,ntot2,etare,sigma,xred,ierr)
      if (ierr .ne. 0) go to 700
      call rebak (maxr2,ntot2,ered2,etaim,ntot2,xred)
      tim(3) = tim(3) + cpulft(1) - time
      do 406 i=1,nstart
      etare(i) = 1.0d0 / etare(ntot2+1-i)
      call dcopy(ntot2,xred(1,ntot2+1-i),1,xred(1,i),1)
406   continue
      nev2 = nstart
c      nev3 = nstart - nevlo1
c
65    call outeig (itable,iwr,etare(nevlo1+1),neig,niter,ntot,nevlo1,1)
      if (opg_root()) then
      open (ieigen,file='roots',status='unknown')
      endif
      call outeig (itable,ieigen,dummy,neig,niter,ntot,nevlo1,2)
      if (opg_root()) close (ieigen)
c
c   new trial vectors to be expected
c
      ntot0 = ntot
      ntot1 = ntot + 1
      do 150 j=nevlo(isym),nev2
      if (iconv(j-nevlo1) .ne. 0) then
         resarr(j-nevlo1) = -1.0d0
         go to 150
      end if
c
c   compute v(j) for non-converged eigenvalue no j
c
      do 110 i=1,isize
      v1(i,j) = 0.0d0
      v2(i,j) = 0.0d0
      do 110 k=1,ntot0
      t = b1(i,k)
      u = b2(i,k)
      v = xred(k,j)
      w = xred(k+ntot0,j)
      v1(i,j) = v1(i,j) + t * v + u * w
      v2(i,j) = v2(i,j) + u * v + t * w
110   continue
c
c   compute residue vector
c
      time = cpulft(1)
      call rpmult (.true.,jstart(isym),npair(isym),isize,num8,
     +   v1(1,j),v2(1,j),h1,h2,scr)
      tim(2) = tim(2) + cpulft(1) - time
      freq = etare(j)
_IFN1(civu)      call vsmsb(v1(1,j),1,freq,h1,1,h1,1,isize)
_IFN1(civu)      call vsmsb(v2(1,j),1,-freq,h2,1,h2,1,isize)
_IF1(civu)      call gtrian(isize,freq,h1,h1,v1(1,j))
_IF1(civu)      call gtrian(isize,-freq,h2,h2,v2(1,j))
      if (odebug(47)) then
         write (iwr,61) j
61       format ('residue vector no.',i3,':')
         call outvec (h1,isize,'first  part')
         call outvec (h2,isize,'second part')
      end if
c
c   compute euclidean norm of residue vector
c
      t = dsqrt(dnrsq(isize,h1,1)+dnrsq(isize,h2,1))
      resarr(j-nevlo1) = t
c
c   threshold ?
c
      if (t .le. eps) iconv(j-nevlo1) = 1
      if (     (iconv(j-nevlo1) .eq. 1)
     +  .or. (ntot .ge. maxr)
     +  .or. (ntot .ge. isize)       ) go to 150
      ntot = ntot + 1
c
c   update residue vector 
c
      do 115 i = 1 , isize
         b1(i,ntot) =  freq - diag(i)
         b2(i,ntot) = -freq - diag(i)
115   continue
      call vdivz (h1,1,b1(1,ntot),1,1.0d0,b1(1,ntot),1,isize)
      call vdivz (h2,1,b2(1,ntot),1,1.0d0,b2(1,ntot),1,isize)
      if (odebug(68)) then
         print*,'updated residue vector no. ',j
         call outvec (b1(1,ntot),isize,'first part')
         call outvec (b2(1,ntot),isize,'second part')
      endif
c
c   gram-schmidt orthogonalisation
c
      do 120 k=1,ntot-1
      t = - ddot(isize,b1(1,ntot),1,b1(1,k),1) 
     +  - ddot(isize,b2(1,ntot),1,b2(1,k),1)
      call daxpy(isize,t,b1(1,k),1,b1(1,ntot),1)
      call daxpy(isize,t,b2(1,k),1,b2(1,ntot),1)
      t = - ddot(isize,b1(1,ntot),1,b2(1,k),1) 
     +  - ddot(isize,b2(1,ntot),1,b1(1,k),1)
      call daxpy(isize,t,b2(1,k),1,b1(1,ntot),1)
120   call daxpy(isize,t,b1(1,k),1,b2(1,ntot),1)
      t = dsqrt(dnrsq(isize,b1(1,ntot),1)+dnrsq(isize,b2(1,ntot),1))
      if (1.0d0+t .eq. 1.0d0) then
         write (iwr,130) j,niter
130      format ('discard new trial vector no. ',i2,' in iteration',
     +           ' step ',i2/'because of linear dependence from',
     +           ' the previous trial vectors.')
         ntot = ntot - 1
      else
c
c   perform symmetric orthonormalisation
c
         call symort (isize,b1(1,ntot),b2(1,ntot),sytol,oldep)
         if (oldep) then
            write (iwr,140) j,niter
140         format ('discard new trial vector no. ',i2,' in',
     +              ' iteration step ',i2/'because b and b#',
     +              ' are (nearly) linearly dependent.')
            ntot = ntot - 1
         end if
      end if
150   continue
      call outres (iwr,inorms,resarr(nevlo(isym)-nevlo1),nevlo1,
     +             neig,niter,1)
      if (opg_root()) then
       open (iresid,file='residues',status='unknown')
      endif
      call outres (iresid,inorms,dummy,nevlo1,neig,niter,2)
      if (opg_root()) close (iresid)
      nconv = isum(neig,iconv,1)
      if (nconv .eq. neig) go to 1000
c
c   if at least one new trial vector has been added,
c   start next iteration step
c
      if (ntot .ge. ntot1) go to 10
      if (ntot .ge. maxr) then
         do 160 j=1,nevlo1
         do 160 i=1,isize
         v1(i,j) = 0.0d0
         v2(i,j) = 0.0d0
         do 160 k=1,ntot0
         t = b1(i,k)
         u = b2(i,k)
         v = xred(k,j)
         w = xred(k+ntot0,j)
         v1(i,j) = v1(i,j) + t * v + u * w
         v2(i,j) = v2(i,j) + u * v + t * w
160      continue
         return
      end if
      write (iwr,170) niter
170   format ('no new trial vector has been added in iteration step ',
     +   i2, ':'/
     +  'unable to finish iterative calculation of rpa eigenvalues.')
      ierr = 1
      go to 1000
c
c   error messages
c
400   write (iwr,410) niter,ierr
410   format (/'convergence problems with QZIT in iteration step',
     +        i3,': ierr = ',i3/)
      go to 1000
c
500   write (iwr,510) niter
510   format (/'in iteration step ',i2,
     +         'there are no real eigenvalues.'/)
      ierr = 1
      go to 1000
c
600   write (iwr,610) niter,ierr
610   format (/'problem in iteration step ',i2,
     +         ': matrix ERED2 is not positive',
     +       ' definite.'/'ierr = ',i4)
      go to 1000
c
700   write (iwr,710) niter,ierr
710   format (/'convergence problems with TQL2 in iteration step',
     +        i2,': ierr = ',i2/)
c
1000  niter = min(niter,maxit)
c
      return
      end
      subroutine rpanal (ipo,nbasis,energy,v1,v2,ycomp,zcomp,
     +                   dip,ibleig,iconv,otda,oincor,mode,qanal)
c
c----------------------------------------------------------------------
c   Performs computation of oscillator strength, analysis of the 
c   eigenvectors and output of RPA/TDA results
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit REAL (a-h,p-z), integer (i-n), logical (o)
      character*1 sgn(3)
      character*3 symbol(14),zymbol(14),method
      character*6 symtex(14),zymtex(14)
      character*80 dmpfil,rstfil
c
INCLUDE(common/sizes)
INCLUDE(common/timeperiods)
c
      common/blkin /gin(511)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      common/atmol3/minv(2),mouta,moutb
INCLUDE(common/machin)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr
      common/junk2 /ix(340),jx(340),kx(340),lx(340)
INCLUDE(common/rpacom)
INCLUDE(common/rpfile)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),
     +              maxvec,mvec,mymaxv
      common/blkin1   /evalue(maxorb),hocc(maxorb),etot,nbas2,
     +               newbas,ncol,ivalue,ioccup,ispa,inx(maxorb),
     +               ipoint(maxorb),isit(maxorb),itogrp(maxat),
     +               ilabel(maxorb),iok
      character*8 zlabel,tlabel
      common/junkc/zlabel(maxorb),tlabel(maxorb)
c
      common/table /dmpfil,rstfil,rlines(2)
c
cjmht      common/block /iky(maxorb+1),ilifq(maxorb+1)
INCLUDE(common/blocko)
      common/fsave/ofiles(2)
      dimension energy(*),v1(*),v2(*),ycomp(*),zcomp(*),
     +          dip(*),ibleig(*),iconv(*),ipo(nbasis,*),qanal(*)
      data symbol/'a''','a"','a1','b1','b2','a2',
     +            'ag','b1u','b2u','b3g','b3u','b2g','b1g','au'/
      data zymbol/'A''','A"','A1','B1','B2','A2',
     +            'Ag','B1u','B2u','B3g','B3u','B2g','B1g','Au'/
      data symtex/'a''','a''''','a_1','b_1','b_2','a_2','a_g',
     +            'b_{1u}','b_{2u}',
     + 'b_{3g}','b_{3u}','b_{2g}','b_{1g}','a_u'/
      data zymtex/'A''','A''''','A_1','B_1','B_2','A_2','A_g',
     +            'B_{1u}','B_{2u}',
     +  'B_{3g}','B_{3u}','B_{2g}','B_{1g}','A_u'/
      data ev,c43/27.21165d0,1.33333333333333333333d0/
      data pageth/45.0d0/
c
      call start_time_period(TP_RPANAL)
c
      if (mode.eq.1) then
c
c       initialize tex file(s)
c
        if (modus .le. 2) then
          if (opg_root()) then
            write (itdtex,1010)
            write (itdtex,1100) 'TDA'
            write (itdtex,1110) epstab
          endif
          rlines(1) = 6.3d0
        endif
        if (modus .ge. 2) then
          if (opg_root()) then
            write (irptex,1010)
            write (irptex,1100) 'RPA'
            write (irptex,1110) epstab
          endif
          rlines(2) = 6.3d0
        endif
c
c       initialize spectrum file(s)
c
        if (opg_root()) then
          if (modus .le. 2) write (itdspc,1130) zymbol(nirs+1),
     +       0.0d0,0.0d0
          if (modus .ge. 2) write (irpspc,1130) zymbol(nirs+1),
     +       0.0d0,0.0d0
          if (nev(1) .eq. 0) then
            if (modus .le. 2) write (itdspc,*)
            if (modus .ge. 2) write (irpspc,*)
          end if
        end if
c
c       retrieve information from symmetry assignment
c
        call resyas
c
      else if (mode.eq.2) then
c
c       write table of results
c
        if (otda) then
          itex = itdtex
          ispec = itdspc
          method = 'TDA'
          il = 1
        else
          itex = irptex
          ispec = irpspc
          method = 'RPA'
          il = 2
        end if
        if (opg_root()) write (itex,2010)
        rlines(il) = rlines(il) + 0.3d0
c
        write (iwr,2020) method,isym,zymbol(nirs+isym),epstab
        k1 = 1
        k2 = k1 + isize
c
c       get SCF MO's for transition density matrix transformation
c
        if(ofiles(1)) then
          imo=1
          iscf = imo+nbasq
          idmo = iscf+nbasq
          idao = idmo+nbatri
          iwork = idao+nbatri
          call secget (mouta,3,iblko)
          iblkq = iblko + lensec(mach(8)) + lensec(mach(9)) + 1
          call rdedx (qanal(imo),nbasq,iblkq,idaf)
          if (odebug(5))
     +      call outsqr (qanal(imo),nbasis,nbasis,nbasis,
     +                   'orig mo coeffs')

          do 15 i = 1, nbasis+1
            ilifq(i) = (i-1)*nbasis
15        continue

          call tdown (qanal(iscf),ilifq,qanal(imo),ilifq,nbasis)
        endif
c
c       loop over the excited states that have been computed
c
        do 630 j = 1 , neig
          if (iconv(j) .ne. 1) then
            write (iwr,2035) j
            go to 625
          end if
          if (.not. oincor) then
            if (otda) then
               call rdedx (v1,isize ,ibleig(j),num8)
            else
               call rdedx (v1,isize2,ibleig(j),num8)
            end if
          end if
c
c         normalization
c
          if (otda) then
            fac = 1.0d0
          else
            if (oincor) then
               ycomp(j) = dnrsq(isize,v1(k1),1) 
               zcomp(j) = dnrsq(isize,v2(k1),1)
            else
               ycomp(j) = dnrsq(isize,v1(k1),1) 
               zcomp(j) = dnrsq(isize,v1(k2),1)
            end if
            fac = ycomp(j) - zcomp(j)
            if (fac+1.0d0 .eq. 1.0d0) then
              write (iwr,210) j
              go to 625
            end if
            if (oincor) then
               call vsub(v1(k1),1,v2(k1),1,v1(k1),1,isize)
            else
               call vsub(v1(k1),1,v1(k2),1,v1(k1),1,isize)
            end if
          end if
c
c         output transition density matrix
c
          if(ofiles(1)) call outtdm(ipo,nbasis,method,qanal(iscf),
     +      v1(k1),isize,j,fac,qanal(idmo),qanal(idao),qanal(iwork))
c
c         calculate oscillator strength
c
          fe = 0.0d0
          idip = 1
          do 260 icoord = 1 , 3
            t = ddot(isize,dip(idip),1,v1(k1),1)
            fe = fe + t * t
            idip = idip + nvc
260       continue
          fe = c43 * energy(j) * fe / fac
          if (otda) then
            fetda = fetda + fe
          else
            ferpa = ferpa + fe
          end if
c
c         write rpa eigenvector x1-x2 on tm file
c
          k11 = k1 - 1
          if (.not.otda .and. num.le.511) then
            do 10030 ii = 1 , na
               call vclr(gin,1,num)
               do 10010 ia = nc1 , num
                  if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
                    gin(ia) = v1(k11+ipo(ii,ia)) / fac
                  end if
10010          continue
               if (opg_root()) write (itmfil) (gin(jj),jj=1,num)
10030       continue
            call vclr(gin,1,num)
            do 10040 ia = nc1 , num
               if (opg_root()) write (itmfil) (gin(jj),jj=1,num)
10040       continue
          end if
c
          t = dnrm2(isize,v1(k1),1)
          if (1.0d0+t .eq. 1.0d0) go to 625
          call dscal(isize,1.0d0/t,v1(k1),1)
c
c         store components of eigenvector j with modulus .ge. epstab
c
          n = 0
          do 510 ii = 1 , na
            do 500 ia = nc1 , num
               if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
                 t = v1(k11+ipo(ii,ia))
                 if (dabs(t) .ge. epstab) then
                   n = n + 1
                   if (n .gt. 511) then
                        write (iwr,499) j
                        n = 511
                        go to 515
                   end if
                   gin(n) = t
                   kx(n) = (ii-1)*num + ia
                 end if
               end if
500         continue
510       continue
515       continue
c
c         sort the components
c
          do 530 i = 1 , n
            k = i
            if (n .gt. 1) then
              cmax = dabs(gin(i))
              do 520 i1 = i+1 , n
                if (dabs(gin(i1)) .gt. cmax) then
                  k = i1
                  cmax = dabs(gin(i1))
                end if
520           continue
            end if
            t = gin(i)
            gin(i) = gin(k)
            gin(k) = t
            ihelp = kx(i)
            kx(i) = kx(k)
            kx(k) = ihelp
530       continue
          if (gin(1) .lt. 0.0d0) 
     +      call dscal(n,-1.0d0,gin,1)
c
c         write line of table in output
c
          n2 = min(n,3)
          j1 = j + nevlo1
          if (isym .eq. 1) j1 = j1 + 1
          sgn(1) = ' '
          do 540 i = 2 , n2
            if (gin(i) .gt. 0.0d0) sgn(i) = '+'
            if (gin(i) .lt. 0.0d0) sgn(i) = '-'
540       continue
          do 550 i = 1 , n2
            ix(i) = (kx(i)-1) / num + 1
            jx(i) = kx(i) - (ix(i)-1)*num
550       continue
          if (iok .eq. 2) then
            write (iwr,560) j1,zymbol(nirs+isym),energy(j)*ev,fe,
     +               (sgn(i),dabs(gin(i)),
     +                isymao(ix(i)),symbol(nirs+isymmo(ix(i))),
     +                isymao(jx(i)),symbol(nirs+isymmo(jx(i))),i=1,n2)
          else
            write (iwr,570) j1,zymbol(nirs+isym),energy(j)*ev,fe,
     +                    (sgn(i),dabs(gin(i)),
     +                     ilabel(ix(i)),zlabel(ix(i))(1:4),
     +                     ilabel(jx(i)),zlabel(jx(i))(1:4),i=1,n2)
          end if
          if (opg_root()) 
     +       write (ispec,1130) zymbol(nirs+isym),energy(j),fe
          n1 = 4
          n2 = 6
580       if (n .lt. n1) go to 2100
          n2 = min(n2,n)
          do 590 i = n1 , n2
            if (gin(i) .gt. 0.0d0) sgn(i-n1+1) = '+'
            if (gin(i) .lt. 0.0d0) sgn(i-n1+1) = '-'
590       continue
          do 600 i = n1 , n2
            ix(i) = (kx(i)-1) / num + 1
            jx(i) = kx(i) - (ix(i)-1)*num
600       continue
          if (iok .eq. 2) then
            write (iwr,610) (sgn(i-n1+1),dabs(gin(i)),
     +              isymao(ix(i)),symbol(nirs+isymmo(ix(i))),
     +              isymao(jx(i)),symbol(nirs+isymmo(jx(i))),i=n1,n2)
          else
            write (iwr,620) (sgn(i-n1+1),dabs(gin(i)),
     +                      ilabel(ix(i)),zlabel(ix(i)),
     +                      ilabel(jx(i)),zlabel(jx(i)),i=n1,n2)
          end if
          n1 = n1 + 3
          n2 = n2 + 3
          go to 580
c
c         write line of tex table
c
2100      continue
          if (opg_root()) then
            if (n .gt. 0) then
              if (iok .eq. 2) then
                write (itex,2110) j1,zymtex(nirs+isym),energy(j)*ev,fe,
     +                gin(1),
     +                isymao(ix(1)),symtex(nirs+isymmo(ix(1))),
     +                isymao(jx(1)),symtex(nirs+isymmo(jx(1)))
              else
                write (itex,2115) j1,zymtex(nirs+isym),energy(j)*ev,fe,
     +           gin(1),
     +           ilabel(ix(1)),tlabel(ix(1)),
     +           ilabel(jx(1)),tlabel(jx(1))
              end if
            else
              write (itex,2111) j1,zymtex(nirs+isym),energy(j)*ev,fe,
     +        gin(1)
            end if
          endif
          rlines(il) = rlines(il) + 1.0d0
          do 2130 i = 2 , n
            if (gin(i) .gt. 0.0d0) sgn(1) = '+'
            if (gin(i) .lt. 0.0d0) sgn(1) = '-'
            if (opg_root()) then
              if (iok .eq. 2) then
                 write (itex,2120) sgn(1),dabs(gin(i)),
     +                      isymao(ix(i)),symtex(nirs+isymmo(ix(i))),
     +                      isymao(jx(i)),symtex(nirs+isymmo(jx(i)))
              else
                 write (itex,2125) sgn(1),dabs(gin(i)),
     +           ilabel(ix(i)),tlabel(ix(i)),
     +           ilabel(jx(i)),tlabel(jx(i))
              end if
            endif
            rlines(il) = rlines(il) + 1.0d0
2130      continue
          if (opg_root()) write (itex,2010)
          rlines(il) = rlines(il) + 0.3d0
          if (rlines(il) .ge. pageth) then
            if (opg_root()) write (itex,2140)
            rlines(il) = 3.3d0
          end if
c
625       if (oincor) k1 = k1 + isize
c
630     continue

c
        if (opg_root()) write (ispec,*)
c
c       y and z components
c
        if (.not. otda) then
          write (iwr,380)
          do 420 j = 1 , neig
            if (iconv(j) .ne. 1) go to 420
            yznrsq = ycomp(j) + zcomp(j)
            if (1.0d0+yznrsq .eq. 1.0d0) go to 420
            yscale = ycomp(j) / yznrsq
            zscale = zcomp(j) / yznrsq
            write (iwr,390) j,yscale,zscale
420       continue
        end if
c
c       detailed analysis of excitations
c
        if (.not. oanal) goto 9999
        mx = 0
        my = mx + nvc
        mz = my + nvc
        write (iwr,320)
        k1 = 1
        do 360 j = 1 , neig
          if (iconv(j) .ne. 1) go to 355
          if (.not. oincor) then
            if (otda) then
               call rdedx (v1,isize ,ibleig(j),num8)
            else
               call rdedx (v1,isize2,ibleig(j),num8)
               call vsub(v1(k1),1,v1(k2),1,v1(k1),1,isize)
            end if
          end if
          t = dnrm2(isize,v1(k1),1)
          if (1.0d0+t .eq. 1.0d0) then
            write (iwr,321) j
            go to 355
          end if
          call dscal(isize,1.0d0/t,v1(k1),1)
          total = 0.0d0
          do 345 ii = 1 , na
            do 340 ia = nc1 , num
               if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
                 k = ipo(ii,ia)
                 t = v1(k1-1+k)
                 if (t*t .ge. epsana) then
                   write (iwr,330) j,ii,ia,isymao(ii),
     +              symbol(nirs+isymmo(ii)),
     +              isymao(ia),
     +              symbol(nirs+isymmo(ia)),
     +              t,t*t,dip(mx+k),dip(my+k),dip(mz+k)
                   total = total + t * t
                 end if
               end if
340         continue
345       continue
          write (iwr,350) total
355       if (oincor) k1 = k1 + isize
360     continue
c
      else ! (mode.eq.3)
c
c       closing statements
c
        if (modus .le. 2) then
          if (opg_root()) write (itdtex,2010)
          rlines(1) = rlines(1) + 0.3d0
          if (opg_root()) then
            write (itdtex,3010)
            close (itdtex)
c
            close (itdspc)
c
          endif
        end if
        if (modus .ge. 2) then
          if (opg_root()) write (irptex,2010)
          rlines(2) = rlines(2) + 0.3d0
          if (opg_root()) then
            write (irptex,3010)
            close (irptex)
c
            close (irpspc)
c
          endif
        end if
      endif
 9999 continue
      call end_time_period(TP_RPANAL)
      return
c
350   format (/44x,f10.4,2x,'total'/)
330   format (i3,i5,' -->',i3,4x,i2,a3,' -->',i3,a3,
     +        2f10.4,4x,3f14.4)
321   format ('cannot analyse eigenvector ',i3)
320   format (/'detailed analysis of excitations:'//
     +         '(c refers to component of normalized vector y-z)'//
     +         'no',5x,'i',6x,'a',6x,'i',9x,'a',9x,'c',7x,'c**2',14x,
     +         'my(x)',9x,'my(y)',9x,'my(z)'/)
390   format (i3,5x,f8.4,9x,f8.4)
380   format (/'analysis of y and z components:'//
     +         ' no   y contribution   z contribution'/)
620   format (32x,3(a1,f5.2,' (',i2,a4,' -->',i3,a4,') '))
610   format (32x,3(a1,f5.2,' (',i2,a3,' -->',i3,a3,') '))
570   format (i3,a3,2f10.4,6x,3(a1,f5.2,' (',i2,
     +           a4,' -->',i3,a4,') '))
560   format (i3,a3,2f10.4,6x,3(a1,f5.2,' (',i2,
     +        a3,' -->',i3,a3,') '))
499   format (/'WARNING: too many large components in',
     +        ' eigenvector ',i3,' -- skipping the rest'/)
210   format ('eigenvector no. ',i3,' cannot be scaled.')
2035  format (/'eigenvector ',i2,' has not converged and',
     +         ' cannot be analysed.'/)
2020  format (//34('=')/1x,a3,' results for symmetry',i2,
     +        ' (',a3,')'/34('=')//
     +          ' state    energy  oscillator      leading excitations'/
     +  '           (eV)    strength       (|c| > ',f5.2,')'/107('-'))
1130  format (4x,'1',a3,2f11.6)
c
_IF(doublebackslash)
1010  format ('\\documentstyle[11pt]{article}'/
     +        '\\pagestyle{empty}'/
     +        '\\topmargin=-1.8cm'/
     +        '\\oddsidemargin=-0.5cm'/
     +        '\\evensidemargin=-0.5cm'/
     +        '\\textwidth=17cm'/
     +        '\\textheight=26cm'/
     +        '\\parindent=0cm'/
     +        '\\begin{document}')
1100  format ('{\\Large Table of ',a3,' results for $\\ldots$}'/
     +        '\\vspace*{8mm}'/)
1110  format ('\\begin{tabular}{|c|c|c|l|} \\hline \\hline & & & \\\\'/
     + 'State & $T_e$ [eV] & $f_e$ & Leading excitations ($|c|\\geq ',
     +  f4.2,'$) \\\\'/' & & & \\\\ \\hline')
2010  format ('\\hline')
2140  format ('\\end{tabular}'/'\\newpage'/
     +        '\\begin{tabular}{|c|c|c|l|} \\hline \\hline & & & \\\\'/
     +        'State & $T_e$ $[eV]$ & $f_e$ & $\\Psi\\approx$ \\\\'/
     +        ' & & & \\\\ \\hline \\hline')
2125  format (' & & & $ ',a1,f4.2,' \\;(\\:',i3,a8,'\\,\\to\\,',
     +        i3,a8,'\\:) $ \\\\')
2120  format (' & & & $ ',a1,f4.2,' \\;(\\:',i3,a6,'\\,\\to\\,',
     +        i3,a6,'\\:) $ \\\\')
2111  format ('$ ',i3,'\\,^1\\!',a6,' $ & $ ',f5.2,' $ & $ ',
     +         f6.3,' $ & \\\\')
2115  format ('$ ',i3,'\\,^1\\!',a6,' $ & $ ',f5.2,' $ & $ ',
     +         f6.3,' $ & $ ',
     +            '\\hspace{0.8em} ',f4.2,' \\;(\\:',i3,a8,'\\,\\to\\,',
     +            i3,a8,'\\:) $ \\\\')
2110  format ('$ ',i3,'\\,^1\\!',a6,' $ & $ ',f5.2,' $ & $ ',
     +         f6.3,' $ & $ ',
     +         '\\hspace{0.8em} ',f4.2,' \\;(\\:',i3,a6,'\\,\\to\\,',
     +         i3,a6,'\\:) $ \\\\')
3010  format ('\\end{tabular}'/'\\end{document}')
_ELSE
1010  format ('\documentstyle[11pt]{article}'/
     +        '\pagestyle{empty}'/
     +        '\topmargin=-1.8cm'/
     +        '\oddsidemargin=-0.5cm'/
     +        '\evensidemargin=-0.5cm'/
     +        '\textwidth=17cm'/
     +        '\textheight=26cm'/
     +        '\parindent=0cm'/
     +        '\begin{document}')
1100  format ('{\Large Table of ',a3,' results for $\ldots$}'/
     +        '\vspace*{8mm}'/)
1110  format ('\begin{tabular}{|c|c|c|l|} \hline \hline & & & \\'/
     + 'State & $T_e$ [eV] & $f_e$ & Leading excitations ($|c|\geq ',
     +  f4.2,'$) \\'/' & & & \\ \hline')
2010  format ('\hline')
2140  format ('\end{tabular}'/'\newpage'/
     +        '\begin{tabular}{|c|c|c|l|} \hline \hline & & & \\'/
     +        'State & $T_e$ $[eV]$ & $f_e$ & $\Psi\approx$ \\'/
     +        ' & & & \\ \hline \hline')
2125  format (' & & & $ ',a1,f4.2,' \;(\:',i3,a8,'\,\to\,',
     +        i3,a8,'\:) $ \\')
2120  format (' & & & $ ',a1,f4.2,' \;(\:',i3,a6,'\,\to\,',
     +        i3,a6,'\:) $ \\')
2111  format ('$ ',i3,'\,^1\!',a6,' $ & $ ',f5.2,' $ & $ ',
     +         f6.3,' $ & \\')
2115  format ('$ ',i3,'\,^1\!',a6,' $ & $ ',f5.2,' $ & $ ',
     +         f6.3,' $ & $ ',
     +            '\hspace{0.8em} ',f4.2,' \;(\:',i3,a8,'\,\to\,',
     +            i3,a8,'\:) $ \\')
2110  format ('$ ',i3,'\,^1\!',a6,' $ & $ ',f5.2,' $ & $ ',
     +         f6.3,' $ & $ ',
     +         '\hspace{0.8em} ',f4.2,' \;(\:',i3,a6,'\,\to\,',
     +         i3,a6,'\:) $ \\')
3010  format ('\end{tabular}'/'\end{document}')
_ENDIF
      end
      block data rpa_data
      implicit none
INCLUDE(common/rpfile)
INCLUDE(common/polars)
      integer isym,isize,neig,iperm,isspac,
     +              maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),
     +              maxvec,mvec,mymaxv
      data itdspc,irpspc,itdtex,irptex,itable,inorms,ivect,idump,irest
     +  /  10,    11,    12,    13,    14,    15,    16,   17,   18   /
      data iperm/1,2,3,4,5,6,7,8,
     1           2,1,4,3,6,5,8,7,
     2           3,4,1,2,7,8,5,6,
     3           4,3,2,1,8,7,6,5,
     4           5,6,7,8,1,2,3,4,
     5           6,5,8,7,2,1,4,3,
     6           7,8,5,6,3,4,1,2,
     7           8,7,6,5,4,3,2,1/
      data ipola/1,2,3,2,4,5,3,5,6/
      end
      subroutine rpinit
c
c---------------------------------------------------------------------
c   Initiates the RPA calculation by reading in the directives
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit REAL (a-h,p-w), integer (i-n), logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*80 dmpfil,rstfil
INCLUDE(common/sizes)
      common/blkin /gin(511)
INCLUDE(common/blocko)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall,ioneel,jasym
INCLUDE(common/iofile)
INCLUDE(common/dump3)
      common/junk3 /ipair,icore(10),ivirt(10)
INCLUDE(common/polars)
INCLUDE(common/rpacom)
INCLUDE(common/rpfile)
      common/runlab/zcom(19),ztitle(10),zanam(maxat),zlabel(maxorb),
     +              ztagg(maxat),zsymm,zgroup,zscftp,zruntp,
     +              zguess,zconf,
     +              zstate,zorb0(maxorb),zpseud(maxat)
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),
     +              maxvec,mvec,mymaxv
      common/table /dmpfil,rstfil,rlines(2)
INCLUDE(common/work)
INCLUDE(common/workc)
INCLUDE(common/restrj)
INCLUDE(common/direc)
INCLUDE(common/machin)
      dimension yrpa(30)
      data yrpa/'symm','maxr','maxi','begi','pola',  !  1 -  5
     +          'vecr','anal','qz  ','thre','eige',  !  6 - 10
     +          'eqsy','tabl','end ','tda ','pref',  ! 11 - 15
     +          'use ','keep','debu','dump','rest',  ! 16 - 20
     +          'do  ','test','maxv','old ','star',  ! 21 - 25
     +          'nosk','----','----','----','jaro'/  ! 26 - 30
c
c   set up defaults
c
      modus = 3
      maxr = 50
      maxit = 30
      maxrr = maxr
      maxitt = maxit
      eps = 0.001d0
      epss = eps
      epsana = 0.01d0
      epstab = 0.2d0
      nbegin = 1
      oanal = .false.
      ouse = .false.
      othrsh = .false.
      oreduc = .true.
      if (odirpa) oreduc = .false.
      okeep = .false.
      nfreq = 0
      isspac(3) = 300000
      isspac(4) = isspac(3) / 20
      mymaxv = 0
      ipref = 7
      onew = .true.
      ldump = 0
      lrest = 0
      iorbit = 0
      igin = 0
      isold = 0
      irold = 0
      iblo13 = 0
      oskip = .true.
      ototal = .false.
      oall = .false.
      call setsto(24,0,nevlo)
      do 500 i=1,100
500   odebug(i) = .false.
      jasym = 0
c
c   read in directives
c
100   call input
      call inpa4 (ytext)
      i = locatc (yrpa,30,ytext)
      if (i .gt. 0) then
         go to ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     +        11,12,13,14,15,16,17,18,19,20,
     +          21,22,23,24,25,26,27,28,29,30),i
      else
         jrec = jrec - 1
_IF(rpagrad)
         go to 50
_ELSE
         ii = locatc(ydd(101),limit2,ytext)
         if (ii.ne.0) then
           go to 50
         else
            call caserr(
     +       'unrecognised directive or invalid directive sequence')
         end if
_ENDIF
      end if
c
c   symmetry and number of eigenvalues wanted
c
1     call inpi (nrsym)
      call inpi (nevlo(nrsym))
      call inpa4 (ytext)
      if (ytext .ne. 'to') call caserr('error in e-value specification')
      call inpi (nevhi(nrsym))
      go to 100
c
c   maximal size of reduced matrices for eigenvalue algorithm (maxr)
c                                  and linear equation solver (maxrr)
c
2     call inpi (maxr)
      if (jump .eq. 3) call inpi (maxrr)
      go to 100
c
c   maximal number of iterations for eigenvalue algorithm (maxit)
c                              and linear equation solver (maxitt)
c
3     call inpi (maxit)
      if (jump .eq. 3) call inpi (maxitt)
      go to 100
c
c   number of orbital where the symmetry counting begins
c
4     call inpi (nbegin)
      go to 100
c
c   polarizabilities
c
5     if (nfreq .gt. mxfreq) go to 100
      do 502 i=2,jump
      nfreq = nfreq + 1
      if (nfreq .le. mxfreq) go to 502
      write (iwr,501) mxfreq
501   format (/
     +  'WARNING: Dynamic polarizabilities will be evaluated only for',
     +  ' the first ',i6,' specified frequencies.')
      go to 100
502   call inpf (freq(nfreq))
      go to 100
c
c   vecrpa
c
6     call inpi (ionsec)
      go to 100
c
c   key word 'anal' for specification of analysis directive and
c   definition of the corresponding threshold
c
7     if (othrsh) then
         call inpf (epsana)
      else
         oanal = .true.
      end if
      go to 100
c
c   QZ algorithm (no reduction of the small eigenvalue problems)
c
8     oreduc = .false.
      go to 100
c
c   thresholds
c
9     othrsh = .true.
      go to 100
c
c   ... for iterative eigenvalue search
c
10    if (othrsh) then
         call inpf (eps)
      else
         call caserr ('keyword threshold was omitted')
      end if
      go to 100
c
c   ... for polarizability calculation
c
11    if (othrsh) then
         call inpf (epss)
      else
         call caserr ('keyword threshold was omitted')
      end if
      go to 100
c
c   ... for table of results
c
12    if (othrsh) then
         call inpf (epstab)
      else
         call caserr ('keyword threshold was omitted')
      end if
      go to 100
c
c   end of threshold block
c
13    othrsh = .false.
      go to 100
c
c   tda
c
14    if (jump .eq. 1) then
         modus = 2
      else
         call inpa4 (ytext)
         if (ytext .ne. 'only') go to 80
         modus = 1
      end if
      go to 100
c
c   prefactor exponent for integrals
c
15    call inpi (ipref)
      go to 100
c
c   use tda vectors as starting vectors for rpa calculation
c
16    call inpa4 (ytext)
      if (ytext .ne. 'tda') go to 80
      call inpa4 (ytext)
      if (ytext .ne. 'vect') go to 80
      ouse = .true.
      go to 100
c
c   keep rpa eigenvectors on fortran file
c
17    okeep = .true.
      go to 100
c
c   debugging options
c
18    call inpi (idebug)
      odebug(idebug) = .true.
      if (idebug .eq. 21) call inpi (jasym)
      go to 100
c
c   dump file
c
19    ii = istrt(2)
      do 190 i = 1 , inumb(2)
         dmpfil(i:i) = char1(ii:ii)
         ii = ii + 1
190   continue
      ldump = inumb(2)
      go to 100
c
c   restore file
c
20    call inpa4 (ytext)
      if (ytext .eq. 'tda') then
         mrest = 1
      else if (ytext .eq. 'rpa') then
         mrest = 2
      else
         go to 80
      end if
      ii = istrt(3)
      do 200 i = 1 , inumb(3)
         rstfil(i:i) = char1(ii:ii)
         ii = ii + 1
200   continue
      lrest = inumb(3)
      go to 100
c
c   compute total oscillator strength
c
21    call inpa4 (ytext)
      if (ytext .ne. 'tota') go to 80
      ototal = .true.
      if (jump .eq. 2) go to 100
      call inpa4 (ytext)
      if (ytext .ne. 'all') go to 80
      oall = .true.
      go to 100
c
c   u test
c
22    call inpi(mspec)
      call inpi(nspec)
      go to 100
c
c   maxvec
c
23    call inpi (mymaxv)
      go to 100
c
c   old
c
c   this directive invokes the old ubuild routine in direct rpa
c   which computes the fock matrices one by one. it is much 
c   slower than the default and should only be used for reasons 
c   of comparison.
c
24    onew = .false.
      go to 100
c
25    if (isold.eq.0 .and. irold.eq.0) then
         iunit = 14
         iblo13 = 1
         call vclr(gin,1,511)
      end if
      npairs = (jump - 3) / 3
      if (npairs*3+3 .ne. jump) 
     +   call caserr ('invalid number of fields in start directive')
      call inpi (isym)
      call inpi (iroot)
      if (isym.ne.isold .or. iroot.ne.irold) then
         igin = igin + 1
         gin(igin) = 0.0d0
         isold = isym
         irold = iroot
         if (igin .ge. 511) then
            call wrt3 (gin,511,iblo13,iunit)
            iblo13 = iblo13 + 1
            call vclr(gin,1,511)
            igin = 0
         end if
         igin = igin + 1
         gin(igin) = dble(isym)
         if (igin .ge. 511) then
            call wrt3 (gin,511,iblo13,iunit)
            iblo13 = iblo13 + 1
            call vclr(gin,1,511)
            igin = 0
         end if
         igin = igin + 1
         gin(igin) = dble(iroot)
         if (igin .ge. 511) then
            call wrt3 (gin,511,iblo13,iunit)
            iblo13 = iblo13 + 1
            call vclr(gin,1,511)
            igin = 0
         end if
      end if
      do 250 np = 1 , npairs
         call inpf (coeff)
         call inpi (ihole)
         call inpi (ipart)
         igin = igin + 1
         gin(igin) = coeff
         if (igin .ge. 511) then
            call wrt3 (gin,511,iblo13,iunit)
            iblo13 = iblo13 + 1
            call vclr(gin,1,511)
            igin = 0
         end if
         igin = igin + 1
         gin(igin) = dble(ihole)
         if (igin .ge. 511) then
            call wrt3 (gin,511,iblo13,iunit)
            iblo13 = iblo13 + 1
            call vclr(gin,1,511)
            igin = 0
         end if
         igin = igin + 1
         gin(igin) = dble(ipart)
         if (igin .ge. 511) then
            call wrt3 (gin,511,iblo13,iunit)
            iblo13 = iblo13 + 1
            call vclr(gin,1,511)
            igin = 0
         end if
250   continue
      go to 100
c
c  noskip
c
c  when this directive is specified, the product vectors of the
c  trial vectors in the first iteration of a restart run are
c  computed via routine rpmuld although this is not necessary.
c  the default is to reconstruct these vectors from the previous
c  product vectors. this directive was kept for reasons of
c  debugging (if necessary).
c
26    oskip = .false.
      go to 100
c
27    go to 100
28    go to 100
29    go to 100
c
c   jaroslav's small matrix
c
30    do 291 i = 2 , jump
         iorbit = iorbit + 1
         call inpi (iky(iorbit))
291   continue
      go to 100
c
c   return to routine "start"
c
50    iky(maxorb+1) = iorbit
      if (iblo13.ne.0 .and. igin.ne.0) call wrt3 (gin,511,iblo13,iunit)
      return
80    call caserr('unrecognised directive or faulty directive ordering')
      return 
      end
      subroutine rpmain (ipo,nbasis,dip,diag,q,qanal,nmaxly)
c
c---------------------------------------------------------------------
c   Performs Random Phase & Tamm-Dancoff calculation
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------
c
      implicit REAL (a-h,p-z), integer (i-n), logical (o)
INCLUDE(common/sizes)
      character*2 char2i
      character*4 name
      common/blkin /gin(511)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      parameter (nvcmax=9953)
      common/dipblk/iconv(nvcmax),index(nvcmax)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr
INCLUDE(common/polars)
INCLUDE(common/rpacom)
INCLUDE(common/rpfile)
_IF(rpagrad)
INCLUDE(common/sector)
_ENDIF
      common/symmet/isym,isize,neig,iperm(8,8),maxspa,intspa,i,need0
     + ipadsymmet(10)
      dimension diag(nvc),dip(3*nvc),q(*),ipo(nbasis,*),qanal(*)
_IF(rpagrad)
c
chvd  Variable for storing Epole (real and imaginary component)
c
      REAL epole(2)
c
chvd  The dumpfile section type for the RPA results section
chvd  (Probably it would be good if the section number could be
chvd   changed from within the input. For now we use a fixed number.)
c
      integer irpatyp, irpasec
      parameter(irpatyp=120)
      parameter(irpasec=351)
c
_ENDIF
      name(1:2) = 'tm'

c   calculate space for iteration and
c   adjust maximal size of reduced matrices if necessary

      maxspa = nmaxly
      maxscr = nmaxly - 200000

      iwarn = 0
      ineed0 = 4 * nvc + lenint(nbasis*nbasis)
50    maxr2 = maxr + maxr
      maxrs = maxr * maxr
      maxr2s = 4 * maxrs
      ineed = ineed0 + 2*maxrs + maxr + 3*maxr2s + 3*maxr2
      k = 0
      do 60 i = 1 , nirr
      kk = 2*icount(i)*(maxr+nevhi(i)+2) + nev(i)
60    k = max(k,kk)
      ineed = ineed + k
      if (ineed .gt. maxscr) then
         iwarn = 1
         maxr = maxr - 1
         go to 50
      end if
      if (iwarn .eq. 1) write (iwr,70) maxr
70    format (/
     +     'WARNING: maximal size of reduced matrices was chosen too',
     +     ' large'/9x,'and has now been adjusted to ',i4)
      intspa = maxspa - ineed
      write (iwr,80) maxspa,ineed,intspa
80    format (/'workspace statistics:'/21('-')/
     +  i14,' words of memory required'/
     +  i9,' words needed for iterative algorithm'/
     +  i14,' words left as scratch space for multiplication routines'/)

c   collect rpa integrals and construct A & B matrices

      call rpsort (q,ipo,nbasis,dip)
      call vclr(q,1,maxspa-ineed0)
_IF(rpagrad)
c
chvd  Set up dumpfile I/O section to store the RPA results.
chvd  The results stored are ipo, Z, Y, re(E), im(E)
chvd  where ipo is stored only once and the other data are stored
chvd  for every state.
c
      irpalen = 0
      do isym = 1, nirr
         irpalen = irpalen
     +           + nev(isym)*(2*lensec(icount(isym))+lensec(2))
     +           + lensec(lenint(nbasis*nbasis))
      enddo
      call secput(irpasec,irpatyp,irpalen,irpablk)
      call wrt3i(ipo,nbasis*nbasis,irpablk,numdu)
_ENDIF
c
c   perform symmetry blocked random phase calculation
c
      do 1000 isym = 1 , nirr
_IF(rpagrad)
c
chvd     Set nrpalen to the correct value. This is needed to keep the
chvd     dumpfile in sync in case a set of eigenvectors didn't converge.
c
         nrpalen = lensec(lenint(nbasis*nbasis))
         do jsym = 1, isym-1
            nrpalen = nrpalen
     +              +nev(jsym)*(2*lensec(icount(jsym))+lensec(2))
         enddo
c
_ENDIF
         isize = icount(isym)
         neig = nev(isym)
         if (neig .gt. 0) write(iwr,120) isym,isize
120      format (///44('-')/'dimension of rpa matrix for symmetry ',
     +           i1,': ',i4/44('-')/)
         call flushn(iwr)
         if (isize .eq. 0) go to 1000

c   tm file

         itmfil = 40 + isym
         name(3:4) = char2i(itmfil)
         if (opg_root()) then
            open (itmfil,file=name,form='unformatted',status='unknown')
         endif
         i1 = 1
         i2 = i1 + nbasq
         call tmfil2 (q(i1),q(i2))
c
         nstart = nevhi(isym)
         nevlo1 = nevlo(isym) - 1
         nvcmxr = isize * maxr
         nvcnev = isize * neig
         nvcnst = isize * nstart
         isize2 = isize + isize
         isize3 = isize + isize2
c
c   define pointers for partitioning of q
c
         ib1    = 1
         ib2    = ib1 + nvcmxr
         isred1 = ib2 + nvcmxr
         iered1 = isred1 + maxrs
         ibred  = iered1 + maxrs
         isred2 = ibred + maxr
         iered2 = isred2 + maxr2s
         ietare = iered2 + maxr2s
         ietaim = ietare + maxr2
         isigma = ietaim + maxr2
         ixred  = isigma + maxr2
         ih1    = ixred + maxr2s
         ih2    = ih1 + isize
         iv1    = ih2 + isize
         iv2    = iv1 + nvcnst
         ires   = iv2 + nvcnst
         iscr   = ires + neig
c
c   get diagonal of matrix A for iterative algorithm
c
         call getdia (.false.,ipo,nbasis,diag,index,q)
         call vclr(q,1,maxspa-ineed0)
         if (neig .eq. 0) then
            if (oall) then
               go to 660
            else
               go to 900
            end if
         end if
c
c   choose unit vectors e(INDEX(1)) to e(INDEX(NSTART)) 
c   as starting vectors
c
         call setsto(nvcmax,0,iconv)
         call vclr(tim,1,4)
         k = ib1 - 1
         do 1300 i = 1 , nstart
            q(k+index(i)) = 1.0d0
            k = k + isize
1300     continue
c
c   special starting vectors ?
c
         if (iblo13 .eq. 0) go to 790
         iunit = 14
         jblo13 = 0
1200     jblo13 = jblo13 + 1
         call rdedx (gin,511,jblo13,iunit)
         inew = 0
         do 1210 m = 1 , 511
            if (gin(m) .eq. 0.0d0) then
               inew = 1
               go to 1210
            end if
            if (inew .eq. 1) then
               jsym = dnint(gin(m))
               inew = 2
               go to 1210
            end if
            if (inew .eq. 2) then
               jroot = dnint(gin(m))
               inew = 3
               if (jsym .eq. isym) 
     +         call vclr(q(ib1+(jroot-1)*isize),1,isize)
               go to 1210
            end if
            if (inew .eq. 3) then
               inew = 4
               if (jsym .ne. isym) go to 1210
               temp = gin(m)
               go to 1210
            end if
            if (inew .eq. 4) then
               inew = 5
               if (jsym .ne. isym) go to 1210
               iitemp = dnint(gin(m))
               go to 1210
            end if
            if (inew .eq. 5) then
               inew = 3
               if (jsym .ne. isym) go to 1210
               iatemp = dnint(gin(m))
               k = ipo(iitemp,iatemp)
               q(ib1-1+(jroot-1)*isize+k) = temp
               go to 1210
            end if
1210     continue
         if (jblo13 .lt. iblo13) go to 1200
         if (odebug(74)) then
            do 1250 j = 1 , nstart
               write (iwr,1220) isym,j
1220       format ('symmetry ',i1,', special starting vector ',i3,':')
               call outvec (q(ib1+(j-1)*isize),isize,' ')
1250        continue
         end if

c   orthonormalize special starting vectors

         ibnew1 = ib1
         do 1280 j = 1 , nstart
            ibrun1 = ib1
            do 1270 k = 1 , j-1
            t = - ddot(isize,q(ibnew1),1,q(ibrun1),1)
            call daxpy(isize,t,q(ibrun1),1,q(ibnew1),1)
            ibrun1 = ibrun1 + isize
1270        continue
            t = dnrm2(isize,q(ibnew1),1)
            if (1.0d0+t .eq. 1.0d0) then
               call caserr ('special trial vector must be discarded.')
            end if
            t = 1.0d0 / t
            call dscal(isize,t,q(ibnew1),1)
            ibnew1 = ibnew1 + isize
1280     continue

790      if (modus .eq. 3) go to 139

c   perform TDA calculation

         itda = 0
         niter = 0
         write (iwr,131) 'TDA'
131      format (///'iterative ',a3,' procedure:'/24('*')//)
         call vfill(1.0d0,q(ietaim),1,neig)
         tim(1) = cpulft(1)
         ierr = 0
         if(opg_root()) then
          rewind itable
          rewind inorms
         endif
1309     call tdalgo (diag,q(ib1),q(iered1),q(ietare),
     +                q(ixred),q(ih1),q(iv1),q(ires),q(iscr))
         if (ierr .ne. 0) go to 670
         nconv = isum(neig,iconv,1)
         if (nconv .eq. neig) go to 200
         if (niter .lt. maxit) go to 1311
         write (iwr,1310)
1310     format 
     +    (/'maximum number of iterations reached, no convergence.'/)
         go to 670
c
c   restart iteration process with lowest eigenvectors of last iteration
c
1311     call dcopy(nvcnst,q(iv1),1,q(ib1),1)
         nstar1 = nstart
         k = 1
         ib = ib1
1320     ibrun = ib1
         do 1330 j = 1 , k-1
            call orth (isize,q(ib),q(ibrun))
            ibrun = ibrun + isize
1330     continue
         t = dnrm2(isize,q(ib),1)
         if (1.0d0+t .eq. 1.0d0) then
            write (iwr,1340)
1340        format 
     + ('new trial vector must be discarded in restart run of tdalgo.')
            nstar1 = nstar1 - 1
            ii = nstar1 * isize
            call dswap(isize,q(ib),1,q(ib1+ii),1)
            go to 1320
         else
            call dscal(isize,1.0d0/t,q(ib),1)
            k = k + 1
            ib = ib + isize
         end if
         if (k .gt. nstar1) go to 1309
         go to 1320

c   perform RPA calculation

139      itda = 1
         ioff = 0
         niter = 0
         write (iwr,131) 'RPA'
         tim(1) = cpulft(1)
         ierr = 0
         if(opg_root())then
          rewind itable
          rewind inorms
         endif
150      call rpalgo (diag,q(ib1),q(ib2),q(isred1),q(iered1),q(ibred),
     +           q(isred2),q(iered2),q(ietare),q(ietaim),q(isigma),
     +           q(ixred),q(ih1),q(ih2),q(iv1),q(iv2),q(ires),q(iscr))
         if (ierr .ne. 0) go to 670
         nconv = isum(neig,iconv,1)
         if (nconv .eq. neig) go to 200
         if (niter .lt. maxit) go to 151
         write (iwr,1310)
         go to 670

c   restart iteration process with lowest eigenvectors of last iteration

151      call dcopy(nvcnst,q(iv1),1,q(ib1),1)
         call dcopy(nvcnst,q(iv2),1,q(ib2),1)
         nstar1 = nstart
         k = 1
         ibnew1 = ib1
         ibnew2 = ib2
170      ibrun1 = ib1
         ibrun2 = ib2
         do 180 j = 1 , k-1
            t = - ddot(isize,q(ibnew1),1,q(ibrun1),1)
     +          - ddot(isize,q(ibnew2),1,q(ibrun2),1)
            call daxpy(isize,t,q(ibrun1),1,q(ibnew1),1)
            call daxpy(isize,t,q(ibrun2),1,q(ibnew2),1)
            t = - ddot(isize,q(ibnew1),1,q(ibrun2),1)
     +          - ddot(isize,q(ibnew2),1,q(ibrun1),1)
            call daxpy(isize,t,q(ibrun2),1,q(ibnew1),1)
            call daxpy(isize,t,q(ibrun1),1,q(ibnew2),1)
            ibrun1 = ibrun1 + isize
            ibrun2 = ibrun2 + isize
180      continue
         t = dsqrt(dnrsq(isize,q(ibnew1),1)+dnrsq(isize,q(ibnew2),1))
         if (1.0d0+t .eq. 1.0d0) then
            write (iwr,190)
190         format 
     +  ('new trial vector must be discarded in restart run of rpalgo.')
            nstar1 = nstar1 - 1
            ii = nstar1 * isize
            call dswap(isize,q(ibnew1),1,q(ib1+ii),1)
            call dswap(isize,q(ibnew2),1,q(ib2+ii),1)
            go to 170
         else
            t = 1.0d0 / t
            call dscal(isize,t,q(ibnew1),1)
            call dscal(isize,t,q(ibnew2),1)
            k = k + 1
            ibnew1 = ibnew1 + isize
            ibnew2 = ibnew2 + isize
         end if
         if (k .gt. nstar1) go to 150
         go to 170
200      tim(4) = cpulft(1) - tim(1)
_IF(rpagrad)
c
chvd     The current eigenvectors are converged so quickly store the RPA
chvd     results on the DUMPFILE
c
         do i = 0, nev(isym)-1
            epole(1) = q(ietare+i+nevlo1)
            epole(2) = q(ietaim+i+nevlo1)
            call wrt3(epole,2,irpablk+nrpalen,numdu)
            nrpalen = nrpalen+lensec(2)
            call wrt3(q(iv1+(i+nevlo1)*isize),isize,irpablk+nrpalen,
     +                numdu)
            nrpalen = nrpalen+lensec(isize)
            call wrt3(q(iv2+(i+nevlo1)*isize),isize,irpablk+nrpalen,
     +                numdu)
            nrpalen = nrpalen+lensec(isize)
         enddo
_ENDIF
c
c   output of convergence behaviour
c
         write (iwr,132)
132      format ('iteration  reduced    eigenvalues'/
     +           'number     dimension'/89('-'))
         call outeig (itable,iwr,dummy,neig,niter,ntot,nevlo1,2)
         if (niter .eq. 1) write (iwr,100)
         if (niter .gt. 1) write (iwr,201) niter
100      format (//89('*')/
     +   '*',10x,'all eigenvalues have converged after one iteration.',
     +       26x,'*'/89('*')//)
201      format (//89('*')/'*',
     +             10x,'all eigenvalues have converged after ',i2, 
     +             ' iterations.',26x,'*'/89('*')//)

         write (iwr,133)
133      format (/'iteration             norms of residue vectors'
     +           /'number'/80('-'))
         call outres (iwr,inorms,dummy,nevlo1,neig,niter,2)

c   write table of results
 
         ie = ietare + nevlo1 + ioff
         ishft = nevlo1 * isize
         otda = itda.eq.0
         oincor = .true.
         call rpanal (ipo,nbasis,q(ie),q(iv1+ishft),q(iv2+ishft),
     +                q(ietaim),q(isigma),dip(iofsym(isym)+1),
     +                dummy,iconv,otda,oincor,2,qanal)
c
c   optionally write eigenvectors on fortran file 'vectors'
c
         if (okeep) then
            write (ivect) isym,neig
            do 650 j = 1 , neig
               write (ivect) (q(iv1+(j-1)*isize+k-1),k=1,isize)
               write (ivect) (q(iv2+(j-1)*isize+k-1),k=1,isize)
650         continue
         end if
c
c   timing analysis
c
670      write (iwr,205) tim(2),tim(3),tim(4)-tim(2)-tim(3),tim(4)
205      format (/'timing analysis:'/
     +   'multiplication of E with trial vectors: ',f14.4,' seconds'/
     +   'solution of reduced eigenvalue problems:',f14.4,' seconds'/
     +   'remaining tasks:                        ',f14.4,' seconds'/
     +   'total iterative procedure:              ',f14.4,' seconds')
c
c   restart ?
c
         if (modus .eq. 1) go to 660
         if (modus .eq .2 .and. itda .eq. 0 .and. .not. ouse) then
            call vclr(q,1,maxspa-ineed0)
            call setsto(nvcmax,0,iconv)
            call vclr(tim,1,4)
            k = ib1 - 1
            do 375 i = 1 , nstart
               q(k+index(i)) = 1.0d0
               k = k + isize
375         continue
            go to 139
         end if
         if (modus .eq. 2 .and. itda .eq. 0 .and. ouse) then
            call dcopy(nvcnst,q(iv1),1,q(ib1),1)
            call vclr(q(ib1+nvcnst),1,maxspa-ineed0-nvcnst)
            call setsto(nvcmax,0,iconv)
            call vclr(tim,1,4)
            go to 139
         end if
c
c   compute total oscillator strength (see C. Fuchs, PhD thesis)
c
660      if (ototal) then
            imy = iofsym(isym) + 1
            do 940 icoord = 1 , 3
               if (dnrm2(isize,dip(imy),1) .gt. tolmy) then
                  if (modus .le. 2) then
                     call rpmult (.false.,jstart(isym),npair(isym),
     +       isize,num8,dip(imy),dip(imy),q(ih1),q(ih2),q(iscr))
                     tottda = tottda + 
     +               ddot(isize,dip(imy),1,q(ih1),1)
                  end if
                  if (modus .ge. 2) then
                     if (odebug(75)) 
     +               call outvec (dip(imy),isize,'vector my')
                     call rpmult (.true.,jstart(isym),npair(isym),
     +       isize,num8,dip(imy),dip(imy),q(ih1),q(ih2),q(iscr))
                     totrpa = totrpa 
     +               + ddot(isize,dip(imy),1,q(ih1),1)
     +               + ddot(isize,dip(imy),1,q(ih2),1)
                  end if
               end if
               imy = imy + nvc
940         continue
         end if
900      continue
         if (nfreq .eq. 0) go to 999
c
c   polarizability calculation
c
         write (iwr,910)
910      format (//122('=')//)
c
         imy    = iofsym(isym) + 1
         ib1    = imy + isize3
         ib2    = ib1 + nvcmxr
         isred  = ib2 + nvcmxr
         iered  = isred + maxrs
         ibred  = iered + maxrs
         iared  = ibred + maxr
         ifred1 = iared + maxr2s
         ifred2 = ifred1 + maxr
         ih1    = ifred2 + maxr2
         ih2    = ih1 + isize
         ic1    = ih2 + isize
         ic2    = ic1 + isize
         iscr   = ic2 + isize

         imy1 = imy
         do 850 icoor1 = 1 , 3
            if (dnrm2(isize,dip(imy1),1) .le. tolmy) go to 845
            do 840 i = 1 , nfreq
c
c   starting vector
c
_IFN1(civu)               call vsadd(diag,1,freq(i),q(ih1),1,isize)
_IF1(civu)               call vsav(isize,freq(i),q(ih1),diag)
               call vdivz (dip(imy1),1,q(ih1),1,1.0d0,q(ib1),1,isize)
               if (freq(i) .eq. 0.0d0) then
                  t = dnrm2(isize,q(ib1),1)
                  if (1.0d0+t .eq. 1.0d0) 
     +            call caserr ('vector f equals zero !')
                  call dscal(isize,1.0d0/t,q(ib1),1)
                  call vclr(q(ib2),1,isize)
               else
_IFN1(civu)               call vsadd(diag,1,-freq(i),q(ih2),1,isize)
_IF1(civu)               call vsav(isize,-freq(i),q(ih2),diag)
                  call vdivz (dip(imy1),1,q(ih2),1,1.0d0,q(ib2),1,isize)
                  call vneg (q(ib2),1,q(ib2),1,isize)
                  call symort (isize,q(ib1),q(ib2),sytol,oldep)
                  if (oldep) call caserr ('oldep = .true. !')
               end if
c
c   solve linear response equation system
c
               write (iwr,823) isym,icoor1,freq(i)
823            format (/'Polarizability calculation for symmetry',i2,
     +             ', coordinate',i2,
     +             ', frequency',f9.4,' a.u.'/80('-')//6x,
     +             'iteration   reduced      residue'/18x,'dimension')
               niter = 0
               iconv(1) = 0
               ierr = 0
824            call rpsolv (diag,q(ib1),q(ib2),freq(i),q(isred),
     +            q(iered),q(ibred),
     +            q(iared),dip(imy1),q(ifred1),q(ifred2),
     +            q(ih1),q(ih2),
     +            q(ic1),q(ic2),q(iscr))
               if (ierr .ne. 0) go to 840
               if (iconv(1) .eq. 1) go to 825
               if (niter .gt. maxitt) go to 8245
               call dcopy(isize,q(ic1),1,q(ib1),1)
               if (freq(i) .eq. 0.0d0) then
                  t = dnrm2(isize,q(ib1),1)
                  if (1.0d0+t .eq. 1.0d0) 
     +             call caserr ('restart trial vector equals zero !')
                  call dscal(isize,1.0d0/t,q(ib1),1)
                  call vclr(q(ib2),1,isize)
               else
                  call dcopy(isize,q(ic2),1,q(ib2),1)
                  call symort (isize,q(ib1),q(ib2),sytol,oldep)
                  if (oldep) call caserr 
     +            ('symort problems in restart run of rpsolv')
               end if
               go to 824
8245           write (iwr,1310)
               go to 840
825            imy2 = imy1
c
c   contribution to polarizability tensors
c
               do 830 icoor2 = icoor1 , 3
                  k = ipola(icoor1,icoor2)
                  pola(k,i) = pola(k,i) - 2.0d0 * 
     +            (ddot(isize,q(ic1),1,dip(imy2),1)
     +          -  ddot(isize,q(ic2),1,dip(imy2),1))
                  imy2 = imy2 + nvc
830            continue
c
840         continue
845         imy1 = imy1 + nvc
850      continue
c
c   end loop over symmetries
c
999      if (opg_root()) then
          close (itmfil)
         endif
1000  continue
c
      return
      end
      subroutine rpmult (orpa,iblo,npair,isize,ned,y,z,v,w,q)

c----------------------------------------------------------------------
c   Computes the matrix-vector product 
c
c               | v |     |  A   B  | | y |
c               |   |  =  |         | |   |
c               | w |     |  B   A  | | z |
c
c   (or simply  v = Ay  if orpa = .false.).
c   A,B are read in from external file NED, starting at block IBLO
c   ISIZE is dimension of A and B
c   NPAIR is the number of rows of A,B that can be held in 
c   temporary storage
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------

      implicit REAL (a-h,p-z), logical (o)
      dimension y(*),z(*),v(*),w(*),q(*)
      call vclr(v,1,isize)
      if (orpa) call vclr(w,1,isize)
      jblo = iblo
      ntot = (isize-1) / npair + 1
      np = npair
      ipair = 0
      do 400 ibox = 1 , ntot
         if (ibox .eq. ntot) np = isize - ipair
         npp = np * isize
         npp2 = npp + npp
         if (orpa) then
            call rdedx (q,npp2,jblo,ned)
         else
            call rdedx_less(q,npp ,jblo,ned)
         end if
         n0 = 0
         do 300 n = 1 , np
            if (orpa) then
               n1 = n0 + npp
               do 100 i = 1 , isize
               v(ipair+n) = v(ipair+n) + q(n0+i) * y(i) + q(n1+i) * z(i)
               w(ipair+n) = w(ipair+n) + q(n0+i) * z(i) + q(n1+i) * y(i)
100            continue
            else
               do 200 i = 1 , isize
                  v(ipair+n) = v(ipair+n) + q(n0+i) * y(i)
200            continue
            end if
            n0 = n0 + isize
300      continue
         jblo = jblo + lensec(npp2)
         ipair = ipair + np
400   continue
      return
      end
      subroutine rpsolv (diag,b1,b2,freq,sred,ered,bred,ared,f,fred1,
     +         fred2,h1,h2,c1,c2,scr)
c
c----------------------------------------------------------------------
c   Solves the linear response equation system (omega * S - E) c = f
c   using the iterative algorithm by Olsen & Joergensen
c   (c) Carsten Fuchs 1991-1993
c----------------------------------------------------------------------
c
      implicit REAL (a-h,p-z), integer (i-n), logical (o)
INCLUDE(common/sizes)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100)
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      parameter (nvcmax=9953)
      common/dipblk/iconv(nvcmax),index(nvcmax)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot,ierr
INCLUDE(common/rpacom)
      common/symmet/isym,isize,neig,ipadsymmet(77)
      dimension diag(isize),b1(isize,maxr),b2(isize,maxr)
      dimension sred(maxr,maxr),ered(maxr,maxr),bred(maxr)
      dimension ared(maxr2,maxr2),f(isize),fred1(maxr),fred2(maxr2)
      dimension h1(isize),h2(isize),c1(isize),c2(isize),scr(*)

      time = cpulft(1)
      ntot = 1
10    niter = niter + 1
      if (niter .gt. maxitt) return

c   calculate the vectors e*b(i) and fill up matrices sred,ered

      fred1(ntot) = ddot(isize,b1(1,ntot),1,f,1) - 
     +              ddot(isize,b2(1,ntot),1,f,1)
      call rpmult (.true.,jstart(isym),npair(isym),isize,num8,
     +   b1(1,ntot),b2(1,ntot),h1,h2,scr)
      bred(ntot) = ddot(isize,b1(1,ntot),1,h2,1) + 
     +             ddot(isize,b2(1,ntot),1,h1,1)
      do 20 j=1,ntot
      ered(j,ntot) = ddot(isize,b1(1,j),1,h1,1)
     +   + ddot(isize,b2(1,j),1,h2,1)
      sred(j,ntot) = ddot(isize,b1(1,ntot),1,b1(1,j),1)
     +   - ddot(isize,b2(1,ntot),1,b2(1,j),1)
      if (j .eq. ntot) go to 20
      ered(ntot,j) = ddot(isize,b1(1,j),1,h2,1)
     +   + ddot(isize,b2(1,j),1,h1,1)
      sred(ntot,j) = ddot(isize,b1(1,ntot),1,b2(1,j),1)
     +   - ddot(isize,b2(1,ntot),1,b1(1,j),1)
20    continue

c   construct and solve reduced equation system

      do 200 i=1,ntot
      fred2(i)      =  fred1(i)
200   fred2(ntot+i) = -fred1(i)
      call redeqs (maxr,ntot,freq,sred,ered,bred,ared,fred2,scr,ierr)
      if (ierr .ne. 0) go to 400

c   expand solution fred2

      do 40 i=1,isize
      c1(i) = 0.0d0
      c2(i) = 0.0d0
      do 40 k=1,ntot
      c1(i) = c1(i) + b1(i,k) * fred2(k) + b2(i,k) * fred2(k+ntot)
      c2(i) = c2(i) + b2(i,k) * fred2(k) + b1(i,k) * fred2(k+ntot)
40    continue

c   compute residue vector

      call rpmult (.true.,jstart(isym),npair(isym),isize,num8,
     +   c1,c2,h1,h2,scr)
      do 50 i=1,isize
      h1(i) =   freq * c1(i) - h1(i) - f(i)
50    h2(i) = - freq * c2(i) - h2(i) + f(i)

c   compute euclidean norm of residue vector

      t = dsqrt(dnrsq(isize,h1,1)+dnrsq(isize,h2,1))
      write (iwr,51) niter,ntot,t
51    format (2(8x,i3),6x,f14.9)

c   threshold ?

      if (t .le. epss) go to 110
      if (ntot.ge.maxr .or. ntot.ge.isize) return
      ntot = ntot + 1

c   update residue vector 

      do 55 i = 1 , isize
         b1(i,ntot) =  freq - diag(i)
55    continue
      if (freq .eq. 0.0d0) go to 57
      do 56 i = 1 , isize
         b2(i,ntot) = -freq - diag(i)
56    continue
57    call vdivz (h1,1,b1(1,ntot),1,1.0d0,b1(1,ntot),1,isize)
      if (freq .ne .0.0d0) 
     + call vdivz (h2,1,b2(1,ntot),1,1.0d0,b2(1,ntot),1,isize)

c   gram-schmidt orthogonalisation

      if (freq .eq. 0.0d0) go to 90
      do 60 k=1,ntot-1
      t = - ddot(isize,b1(1,ntot),1,b1(1,k),1)
     +    - ddot(isize,b2(1,ntot),1,b2(1,k),1)
      call daxpy(isize,t,b1(1,k),1,b1(1,ntot),1)
      call daxpy(isize,t,b2(1,k),1,b2(1,ntot),1)
      t = - ddot(isize,b1(1,ntot),1,b2(1,k),1)
     +  - ddot(isize,b2(1,ntot),1,b1(1,k),1)
      call daxpy(isize,t,b2(1,k),1,b1(1,ntot),1)
60    call daxpy(isize,t,b1(1,k),1,b2(1,ntot),1)
      t = dsqrt(dnrsq(isize,b1(1,ntot),1)+dnrsq(isize,b2(1,ntot),1))
      if (1.0d0+t .eq. 1.0d0) then
         write (iwr,70) niter
70       format (
     +  'unable to finish iterative solution after iteration step ',i2/
     +  'since new trial vector is linear dependent.')
         ierr = 1
         return
      else
c
c   perform symmetric orthonormalisation
c
         call symort (isize,b1(1,ntot),b2(1,ntot),sytol,oldep)
         if (oldep) then
            write (iwr,80) niter
80          format (
     +      'unable to finish iterative solution after iteration step ',
     +      i2/'because b and b# are (nearly) linearly dependent.')
            ierr = 1
            return
         end if
      end if
      go to 10
90    continue
      do 100 k=1,ntot-1
100   call orth (isize,b1(1,ntot),b1(1,k))
      t = dnrm2(isize,b1(1,ntot),1)
      if (1.0d0+t .eq. 1.0d0) then
         write (iwr,70) niter
         ierr = 1
         return
      else
         call dscal(isize,1.0d0/t,b1(1,ntot),1)
         call vclr(b2(1,ntot),1,isize)
      end if
      go to 10

c   convergence reached

110   iconv(1) = 1
      if (niter .eq. 1) write (iwr,280) cpulft(1)-time
      if (niter .gt .1) write (iwr,290) niter,cpulft(1)-time
280   format (/'convergence reached after one iteration in ',f9.4,
     +         ' seconds.'/)
290   format (/'convergence reached after',i3,' iterations in ',f9.4,
     +         ' seconds.'/)
      return
c
c   error message
c
400   write (iwr,410) niter,ierr
410   format (/'iteration step ',i2,
     +         ': condition number too small: ierr = ',i3/)
c
      return
      end
      subroutine rpsort (q,ipo,nbasis,dip)

c   Collects and reorders the RPA integrals (ai|bj),(ab|ij) and, in the
c   end, constructs the RPA matrices A and B

      implicit REAL (a-h,p-z), integer (i-n), logical (o)
_IFN1(ck)      integer*2 itemp(4)
      logical btest
_IF1(k)      integer itemp(4)
      character*8 title
      character*3 symb(16)
INCLUDE(common/sizes)
      common/blkin /gin(511)
      common/blkin1/evalue(maxorb),hocc(maxorb),etot,nbas2,
     +            newbas,ncol,ivalue,ioccup,ispa,inx(maxorb),
     +            ipoint(maxorb),isit(maxorb),itogrp(maxat),
     +            ilabel(maxorb),iok
INCLUDE(common/blocko)
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall,ioneel,jasym
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      parameter (nvcmax=9953)
      common/dipblk/ifill(nvcmax),index(nvcmax)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/atmblk)
      common/junk2 /i205(4,340)
      common/junko /cspace(30),irspa(700),irspl(40),irspb(1590),
     +              irspc(336),irspd(8)
      common/restri/nfils(63),lda(508),isect(508),ldx(508)
      common/polars/nfreq
      common/restar/nprint(142),n6file,n6tape(20),n6blk(20),n6last(20)
      common/restrz/title(50)
INCLUDE(common/rpacom)
      common/symmet/isym,isize,neig,iperm(8,8),maxspa,intspa,ineed0,
     + ipadsymmet(10)
      common/scrtch/a(200,200),eta(200),h(200),x(200,200),s(200,200),
     +              etaim(200),sigma(200)
c 
      dimension q(*),l6blk(20),iad(980),ibuf(980),jbuf(980)
      dimension dip(*),ipo(nbasis,*)
      dimension iipt(100),iapt(100)
c
_IFN1(ck)      equivalence (pack,itemp(1))
      equivalence (mword,gin(511))
_IF(linux)
      external fget
c
_ENDIF
      data m21/21/
      data symb/'a''','a"','a1','b1','b2','a2',
     +    'ag','b1u','b2u','b3g','b3u','b2g','b1g','au','TDA','RPA'/
_IF(cray,t3d)
      data mask /37777777777b/
_ENDIF
c 
      time = cpulft(1)
      ned4 = 5
      call secget (isect(501),m21,ibl3rs)
      call rdchr (title,ldx(isect(501)),ibl3rs,idaf)
      call reads (cspace,lda(isect(501)),idaf)
      call icopy (203,irspa,1,nprint,1)
      do 10 i = 1 , n6file
         l6blk(i) = n6blk(i) - n6last(i)
10    continue
      call setsto(8,0,npair)
      call setsto(16,0,nbox)
      call setsto(980,0,ibuf)
      call setsto(980,0,jbuf)
      call vclr(q,1,intspa)
      do 11 i = 1 , nirr
         if (icount(i) .gt. 0) then
            if (icount(i) .gt. nvcmax) then
               write (iwr,1001) i
1001         format ('too many excitation pairs for symmetry ',i1//
     +        'this error can easily be fixed -- contact author'//)
               call caserr ('error in rpa program')
            end if
            npair(i) = intspa / (3*icount(i))
            nbox(i) = (icount(i)-1) / npair(i) + 1
         end if
11    continue
c 
      nbox0(1) = 1
      do 12 i = 2 , nirr
         nbox0(i) = nbox0(i-1) + nbox(i-1)
12    continue
c
      ntot = isum(nirr,nbox,1)
      if (ntot*511 .gt. maxspa-ineed0) then
         write (iwr,13)
13       format ('insufficient space for rpa integral sorting')
         call caserr ('error in rpa program')
      end if
      write (iwr,133)
133   format (122('=')//'rpsort parameters:'//
     +  'irrep   size of matrices A,B   no of pairs that can be',
     +  '    no of boxes'/
     +  '                               held in scratch space'/
     +  '----------------------------------------------------',
     +  '-----------------')
      do 134 i = 1 , nirr
         write (iwr,135) i,icount(i),npair(i),nbox(i)
134   continue
135   format (2x,i1,6x,i10,14x,i10,15x,i5)
      write (iwr,136)
136   format (//)

      iorbit = iky(maxorb+1)
      if (iorbit .gt. 0 .and. odebug(21)) then
         ndim = 0
         irun = 0
800      irun = irun + 1
         if (irun .gt. iorbit) go to 810
         ndim = ndim + 1
         iipt(ndim) = iky(irun)
         irun = irun + 1
         iapt(ndim) = iky(irun)
         go to 800
810      continue
         if (ndim .gt. 100) call caserr 
     +                     ('small rpa matrix is too big !')
      end if
c
      call setsto(ntot,0,ifill)
      call setsto(ntot,0,iad)
      iblo4 = 1
      mtot = 0
      mtorpa = 0
      nconta = 0
c
c   now read in integral file ed6
c
      call setsto(1360,0,i205)
c
      do 100 ifile=1,n6file
      ned = n6tape(ifile)
      call search (n6blk(ifile),ned)
      call find (ned)
      jj = l6blk(ifile)
20    jj = jj + 1
_IF1(c)      call get (gin,nw)
_IFN1(c)      call fget (gin,nw,ned)
      if (nw .le. 0) go to 100
      call unpack(gin(num2e+1),lab816,i205,numlab)
      mtot = mtot + mword
      if (odebug(71)) print*,'mword = ',mword
_IF1(c)      if (jj.ne.0) call find(ned)
      do 60 irun=1,mword
_IF(littleendian)
      indxj = i205(1,irun)
      indxi = i205(2,irun)
      indxl = i205(3,irun)
      indxk = i205(4,irun)
_ELSE
      indxi = i205(1,irun)
      indxj = i205(2,irun)
      indxk = i205(3,irun)
      indxl = i205(4,irun)
_ENDIF
      if (odebug(71)) write (iwr,33) 
     +   indxi,indxj,indxk,indxl,gin(irun)
33    format (4i3,f14.9)
      if (iperm(isymmo(indxi),isymmo(indxj)) .ne. 
     +    iperm(isymmo(indxk),isymmo(indxl))) 
     +    nconta = nconta + 1
      if (indxi.le.na .or.  indxl.gt.na) go to 60
      if (indxj.gt.na .and. indxk.gt.na) go to 60
      if (indxj.le.na .and. indxk.le.na) go to 60
      mtorpa = mtorpa + 1
      if (indxj .gt. na) go to 40

c   cases (ai|bj), (bj|ai)
c          --          --

      isym = iperm(isymmo(indxi),isymmo(indxj))
      if (iperm(isymmo(indxk),isymmo(indxl)) .ne. isym) 
     +        go to 60
      ipos = ipo(indxi,indxj)
      jpos = ipo(indxk,indxl)
      call inbuck (iad,iblo4,irun,isym,ipos,jpos,1,q)
      call inbuck (iad,iblo4,irun,isym,jpos,ipos,1,q)

c   cases (aj|bi), (bi|aj)
c          -   -     - -

      isym = iperm(isymmo(indxi),isymmo(indxl))
      if (iperm(isymmo(indxj),isymmo(indxk)) .ne. isym) 
     +        go to 60
      ipos = ipo(indxi,indxl)
      jpos = ipo(indxk,indxj)
      call inbuck (iad,iblo4,irun,isym,ipos,jpos,3,q)
      call inbuck (iad,iblo4,irun,isym,jpos,ipos,3,q)
      go to 60

c   cases (ab|ij), (ba|ji)
c          -  -      -  -

40    isym = iperm(isymmo(indxi),isymmo(indxk))
      if (iperm(isymmo(indxj),isymmo(indxl)) .ne. isym) 
     +        go to 60
      ipos = ipo(indxi,indxk)
      jpos = ipo(indxj,indxl)
      call inbuck (iad,iblo4,irun,isym,ipos,jpos,2,q)
      call inbuck (iad,iblo4,irun,isym,jpos,ipos,2,q)

c   cases (ab|ji), (ba|ij)
c          -   -     - -

      isym = iperm(isymmo(indxi),isymmo(indxl))
      if (iperm(isymmo(indxj),isymmo(indxk)) .ne. isym) 
     +        go to 60
      ipos = ipo(indxi,indxl)
      jpos = ipo(indxj,indxk)
      call inbuck (iad,iblo4,irun,isym,ipos,jpos,2,q)
      call inbuck (iad,iblo4,irun,isym,jpos,ipos,2,q)
60    continue

      if (jj .ne. 0) go to 20
100   continue

c   finally write partially filled blocks onto external file

      do 120 ibox=1,ntot
      j = ifill(ibox)
      if (j .eq. 0) go to 120
      mlower = (ibox-1) * 511
      if (btest(j,0)) then
         kpos = mlower + 340 + j/2 + 1
_IF(cray,t3d)
         q(kpos) = shiftl(mask,32).and.q(kpos)
_ELSE
         pack = q(kpos)
         itemp(3) = 0
         itemp(4) = 0
         q(kpos) = pack
_ENDIF
      end if
      q(mlower+511) = dble(iad(ibox))
      call wrt3 (q(mlower+1),511,iblo4,ned4)
      iad(ibox) = iblo4*1000 + j
      iblo4 = iblo4 + 1
120   continue

      if (mtot .eq. 0) call caserr ('no integrals available')
      ratio = (dble(mtorpa)/dble(mtot)) * 100.0d0
      write (iwr,2010) mtot,nconta,mtot-nconta,mtorpa,ratio
2010  format ('total no of non-zero integrals: ',i10/
     +        'symm. contaminated integrals:   ',i10/
     +        'no of decontaminated integrals: ',i10/
     +        'total no of rpa integrals:      ',i10,1x,
     +        '(',f4.1,' % of all non-zero integrals)')

c   reordering step 2

      iblo7 = 1
      do 700 isym=1,nirr
      if (icount(isym).eq.0 .and. nfreq.eq.0) go to 700
      jstart(isym) = iblo7
      ipair = 0
      np = npair(isym)
      do 650 ibox=1,nbox(isym)
      if (ibox .eq. nbox(isym)) np = icount(isym) - ipair
      npp = np * icount(isym)
      call vclr(q,1,intspa)
      ii = iad(nbox0(isym)-1+ibox)
600   iblo4 = ii / 1000
      if (iblo4 .eq. 0) go to 620
      m = ii - 1000 * iblo4
      call rdedx (gin,511,iblo4,ned4)
      call upak4s (680,gin(num2ep+1),ibuf,jbuf)
      do 610 i=1,m
      ipos = ibuf(i) / 10
      itype = ibuf(i) - 10 * ipos
610   q((itype-1)*npp+(ipos-ipair-1)*icount(isym)+jbuf(i)) = gin(i)
      ii = idint(gin(511))
      go to 600
620   npp2 = npp + npp

c   construct A and B

      do 630 i=1,npp
      q(npp+i) = q(i) + q(i) - q(npp+i)
630   q(npp2+i) = q(npp2+i) - q(i) - q(i)

      do 640 ii=1,na
      do 640 ia=nc1,num
      if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
         k = ipo(ia,ii)
         if (k.gt.ipair .and. k.le.ipair+np) then
            kk = npp + (k-ipair-1)*icount(isym) + k
            q(kk) = q(kk) + evalue(ia) - evalue(ii)

c   prepare "forest"

            if (odebug(20)) then
               innx = iofsym(isym) + k
               inny = innx + nvc
               innz = inny + nvc
               xyz = dip(innx)*dip(innx) + dip(inny)*dip(inny) + 
     +               dip(innz)*dip(innz)
               write (ioneel) ii,ia,q(kk),xyz
            end if

         end if
      end if
640   continue

c   construct jaroslav''s small rpa matrix

      if (iorbit .gt. 0 .and. odebug(21) .and. isym.eq.jasym) then
c###################### construct the matrix ###########################
         do m=1,ndim
            ii=iipt(m)
            ia=iapt(m)
            if (iperm(isymmo(ii),isymmo(ia)) .eq. isym) then
               k=ipo(ii,ia)
               if (k.gt.ipair .and. k.le.ipair+np) then
                  do n=1,ndim
             ij=iipt(n)
             ib=iapt(n)
                     if (iperm(isymmo(ij),isymmo(ib)) .eq. isym) then
                        l=ipo(ij,ib)
                        in=(k-ipair-1)*icount(isym)+l
                        temp=q(npp+in)
                        a(m,n)=temp
                        a(m+ndim,n+ndim)=temp
                        temp=q(npp2+in)
                        a(m,n+ndim)=temp
                        a(m+ndim,n)=temp
                     endif
                  enddo
               endif
            endif
         enddo
c#################### print the matrix A #############################
         write(iwr,1050)isym
1050     format(/'small rpa matrix, symmetry ',i1,' -- part A'/)
         write(iwr,1051)(isymao(iipt(m)),symb(nirs+isymmo(iipt(m))),
     +                  m=1,ndim)
1051     format(17x,100(2x,i2,a3,2x))
         write(iwr,1052)
1052     format(17x,100(3x,'-->',3x))
         write(iwr,1053)(isymao(iapt(m)),symb(nirs+isymmo(iapt(m))),
     +                  m=1,ndim)
1053     format(17x,100(2x,i2,a3,2x))
         write(iwr,1054)
1054     format(/)
         do m=1,ndim
            write(iwr,1060) isymao(iipt(m)),symb(nirs+isymmo(iipt(m))),
     +       isymao(iapt(m)),symb(nirs+isymmo(iapt(m))),
     +       (a(m,n),n=1,ndim)
1060        format(i2,a3,' --> ',i2,a3,2x,100f9.4)
         enddo
c################### print the matrix B #############################
        write(iwr,1061)isym
1061     format(/'small rpa matrix, symmetry ',i1,' -- part B'/)
         write(iwr,1051)(isymao(iipt(m)),symb(nirs+isymmo(iipt(m))),
     +                  m=1,ndim)
         write(iwr,1052)
         write(iwr,1053)(isymao(iapt(m)),symb(nirs+isymmo(iapt(m))),
     +                  m=1,ndim)
         write(iwr,1054)
        do m=1,ndim
           write(iwr,1060) isymao(iipt(m)),symb(nirs+isymmo(iipt(m))),
     +     isymao(iapt(m)),symb(nirs+isymmo(iapt(m))),
     +     (a(m+ndim,n),n=1,ndim)
         enddo
c#################### diagonalize the matrix A ########################
         write (iwr,1069)
1069     format(/'eigenvectors of the matrix A:'/
     +           '-----------------------------'/)
        call tred2(200,ndim,a,eta,h,x)
        call tql2(200,ndim,eta,h,x,ierr)
         do n=1,ndim
           write(iwr,1070)n,eta(n)*27.21165d0
1070        format(/'eigenvalue ',i2,' = ',f9.4/30('-')/'eigenvector:')
            do m=1,ndim
              write(iwr,1071) x(m,n),isymao(iipt(m)),
     +         symb(nirs+isymmo(iipt(m))),
     +         isymao(iapt(m)),
     +         symb(nirs+isymmo(iapt(m)))
1071           format(f9.4,4x,i2,a3,' --> ',i2,a3)
            enddo
        enddo
c                                       |  A  B |
c############### diagonalize the matrix |       | ####################
c                                       | -B -A |
         write (iwr,1072)
1072     format(/'                           |  A  B |'
     +          /'eigenvectors of the matrix |       |:'
     +          /'                           | -B -A |'
     +          /'-------------------------------------'/)
         call vclr(s(1,1),1,200*200)
         do i=1,ndim
            s(i,i)=1.0d0
            s(i+ndim,i+ndim)=-1.0d0
         enddo
         call qzhes (200,2*ndim,a,s,.true.,x)
         call qzit (200,2*ndim,a,s,0.0d0,.true.,x,ierr)
         call qzval (200,2*ndim,a,s,eta,etaim,sigma,.true.,x)
         if (ierr .ne. 0) call caserr ('ierr .ne. 0 in rpsort')
         call qzvec (200,2*ndim,a,s,eta,etaim,sigma,x)
         nrealp = 0
         do i=1,2*ndim
            t = eta(i)
            u = etaim(i)
            v = sigma(i)
            if (1.0d0+u .eq. 1.0d0 .and. 1.0d0+v .ne. 1.0d0) then
               t = t / v
               if (t .ge. 0.0d0) then
                  nrealp = nrealp + 1
                  eta(nrealp) = t
                  call dcopy(2*ndim,x(1,i),1,x(1,nrealp),1)
               endif
            endif
         enddo
         write (iwr,*) 'nrealp = ',nrealp
         if(nrealp.ne.ndim) call caserr ('nrealp .ne. ndim')
         do i=1,nrealp
            k = i
            if (nrealp .gt. 1) then
               emin = eta(i)
               do j=i+1,nrealp
                  if (eta(j) .lt. emin) then
                     k = j
                     emin = eta(j)
                  end if
               enddo
            end if
            t = eta(i)
            eta(i) = eta(k)
            eta(k) = t
            call dswap(2*ndim,x(1,i),1,x(1,k),1)
         enddo
         do n=1,ndim
            write(iwr,1070)n,eta(n)*27.21165d0
            do m=1,ndim
               write(iwr,1073) x(m,n),x(m+ndim,n),
     +      isymao(iipt(m)),symb(nirs+isymmo(iipt(m))),
     +      isymao(iapt(m)),symb(nirs+isymmo(iapt(m)))
1073           format(2f9.4,4x,i2,a3,' --> ',i2,a3)
            enddo
         enddo
      endif

c   write the parts of A and B just constructed on external file

      call wrt3 (q(npp+1),npp2,iblo7,num8)
      ipair = ipair + np
650   iblo7 = iblo7 + lensec(npp2)
700   continue

      write (iwr,200) cpulft(1)-time
200   format ('collection of rpa integrals:',f14.4,' seconds'//122('='))
      if (odebug(20) .or. iorbit.gt.0) stop

      return
      end
      subroutine symort (n,x,y,sytol,oldep)

c---------------------------------------------------------------------
c   Performs the symmetric orthonormalization of the vectors 
c   B=(X,Y) and B#=(Y,X), where X and Y are vectors of dimension N
c   (c) Carsten Fuchs 1990
c---------------------------------------------------------------------

      implicit REAL (a-h,p-z), integer (i-n), logical (o)
      dimension x(*),y(*)
      oldep = .false.
      s = ddot(n,x,1,x,1) + 
     +    ddot(n,y,1,y,1)
      t = ddot(n,x,1,y,1)
      t = t + t
      if (s-dabs(t) .lt. sytol) then
         oldep = .true.
         return
      end if
      alpha = dsqrt( s + t )
      beta  = dsqrt( s - t )
      do 10 i=1,n
      u = x(i)
      v = y(i)
      x(i) = u + v
      y(i) = u - v
      u = x(i) / alpha
      v = y(i) / beta
      x(i) = 0.5d0 * ( u + v )
      y(i) = 0.5d0 * ( u - v )
10    continue
      return
      end
      subroutine tdalgo (diag,b,ared,eta,xred,h1,v1,resarr,scr)

c---------------------------------------------------------------------
c   Computes the lowest NEIG roots of the TDA matrix using the
c   Davidson algorithm
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit REAL (a-h,p-z), integer (i-n), logical (o)
INCLUDE(common/sizes)
      character *1 xn 
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      common/debug /odebug(100),oanal,ouse,oreduc,okeep,onew,oskip,
     +              ototal,oall
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      parameter (nvcmax=9953)
      common/dipblk/iconv(nvcmax),index(nvcmax)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
      common/iter  /ioff,niter,tim(4),isize2,isize3,ntot0,ierr
INCLUDE(common/rpacom)
INCLUDE(common/rpfile)
      common/symmet/isym,isize,neig,ipadsymmet(77)
      dimension diag(isize),h1(isize),v1(isize,nstart)
      dimension b(isize,maxr),ared(maxr,maxr)
      dimension eta(maxr),xred(maxr,maxr),resarr(*),scr(*)
      data xn /'n'/
      ntot1 = 1
      ntot = nstart
10    niter = niter+1
      if (niter .gt. maxit) go to 1000

c   calculate the vectors a * b(i) and fill up matrix ared

      do 30 i=ntot1,ntot
      time = cpulft(1)
      call rpmult (.false.,jstart(isym),npair(isym),isize,num8,
     +   b(1,i),dummy,h1,dummy,scr)
      tim(2) = tim(2)+cpulft(1)-time
      do 30 j=1,i
      ared(j,i) = ddot(isize,b(1,j),1,h1,1)
30    ared(i,j) = ared(j,i)

c   solve reduced eigenvalue problem

      time = cpulft(1)
      call tred2(maxr,ntot,ared,eta,h1,xred)
      call tql2(maxr,ntot,eta,h1,xred,ierr)
      if (ierr .ne. 0) go to 400
      tim(3) = tim(3) + cpulft(1) - time

c   collect lowest positive eigenvalues

      do 40 i=1,ntot
40    if (eta(i) .ge. 0.0d0) go to 50
      ioff = ntot
      go to 60
50    ioff = i - 1
60    nrealp = ntot - ioff
      if (nrealp .eq. 0) go to 500
      if (nrealp .lt. nstart) write (iwr,70) niter,nrealp
70    format ('warning: in iteration step ',i2,' there are only '
     +   ,i2,' positive real eigenvalues.'/'trying to continue ...')
      nev2 = min(nrealp,nstart)
      nev3 = nev2 - nevlo1
      call outeig (itable,iwr,eta(ioff+nevlo1+1),neig,niter,ntot,
     +             nevlo1,1)
c
c   new trial vectors to be expected

      ntot0 = ntot
      ntot1 = ntot + 1
      iev = 1
      do 140 j=ioff+nevlo(isym),ioff+nev2
      if (iconv(iev) .ne. 0) then
         resarr(iev) = -1.0d0
         go to 140
      end if

c   compute v(iev) for non-converged eigenvalue no iev

      call mxmaa (b,1,isize, xred(1,j),1,ntot0, v1(1,j),1,isize,
     +            isize,ntot0,1)

c   compute residue vector

      time = cpulft(1)
      call rpmult (.false.,jstart(isym),npair(isym),isize,num8,
     +  v1(1,j),dummy,h1,dummy,scr)
      tim(2) = tim(2) + cpulft(1) - time
      call daxmyz (isize,eta(j),v1(1,j),h1,h1)

c   compute euclidean norm of residue vector

      t = dnrm2(isize,h1,1)
      resarr(iev) = t

c   threshold ?

      if (t .le. eps) iconv(iev) = 1
      if (iconv(iev).eq.1 .or. ntot.ge.maxr .or. ntot.ge.isize) 
     +        go to 140
      ntot = ntot + 1

c   update residue vector 

      do 115 i = 1 , isize
         b(i,ntot) = eta(j) - diag(i)
115   continue
      call vdivz (h1,1,b(1,ntot),1,1.0d0,b(1,ntot),1,isize)
      call dcopy(isize,b(1,ntot),1,h1,1)

c   gram-schmidt orthogonalisation

      do 120 k=1,ntot-1
120   call orth (isize,b(1,ntot),b(1,k))
      t = dnrm2(isize,b(1,ntot),1)
      if (1.0d0+t .eq. 1.0d0) then
         write (iwr,130) iev,niter
130      format ('discard new trial vector no. ',i2,' in iteration',
     +           ' step ',i2/'because of linear dependence from',
     +           ' the previous trial vectors.')
         ntot = ntot - 1
      else
         call dscal(isize,1.0d0/t,b(1,ntot),1)
      end if
140   iev = iev + 1 
      call outres (iwr,inorms,resarr,nevlo1,neig,niter,1)
      nconv = isum(neig,iconv,1)
      if (nconv .eq. neig) go to 1000
c
c   if at least one new trial vector has been added, start 
c   next iteration step
c
      if (ntot .ge. ntot1) go to 10
      if (ntot .ge. maxr) then
      call dgemm(xn,xn,isize,nevlo1,maxr
     +          ,1.0d0,b,isize,xred,maxr,0.0d0,v1,isize)
        return
      end if
c
c   error messages
c
      write (iwr,160) niter
160   format ('no new trial vector has been added in iteration step ',
     +         i2, ':'/
     +   'unable to finish iterative calculation of tda eigenvalues.')
      ierr = 1
      go to 1000

400   write (iwr,410) niter,ierr
410   format (/'convergence problems with tql2g in iteration step',i3,
     +         ': ierr = ',i3/)
      go to 1000

500   write (iwr,510) niter
510   format (/
     +      'in iteration step ',i2,' there are no real eigenvalues.'/)
      ierr = 1

1000  niter = min(niter,maxit)

      return
      end
      subroutine tmfile (q1,q2)

c---------------------------------------------------------------------
c   Writes the basis set information and the MO coefficients to
c   fortran file itmfil which is used for plotting of transition 
c   density matrices. Later the RPA excited state "wavefunctions"
c   are also written to itmfil.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit REAL (a-h,p-w), integer (i-n), logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxcen3 = maxat * 3, maxterms = 10)
INCLUDE(common/blocko)
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/machin)
INCLUDE(common/nshel)
INCLUDE(common/rpacom)
INCLUDE(common/rpfile)
      common/symmet/isym,isize,neig,iperm(8,8),maxspa,intspa,ineed0,
     + ipadsymmet(10)
      common/tmdata/icen(maxorb),ixexp(maxorb),iyexp(maxorb),
     +              izexp(maxorb),ic(maxorb),
     +              ccoef(maxorb,maxterms),dzeta(maxorb,maxterms)
INCLUDE(common/tran)
INCLUDE(common/atmol3)
      dimension q1(*),q2(*),jshell(5),iexp(3,20),intyp(mxshel)
      data ntype / 5 /
      data jshell / 1 , 3 , 6 , 10 , 4 /
      data iexp / 
     + 0,0,0, 1,0,0, 0,1,0, 0,0,1, 2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1,
     + 0,1,1, 3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 1,2,0, 0,2,1, 1,0,2,
     + 0,1,2, 1,1,1 /
c
      if(opg_root()) then
       rewind itmfil
       write (itmfil) nev(isym),num,num,0,2*na,nat
       write (itmfil) (czan(iat),iat=1,nat)
       write (itmfil) (c(1, iat),iat=1,nat)
       write (itmfil) (c(2, iat),iat=1,nat)
       write (itmfil) (c(3, iat),iat=1,nat)
       write (itmfil) (j,j=1,num)
      endif
c
      do 110 loop = 1 , nshell
         ii = kmax(loop) - kmin(loop) + 1
         do 115 i = 1 , ntype
            if (ii .eq. jshell(i)) then
               intyp(loop) = i
            end if
115      continue
110   continue
c
      n = 0
      do 1000 ii = 1,nshell
         iat = katom(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         do 120 i = mini , maxi
            n = n + 1
            icen(n) = iat
            ixexp(n) = iexp(1,i)
            iyexp(n) = iexp(2,i)
            izexp(n) = iexp(3,i)
            ic(n) = kng(ii)
            i1 = kstart(ii)
            i2 = i1 + kng(ii) - 1
            ityp = intyp(ii)
            k = 0
            do 1520 ig = i1 , i2
               k = k + 1
               if (ityp .eq. 1) then
                  c1 = cs(ig)
               else if (ityp .eq. 2) then
                  c1 = cp(ig)
               else if (ityp .eq. 3) then
                  c1 = cd(ig)
               else if (ityp .eq. 4) then
                  c1 = cf(ig)
               else if (ityp .eq. 5) then
                  c1 = cs(ig)
c                 c3 = cp(ig)
               end if
               ccoef(n,k) = c1
               dzeta(n,k) = ex(ig)
1520        continue
120      continue
1000  continue

c   write centers / x,y,z exponents / number of primitives to tm file
      if(opg_root()) then
       write (itmfil) (icen(j),j=1,num)
       write (itmfil) (ixexp(j),j=1,num)
       write (itmfil) (iyexp(j),j=1,num)
       write (itmfil) (izexp(j),j=1,num)
       write (itmfil) (ic(j),j=1,num)

c   write contraction coefficients and exponents to tm file

       do j = 1 , num
         do k = 1 , ic(j)
           write (itmfil) ccoef(j,k),dzeta(j,k)
         enddo
       enddo
      endif
c
c   read in scf mos
c
      ibase = 0
      call secget (mouta,3,iblko)
      iblkq = iblko + lensec(mach(8)) + lensec(mach(9)) + 1
      call rdedx (q1,nbasq,iblkq,idaf)
      otran = .false.
      call tdown (q2,ilifq,q1,ilifq,num)
c
c   write mos to tm file
c
      if(opg_root()) then
       ibase = 0
       do 8004 i = 1 , num
          write (itmfil) (q2(ibase+j),j=1,num)
          ibase = ibase + num
8004   continue
      endif
c
      return
      end
      subroutine tmfil2 (q1,q2)

c---------------------------------------------------------------------
c   Writes the basis set information and the MO coefficients to
c   fortran file itmfil which is used for plotting of transition 
c   density matrices. Later the RPA excited state "wavefunctions"
c   are also written to itmfil.
c   (c) Carsten Fuchs 1991-1993
c---------------------------------------------------------------------

      implicit REAL (a-h,p-w), integer (i-n), logical (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      parameter (mxcen3 = maxat * 3, maxterms = 10)
INCLUDE(common/blocko)
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/machin)
INCLUDE(common/nshel)
INCLUDE(common/rpacom)
INCLUDE(common/rpfile)
      common/symmet/isym,isize,neig,iperm(8,8),maxspa,intspa,ineed0,
     + ipadsymmet(10)
      common/tmdata/icen(maxorb),ixexp(maxorb),iyexp(maxorb),
     +              izexp(maxorb),ic(maxorb),
     +              ccoef(maxorb,maxterms),dzeta(maxorb,maxterms)
INCLUDE(common/tran)
INCLUDE(common/atmol3)
      dimension q1(*),q2(*),jshell(5),iexp(3,20),intyp(mxshel)
      data ntype / 5 /
      data jshell / 1 , 3 , 6 , 10 , 4 /
      data iexp / 
     + 0,0,0, 1,0,0, 0,1,0, 0,0,1, 2,0,0, 0,2,0, 0,0,2, 1,1,0, 1,0,1,
     + 0,1,1, 3,0,0, 0,3,0, 0,0,3, 2,1,0, 2,0,1, 1,2,0, 0,2,1, 1,0,2,
     + 0,1,2, 1,1,1 /
c
      if(opg_root()) then
      rewind itmfil
      write (itmfil) nev(isym),num,num,0,2*na,nat
      write (itmfil) (czan(iat),iat=1,nat)
      write (itmfil) (c(1, iat),iat=1,nat)
      write (itmfil) (c(2, iat),iat=1,nat)
      write (itmfil) (c(3, iat),iat=1,nat)
      write (itmfil) (j,j=1,num)
c
      do 110 loop = 1 , nshell
         ii = kmax(loop) - kmin(loop) + 1
         do 115 i = 1 , ntype
            if (ii .eq. jshell(i)) then
               intyp(loop) = i
            end if
115      continue
110   continue
c
      n = 0
      do 1000 ii = 1,nshell
         iat = katom(ii)
         mini = kmin(ii)
         maxi = kmax(ii)
         do 120 i = mini , maxi
            n = n + 1
            icen(n) = iat
            ixexp(n) = iexp(1,i)
            iyexp(n) = iexp(2,i)
            izexp(n) = iexp(3,i)
            ic(n) = kng(ii)
            i1 = kstart(ii)
            i2 = i1 + kng(ii) - 1
            ityp = intyp(ii)
            k = 0
            do 1520 ig = i1 , i2
               k = k + 1
               if (ityp .eq. 1) then
                  c1 = cs(ig)
               else if (ityp .eq. 2) then
                  c1 = cp(ig)
               else if (ityp .eq. 3) then
                  c1 = cd(ig)
               else if (ityp .eq. 4) then
                  c1 = cf(ig)
               else if (ityp .eq. 5) then
                  c1 = cs(ig)
c                 c3 = cp(ig)
               end if
               ccoef(n,k) = c1
               dzeta(n,k) = ex(ig)
1520        continue
120      continue
1000  continue

c   write centers / x,y,z exponents / number of primitives to tm file

      write (itmfil) (icen(j),j=1,num)
      write (itmfil) (ixexp(j),j=1,num)
      write (itmfil) (iyexp(j),j=1,num)
      write (itmfil) (izexp(j),j=1,num)
      write (itmfil) (ic(j),j=1,num)

c   write contraction coefficients and exponents to tm file

      do 8000 j = 1 , num
         do 8001 k = 1 , ic(j)
            write (itmfil) ccoef(j,k),dzeta(j,k)
8001     continue
8000  continue
c
      endif
c   read in scf mos
c
      ibase = 0
      call secget (mouta,3,iblko)
      iblkq = iblko + lensec(mach(8)) + lensec(mach(9)) + 1
      call rdedx (q1,nbasq,iblkq,idaf)
      otran = .false.
      call tdown (q2,ilifq,q1,ilifq,num)
      if(opg_root()) then
c
c   write mos to tm file
c
       ibase = 0
       do 8004 i = 1 , num
          write (itmfil) (q2(ibase+j),j=1,num)
          ibase = ibase + num
8004   continue
      endif
c
      return
      end
_IFN(cray,fps,cyber205,ksr,i8)
      subroutine upak4s (nword,iii,i,j)
c---------------------------------------------------------------------
c   Unpacking routine.
c   (c) Carsten Fuchs 1991
c---------------------------------------------------------------------
      integer*2 iii(*),i(2,*),j(2,*)
      k = 0
      do 10 loop = 1 , nword , 4
         k = k + 1
_IF(littleendian)
         i(1,k) = iii(loop)
         j(1,k) = iii(loop+1)
         k = k + 1
         i(1,k) = iii(loop+2)
         j(1,k) = iii(loop+3)
_ELSE
         i(2,k) = iii(loop)
         j(2,k) = iii(loop+1)
         k = k + 1
         i(2,k) = iii(loop+2)
         j(2,k) = iii(loop+3)
_ENDIF
10    continue
      return
      end
_ELSE
      subroutine upak4s (nword,iii,i,j)
c---------------------------------------------------------------------
c   Unpacking routine.
c   (c) Carsten Fuchs 1991
c---------------------------------------------------------------------
      integer iii(*),i(*),j(*),ibuf(680)
      k = 0
      call unpack(iii,16,ibuf,680)
      do 10 loop = 1 , nword , 4
         k = k + 1
         i(k) = ibuf(loop)
         j(k) = ibuf(loop+1)
         k = k + 1
         i(k) = ibuf(loop+2)
         j(k) = ibuf(loop+3)
10    continue
      return
      end
_ENDIF
      subroutine vabs(a,ia,c,ic,n)

c*****  vabs    vector absolute value             math advantage rel 2.0
c
c    ** copyright 1984-1985 quantitative technology corporation **
c
c  call format
c
c       call vabs (a,ia,c,ic,n)
c
c       where,
c
c       a       real input vector.
c
c       ia      integer input stride for vector a.
c
c       c       real output vector.
c
c       ic      integer input stride for vector c.
c
c       n       integer input element count.
c
c
c  description
c
c       this routine computes the absolute value of the elements
c       of the vector a and stores the results in c.
c
c            c(i) = abs(a(i))    for i=1,n
c
c
c  example
c
c       call vabs (a,1,c,1,5)
c
c       input operands:
c
c       a = -1.00
c            2.00
c            3.00
c            4.00
c           -5.00
c
c       output operands:
c
c       c = 1.00
c           2.00
c           3.00
c           4.00
c           5.00
c
c  history
c         1) oct 84     d. cooper       original.
c
c
      integer ia,ic,n,ii,kk,m
      REAL c(*),a(*)
      if (n.le.0) go to 12
      ii = 1
      kk = 1
      do 10 m=1,n
        c(kk) = dabs(a(ii))
        ii = ii + ia
        kk = kk + ic
10    continue
12    return
      end
      subroutine vdivz(a,ia,b,ib,c,d,id,n)

c*****  vdivz   vector divide with zero test math advantage rel 2.0
c
c    ** copyright 1986 quantitative technology corporation **
c
c  call format
c
c       call vdivz (a,ia,b,ib,c,d,id,n)
c
c       where,
c
c       a       real input vector.
c
c       ia      integer input stride for vector a.
c
c       b       real input vector.
c
c       ib      integer input stride for vector b.
c
c       c       real input scalar.
c
c       d       real output vector.
c
c       id      integer input stride for vector d.
c
c       n       integer input element count.
c
c
c  description
c
c       this routine divides the elements of a vector a by
c       the corresponding non-zero elements of a vector b 
c       and stores the results in vector d.  elements of
c       vector d corresponding to zero elements of vector
c       b are set equal to c.
c
c         d(i) = a(i) / b(i)  if b(i) <> 0.0
c         d(i) = c            if b(i) = 0.0  for i=1,n
c
c
c  example
c
c       call vdivz (a,1,b,1,c,d,1,5)
c
c       input operands:
c
c       a = 1.0     b = 0.50     c = 3.33
c           2.0         0.40
c           3.0         0.30
c           4.0         0.00
c           5.0         0.10
c
c       output operands:
c
c       d = 2.00
c           5.00
c          10.00
c           3.33
c          50.00
c
c  history
c         1) jul 86     d. benua      original.
c
c
      integer ia,ib,id,n,ii,jj,kk,m
      REAL a(*),b(*),c,d(*)
      if (n.le.0) go to 22
      ii = 1
      jj = 1
      kk = 1
      do 20 m=1,n
        if (b(jj) .eq. 0.0d0) goto 10
        d(kk) = a(ii) / b(jj)
        goto 12
10      d(kk) = c
12      ii = ii+ia
        jj = jj+ib
        kk = kk+id
20    continue
22    return
      end
      subroutine ver_rpa(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/rpa.m,v $
     +     "/
      data revision /"$Revision: 6277 $"/
      data date /"$Date: 2013-02-10 13:14:47 +0100 (Sun, 10 Feb 2013) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
