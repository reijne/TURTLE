_IF(ga)
      subroutine diis_ga(h0,h1,h2,iblkqa,diff0,idomo,ivmo)
c
c --- subroutine sets up and solves the diis equations
c             direct inversion of iterated subspace
c --- p.pulay chem.phys.lett. 73 (1980) 393
c
c     this version uses ga storage
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension h0(*),h1(*),h2(*)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/scra7)
      common/diisd/st(210),ct(20),r(19),derror,scale(20),
     +             iposit(20),nstore,mp,ondiis
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/harmon)
INCLUDE(common/gadiis)
      common/scftim/tdiag(4),tdiis,tmult(5)
c
c the error matrix is the lower triangle of off diag
c fock matrix elements in the mo basis. this is then transformed back
c to the ao basis and written to ga in square form
c
c arguments:
c
c   h0       input fock, overwritten by result
c   h1       scratch (triangle)
c   h2       scratch (square)
c   iblkqa   address of vectors
c   diff0    diis tester
c
c storage in /diisd/:
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored on ed7, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c        is to be used in the scf program
c accdi1  diis comes in only when diff is less than this
c
c
INCLUDE(common/prints)
c

      if (supressed) call caserr('diis - no gas')

      dumtim=dclock()
      derror=0.0d0
      num0 = newbas0
      l3 = num*num
      idomo = 0
      ivomo = 0

      nmin=3
      nmax=maxdiis
c
c --- should we begin storing diis info ???
c
      al2=nx
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax)nmin=nmax-1
      if(nmin.lt.1)then
c
c too few degrees of freedom to run diis
c
         ondiis = .false.
         return
      endif
c
c ----sort out the vectors, store in ga
c
      call rdedx(h2,l3,iblkqa,idaf)
      call tdown(h2,ilifq,h2,ilifq,num0)
      call load_ga_from_square(ih_vec,h2,num)
c
c store input fock matrix to ga location (ipos,1) and
c (ipos,2) (to be converted to error vector in parallel).
c
      ipos=mod(nstore,nmax)+1
      call load_ga_from_triangle(ihandle(ipos,1),h0,num)
      call load_ga_from_triangle(ihandle(ipos,2),h0,num)
c
c ----- calculate the error vector, 
c       transform to mo basis
c
      call mult2_ga(ihandle(ipos,2), ih_vec, ih_scr, num)
c
c       zero oo and vv blocks, and find absmax of ov
c
      call pg_dscal_patch(ihandle(ipos,2),1,na,1,na,0.0d0)
      call pg_dscal_patch(ihandle(ipos,2),na+1,num0,na+1,num0,0.0d0)
      diff0 = gax_absmax_patch(ihandle(ipos,2),1,na,na+1,num0,
     +                         idomo,ivmo)

      if(diff0.gt.accdi1.or.diff0.eq.0.0d0) then
c
c components of error vector > accdi1, turn off diis
c
         nstore=0
         itot=0
         mp=0
         ondiis=.false.
c        if(opg_root())write(6,*)'diis reject on diff0', diff0
         goto 9999
      endif
c
c ------ transform back to the ao basis
c
      call mult2t_ga(ihandle(ipos,2),ih_vec, ih_scr, num)
      call mult2_ga(ihandle(ipos,2),ih_ov, ih_scr, num)

      nstore=nstore+1
      ipos=mod(nstore-1,nmax)+1

      ipos1=ipos-1
      mp=iky(ipos)
      if(ipos1.ge.1) then
         do 50 i=1,ipos1
            mp=mp+1
            st(mp)=pg_ddot(ihandle(i,2),ihandle(ipos,2))
 50      continue
      endif
c
c store diagonal term
      mp=mp+1
      st(mp)=pg_ddot(ihandle(ipos,2),ihandle(ipos,2))

      
      itot=nstore
      if(nstore.gt.nmax)itot=nmax

      ipos1=ipos+1
      if(ipos1.le.itot) then
         do 110 i=ipos1,itot
            mp=mp+i-1
            st(mp)=pg_ddot(ihandle(i,2),ihandle(ipos,2))
 110     continue
      endif

      if(nstore.le.nmin)then
c
c not enough triangles yet, restore the input vectors and 
c continue
c
         ondiis = .false.
         goto 9999

      endif
c
c --- now solve the diis equations
c
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
c
      call square(h0,st,ndim,ndim)
      do 125 i=1,ndim
 125  scale(i)=dsqrt(st(ikyp(i)))
c
c --- scale the matrix
c
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0(mp1)=h0(mp1)/cij
      sc=h0(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0(k1)=h0(k1)/sc
 123  h0(k2)=h0(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1
      if(oprint(47))write(iwr,126)(st(k),k=1,mp)
      call f04atf(h0,ndim,r,ndim,ct,h1,ndim,h2,h2(ndim+1),ifail)
      if(ifail.ne.0) then
c
c diis failure - exit, but retain the current point
c
         if(oprint(47))write(iwr,129)ifail
 129     format(//1x,'diis failure in f04atf .... ifail= ',i2)

c swap ga handles to bring current h,e to ipos=1
c this has not been checked
c
         itmp = ihandle(ipos,1)
         ihandle(ipos,1) = ihandle(1,1)
         ihandle(1,1) = itmp

         itmp = ihandle(ipos,2)
         ihandle(ipos,2) = ihandle(1,2)
         ihandle(1,2) = itmp

         nstore=1
         itot=1
         mp=iky(ipos)+ipos
         st(1)=st(mp)

         ondiis = .false.
         goto 9999

      endif

      do 114 i=1,ndim
 114     ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(oprint(47))write(iwr,128)(ct(k),k=1,ndim)
 126  format(/1x,'diis matrix'/(1x,10e11.4))
 128  format(/1x,'solution'/(1x,10e11.4))
c
c construct extrapolated fock
c
      call pg_zero(ih_scr)
      do 150 k=1,itot
         call pg_dadd(1.0d0,ih_scr,ct(k),ihandle(k,1),ih_scr)
 150  continue

c
c unpack square ga -> triangular result
c
      call load_triangle_from_ga(h0,ih_scr,num)

      if(.not.ondiis) kcount=0
      ondiis=.true.
c
c exit diis - if we failed to extrapolate, restore the original
c fock matrix
c
9999  continue

      if(.not.ondiis)
     &     call load_triangle_from_ga(h0,ihandle(ipos,1),num)

      tdiis=tdiis+(dclock()-dumtim)
      return
      end
c
      subroutine diisu_ga(h0a,h0b,h1,h2,iblkqa,iblkqb,
     +                    diffa,diffb,
     +                    idomoa,ivmoa,idomob,ivmob)
c
c --- subroutine sets up and solves the diis equations for uhf
c     direct inversion of iterated subspace
c --- p.pulay chem.phys.lett. 73 (1980) 393
c
c     this version uses ga storage
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      dimension h0a(*),h0b(*),h1(*),h2(*)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/scra7)
INCLUDE(common/harmon)
      common/diisd/st(210),ct(20),r(19),derror,scale(20),
     +             iposit(20),nstore,mp,ondiis
INCLUDE(common/scfopt)
INCLUDE(common/iofile)
INCLUDE(common/gadiis)
      common/scftim/tdiag(4),tdiis,tmult(5)
c
c the error matrix is the lower triangle of off diag
c fock matrix elements in the mo basis. this is then transformed back
c to the ao basis and written to ga in square form
c
c arguments:
c
c   h0a      input alpha fock, overwritten by result
c   h0b      input beta fock, overwritten by result
c   h1       scratch (triangle)
c   h2       scratch (square)
c   iblkqa   address of alpha vectors
c   iblkqb   address of beta vectors
c   diffa    diis tester for alpha MOs
c   diffb    diis tester for beta MOs
c
c storage in /diisd/:
c
c st     holds the lower triangle of the diis equations
c ct     holds the solutions
c r      holds the rhs
c nstore the number of arrays stored on ed7, initialise to
c        zero before call to scf routine
c mp     current position in the lower triangle
c        is to be used in the scf program
c accdi1  diis comes in only when diff is less than this
c
c
INCLUDE(common/prints)
c

      if (supressed) call caserr('diis - no gas')

      dumtim=dclock()
      derror=0.0d0
      num0 = newbas0
      l3=num*num
      idomoa = 0
      ivmoa = 0
      idomob = 0
      ivmob = 0

      nmin=3
      nmax=maxdiis
c
c --- should we begin storing diis info ???
c
      al2=nx
      ispace=dsqrt(al2)+0.01d0
      if(nmax.ge.ispace) nmax=ispace-1
      if(nmin.ge.nmax)nmin=nmax-1
      if(nmin.lt.1)then
c
c too few degrees of freedom to run diis
c
         ondiis = .false.
         return
      endif
c
c ----sort out the vectors, store in ga
c ----first the alpha vectors
c
      call rdedx(h2,l3,iblkqa,idaf)
      call tdown(h2,ilifq,h2,ilifq,num0)
      call load_ga_from_square(ih_vec,h2,num)
c
c store input alpha fock matrix to ga location (ipos,1) and
c (ipos,2) (to be converted to error vector in parallel).
c
      ipos=mod(nstore,nmax)+1
      call load_ga_from_triangle(ihandle(ipos,1),h0a,num)
      call load_ga_from_triangle(ihandle(ipos,2),h0a,num)
c
c ----- calculate the error vector, transform to mo basis
c
      call mult2_ga(ihandle(ipos,2), ih_vec, ih_scr, num)
c
c       zero oo and vv blocks, and find absmax of ov
c
      call pg_dscal_patch(ihandle(ipos,2),1,na,1,na,0.0d0)
      call pg_dscal_patch(ihandle(ipos,2),na+1,num0,na+1,num0,0.0d0)
      diffa = gax_absmax_patch(ihandle(ipos,2),1,na,na+1,num0,
     +                         idomoa,ivmoa)
c
c ----now the beta vectors
c
      call rdedx(h2,l3,iblkqb,idaf)
      call tdown(h2,ilifq,h2,ilifq,num0)
      call load_ga_from_square(ih_vecb,h2,num)
c
c store input beta fock matrix to ga location (ipos,3) and
c (ipos,4) (to be converted to error vector in parallel).
c
      call load_ga_from_triangle(ihandle(ipos,3),h0b,num)
      call load_ga_from_triangle(ihandle(ipos,4),h0b,num)
c
c ----- calculate the error vector, transform to mo basis
c
      call mult2_ga(ihandle(ipos,4), ih_vecb, ih_scr, num)
c
c       zero oo and vv blocks, and find absmax of ov
c
      call pg_dscal_patch(ihandle(ipos,4),1,nb,1,nb,0.0d0)
      call pg_dscal_patch(ihandle(ipos,4),nb+1,num0,nb+1,num0,0.0d0)

      diffb = gax_absmax_patch(ihandle(ipos,4),1,nb,nb+1,num0,
     +                         idomob,ivmob)
      diff0 = max(diffa,diffb)
      if(diff0.gt.accdi1.or.diff0.eq.0.0d0) then
c
c components of error vector > accdi1, turn off diis
c
         nstore=0
         itot=0
         mp=0
         ondiis=.false.
c        if(opg_root())write(6,*)'diis reject on diff', diff0
         goto 9999
      endif
c
c ------ transform back to the ao basis - alpha
c
      call mult2t_ga(ihandle(ipos,2),ih_vec, ih_scr, num)
      call mult2_ga(ihandle(ipos,2),ih_ov, ih_scr, num)

c
c ------ transform back to the ao basis - alpha
c
      call mult2t_ga(ihandle(ipos,4),ih_vecb, ih_scr, num)
      call mult2_ga(ihandle(ipos,4),ih_ov, ih_scr, num)

      nstore=nstore+1
      ipos=mod(nstore-1,nmax)+1

      ipos1=ipos-1
      mp=iky(ipos)
      if(ipos1.ge.1) then
         do 50 i=1,ipos1
            mp=mp+1
            st(mp)=pg_ddot(ihandle(i,2),ihandle(ipos,2)) +
     +             pg_ddot(ihandle(i,4),ihandle(ipos,4))
 50      continue
      endif
c
c store diagonal term
      mp=mp+1
      st(mp)=pg_ddot(ihandle(ipos,2),ihandle(ipos,2)) +
     +       pg_ddot(ihandle(ipos,4),ihandle(ipos,4))

      
      itot=nstore
      if(nstore.gt.nmax)itot=nmax

      ipos1=ipos+1
      if(ipos1.le.itot) then
         do 110 i=ipos1,itot
            mp=mp+i-1
            st(mp)=pg_ddot(ihandle(i,2),ihandle(ipos,2)) +
     +             pg_ddot(ihandle(i,4),ihandle(ipos,4))

 110     continue
      endif

      if(nstore.le.nmin)then
c
c not enough triangles yet, restore the input vectors and 
c continue
c
         ondiis = .false.
         goto 9999

      endif
c
c --- now solve the diis equations
c
      ndim=itot+1
      mp=iky(ndim)
      mp1=mp
      do 120 i=1,ndim
      mp1=mp1+1
      st(mp1)=-1.0d0
 120  r(i)=0.0d0
      st(mp1)=0.0d0
      r(ndim)=-1.0d0
c
      call square(h0a,st,ndim,ndim)
      do 125 i=1,ndim
 125  scale(i)=dsqrt(st(ikyp(i)))
c
c --- scale the matrix
c
      scale(ndim)=1.0d0
      mp1=0
      do 122 i=1,ndim
      ci=scale(i)
      do 122 j=1,ndim
      cij=ci*scale(j)
      mp1=mp1+1
 122  h0a(mp1)=h0a(mp1)/cij
      sc=h0a(ndim)
      sc= dabs(sc)
      scale(ndim)=sc
      do 123 i=1,ndim
      k1=(i-1)*ndim+ndim
      k2=(ndim-1)*ndim+i
      h0a(k1)=h0a(k1)/sc
 123  h0a(k2)=h0a(k2)/sc
      r(ndim)=r(ndim)/sc
      ifail=1
      if(oprint(47))write(iwr,126)(st(k),k=1,mp)
      call f04atf(h0a,ndim,r,ndim,ct,h1,ndim,h2,h2(ndim+1),ifail)
      if(ifail.ne.0) then
c
c diis failure - exit, but retain the current point
c
         if(oprint(47))write(iwr,129)ifail
 129     format(//1x,'diis failure in f04atf .... ifail= ',i2)

c swap ga handles to bring current h,e to ipos=1
c this has not been checked
c alpha
         itmp = ihandle(ipos,1)
         ihandle(ipos,1) = ihandle(1,1)
         ihandle(1,1) = itmp

         itmp = ihandle(ipos,2)
         ihandle(ipos,2) = ihandle(1,2)
         ihandle(1,2) = itmp
c beta
         itmp = ihandle(ipos,3)
         ihandle(ipos,3) = ihandle(1,3)
         ihandle(1,3) = itmp

         itmp = ihandle(ipos,4)
         ihandle(ipos,4) = ihandle(1,4)
         ihandle(1,4) = itmp

         nstore=1
         itot=1
         mp=iky(ipos)+ipos
         st(1)=st(mp)

         ondiis = .false.
         goto 9999

      endif

      do 114 i=1,ndim
 114     ct(i)=ct(i)/scale(i)
      derror=ct(ndim)
      if(oprint(47))write(iwr,128)(ct(k),k=1,ndim)
 126  format(/1x,'diis matrix'/(1x,10e11.4))
 128  format(/1x,'solution'/(1x,10e11.4))
c
c construct extrapolated alpha-fock
c
      call pg_zero(ih_scr)
      do  k=1,itot
         call pg_dadd(1.0d0,ih_scr,ct(k),ihandle(k,1),ih_scr)
      enddo
c
c unpack square ga -> triangular result
c
      call load_triangle_from_ga(h0a,ih_scr,num)
c
c construct extrapolated beta-fock
c
      call pg_zero(ih_scr)
      do k=1,itot
         call pg_dadd(1.0d0,ih_scr,ct(k),ihandle(k,3),ih_scr)
      enddo
c
c unpack square ga -> triangular result
c
      call load_triangle_from_ga(h0b,ih_scr,num)

      if(.not.ondiis) kcount=0
      ondiis=.true.
c
c exit diis - if we failed to extrapolate, restore the original
c fock matrix
c
9999  continue

      if(.not.ondiis) then
         call load_triangle_from_ga(h0a,ihandle(ipos,1),num)
         call load_triangle_from_ga(h0b,ihandle(ipos,3),num)
      endif

      tdiis=tdiis+(dclock()-dumtim)
      return
      end
      subroutine declare_diis_storage(ndim,beta)
      implicit none
      integer ndim, igmem_alloc
      logical beta

INCLUDE(common/gadiis)
INCLUDE(common/parcntl)
INCLUDE(common/iofile)
INCLUDE(common/gmempara)

      character*2 l1(maxdiis), l2(maxdiis)
      data l1 / 'f1','f2','f3','f4','f5','f6','f7','f8'/
      data l2 / 'e1','e2','e3','e4','e5','e6','e7','e8'/

      logical pg_create_inf

      integer i, ltri, ibuff1, ibuff2

      if (declared) return

      if(idpdiis .eq. 99999999 .and.
     &   idporth .eq. 99999999 .and.
     &   idpmult2 .eq. 99999999)then
c
         write(iwr,*)'diis ga allocation supressed'

         supressed = .true.
         return

      endif	

c
c maxdiis pairs of buffers for fock and error vectors
c
      do i = 1, maxdiis

         if (.not. pg_create_inf(0,ndim,ndim,l1(i),0,0,
     &        ihandle(i,1),'ga.m','declare_diis_storage',IGMEM_NORMAL ))
     &   then
            call pg_error('failed to create fock ga ',i)
         end if
         if (.not. pg_create_inf(0,ndim,ndim,l2(i),0,0,
     &        ihandle(i,2),'ga.m','declare_diis_storage',IGMEM_NORMAL ))
     &   then
            call pg_error('failed to create error ga ',i)
         end if

         if (beta) then

c maxdiis pairs of buffers for beta fock and error vectors

           if (.not. pg_create_inf(0,ndim,ndim,l1(i),0,0,
     &          ihandle(i,3),'ga.m','declare_diis_storage',
     &          IGMEM_NORMAL )) then
              call pg_error('failed to create beta fock ga ',i)
           end if
           if (.not. pg_create_inf(0,ndim,ndim,l2(i),0,0,
     &          ihandle(i,4),'ga.m','declare_diis_storage',
     &          IGMEM_NORMAL  )) then
              call pg_error('failed to create beta error ga ',i)
           end if

         end if

      enddo
c
c one square buffer for scratch sum/mult2
c
      if (.not. pg_create_inf(0,ndim,ndim,'sqbuf',0,0,
     &     ih_scr,'ga.m','declare_diis_storage', IGMEM_NORMAL)) then
         call pg_error('failed to create buffer ga ',i)
      endif
c
c square buffer for scratch sum/mult2
c
      if (.not. pg_create_inf(0,ndim,ndim,'sqbuf2',0,0,
     &     ih_scr2,'ga.m','declare_diis_storage', IGMEM_NORMAL)) then
         call pg_error('failed to create buffer ga ',i)
      endif
c
c one square buffer for vectors
c
      if (.not. pg_create_inf(0,ndim,ndim,'vec',0,0,
     &     ih_vec,'ga.m','declare_diis_storage', IGMEM_NORMAL)) then
         call pg_error('failed to create vec ga ',i)
      endif
c
c and beta vectors
c
      if (beta) then
        if (.not. pg_create_inf(0,ndim,ndim,'vec_b',0,0,
     &       ih_vecb,'ga.m','declare_diis_storage', IGMEM_NORMAL)) then
           call pg_error('failed to create beta vec ga ',i)
        endif
      endif
c
c buffer for overlap matrix  (store now)
c
      ltri = (ndim+1)*ndim/2

c
      if (.not. pg_create_inf(0,ndim,ndim,'overlap',0,0,
     &     ih_ov,'ga.m','declare_diis_storage', IGMEM_NORMAL)) then
         call pg_error('failed to create overlap ga ',i)
      endif
c
c workspace for transformed overlap (store now)
c
      if (.not. pg_create_inf(0,ndim,ndim,'overlap-trans',0,0,
     &     ih_ovt,'ga.m','declare_diis_storage', IGMEM_NORMAL)) then
         call pg_error('failed to create tr overlap ga ',i)
      endif


      declared = .true.

      end

      subroutine init_diis_storage(ndim,s,st)
      implicit none
      REAL s(*), st(*)
      integer ndim

_IF(charmm,chemshell)
*
* Ensure loading of block data (cray requirement)
*
      external block_diis_storage
_ENDIF

INCLUDE(common/gadiis)
      if (supressed) return
      if (.not. declared) call caserr('init_diis_storage error')
      call load_ga_from_triangle(ih_ov,s,ndim)
      call load_ga_from_triangle(ih_ovt,st,ndim)
      end	

      subroutine destroy_diis_storage (beta)
INCLUDE(common/gadiis)
      logical pg_destroy_inf,o,beta
      external pg_destroy_inf

      character*2 l1(maxdiis), l2(maxdiis)
      data l1 / 'f1','f2','f3','f4','f5','f6','f7','f8'/
      data l2 / 'e1','e2','e3','e4','e5','e6','e7','e8'/

      if (supressed) return
      o = .true.
      do i=1,maxdiis
         o = o .and. pg_destroy_inf(ihandle(i,1),l1(i),'ga.m',
     &                              'destroy_diis_storage')
         o = o .and. pg_destroy_inf(ihandle(i,2),l2(i),'ga.m',
     &                              'destroy_diis_storage')
       if (beta) then
          o = o .and. pg_destroy_inf(ihandle(i,3),l1(i),'ga.m',
     &                               'destroy_diis_storage')
          o = o .and. pg_destroy_inf(ihandle(i,4),l2(i),'ga.m',
     &                               'destroy_diis_storage')
       endif
      enddo
      o = o .and. pg_destroy_inf(ih_ov,'overlap','ga.m',
     &                           'destroy_diis_storage')
      o = o .and. pg_destroy_inf(ih_ovt,'overlap-trans','ga.m',
     &                           'destroy_diis_storage')
      o = o .and. pg_destroy_inf(ih_scr,'sqbuf','ga.m',
     &                           'destroy_diis_storage')
      o = o .and. pg_destroy_inf(ih_scr2,'sqbuf2','ga.m',
     &                           'destroy_diis_storage')
      o = o .and. pg_destroy_inf(ih_vec,'vec','ga.m',
     &                           'destroy_diis_storage')
      if (beta) o = o .and. pg_destroy_inf(ih_vecb,'vec_b','ga.m',
     &                                     'destroy_diis_storage')
      if(.not. o)call caserr('error destrying diis gas')
      declared = .false.
      end

      block data block_diis_storage
INCLUDE(common/gadiis)
      data declared /.false./
      data supressed /.false./
      end
c
c  load a square ga from input triangle
c
      subroutine load_ga_from_triangle(handle,f,ndim)
      implicit none
      integer handle
      REAL f(*)
      integer ndim
      integer ilo, ihi, jlo, jhi, nx, ny, ibuff
      integer ipg_nodeid, igmem_alloc_inf, i1, i2, i, j, ix, jx
INCLUDE(common/vcore)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/mapper)
      character *4  fnm
      character *21 snm
      data fnm/'ga.m'/
      data snm/'load_ga_from_triangle'/

      call pg_distribution(handle, ipg_nodeid(), ilo, ihi, jlo, jhi)

      if(ihi.le.0)return
      if(jhi.le.0)return

      nx = ihi - ilo + 1
      ny = jhi - jlo + 1
      
      ibuff = igmem_alloc_inf(nx*ny,fnm,snm,'buffer',IGMEM_DEBUG)

      do i = 1,nx
         do j = 1,ny
            ix = i + ilo - 1
            jx = j + jlo - 1
            i1 = min(ix,jx)
            i2 = max(ix,jx)
            Q(ibuff + nx*(j-1) + (i-1) ) = f(iky(i2) + i1)
         enddo
      enddo

      call pg_put(handle,ilo,ihi,jlo,jhi,
     &              Q(ibuff),nx)
      
      call gmem_free_inf(ibuff,fnm,snm,'buffer')

      return
      end

c
c  load a square ga from input square
c
      subroutine load_ga_from_square(handle,f,ndim)
      implicit none
      integer handle
      integer ndim
      REAL f(ndim,ndim)
      integer ilo, ihi, jlo, jhi, nx, ny, ibuff
      integer ipg_nodeid, igmem_alloc_inf, i1, i2, i, j, ix, jx
INCLUDE(common/vcore)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/mapper)
      character *4  fnm
      character *19 snm
      data fnm/'ga.m'/
      data snm/'load_ga_from_square'/

      call pg_distribution(handle, ipg_nodeid(), ilo, ihi, jlo, jhi)

      if(ihi.le.0)return
      if(jhi.le.0)return

      nx = ihi - ilo + 1
      ny = jhi - jlo + 1

      ibuff = igmem_alloc_inf(nx*ny,fnm,snm,'buffer',IGMEM_DEBUG)

      do i = 1,nx
         do j = 1,ny
            ix = i + ilo - 1
            jx = j + jlo - 1
            Q(ibuff + nx*(j-1) + (i-1) ) = f(ix,jx)
         enddo
      enddo

      call pg_put(handle,ilo,ihi,jlo,jhi,
     &              Q(ibuff),nx)
      
      call gmem_free_inf(ibuff,fnm,snm,'buffer')

      return
      end

c
c  load triangle from a square ga
c
      subroutine load_triangle_from_ga(f,handle,ndim)
      implicit none
      integer handle
      REAL f(*)
      integer ndim
      integer ilo, ihi, jlo, jhi, nx, ny, ibuff
      integer ipg_nodeid, igmem_alloc, i1, i2, i, j, ltri, ix, jx
INCLUDE(common/vcore)
INCLUDE(common/sizes)
INCLUDE(common/mapper)

      ltri = (ndim+1)*ndim/2

      call pg_distribution(handle, ipg_nodeid(), ilo, ihi, jlo, jhi)

      nx = ihi - ilo + 1
      ny = jhi - jlo + 1

      ibuff = igmem_alloc(nx*ny)

      call pg_get(handle,ilo,ihi,jlo,jhi,
     &              Q(ibuff),nx)
      
      call dcopy(ltri,0.0d00,0,f,1)

      do i = 1,nx
         do j = 1,ny
            ix = i + ilo - 1
            jx = j + jlo - 1
            if(ix.le.jx)then
               f(iky(jx) + ix) = Q(ibuff + nx*(j-1) + (i-1)) 
            endif
         enddo
      enddo

      call pg_dgop(1001, f, ltri, '+')

      call gmem_free(ibuff)
      return
      end
c
c  load square from a square ga
c
      subroutine load_square_from_ga(f,handle,ndim)
      implicit none
      integer handle
c     the pg_synch changes below are needed on the pwr4-based
c     IBM systems, but cause the PSC TCS Alphaserver system to
c     generate wrong results! Included now as comments until
c     properly resolved.
***   integer handle2
***   parameter (handle2=5566)
      integer ndim
      REAL f(ndim,ndim)
      integer ilo, ihi, jlo, jhi, nx, ny, ibuff
      integer ipg_nodeid, igmem_alloc, i1, i2, i, j, lsq, ix, jx
INCLUDE(common/vcore)
INCLUDE(common/sizes)
INCLUDE(common/mapper)

c
c so far ga_get seems faster than dgop
c
      call pg_get(handle,1,ndim,1,ndim,f,ndim)
***   call pg_synch(handle2)
c
c  synchronisation ??

cc      lsq = ndim*ndim
cc      call pg_distribution(handle, ipg_nodeid(), ilo, ihi, jlo, jhi)
cc
cc      nx = ihi - ilo + 1
cc      ny = jhi - jlo + 1
cc
cc      ibuff = igmem_alloc(lsq)
cc
cc      call pg_get(handle,ilo,ihi,jlo,jhi,
cc     &              Q(ibuff),nx)
cc      
cc      call dcopy(lsq,0.0d00,0,f,1)
cc
cc      do i = 1,nx
cc         do j = 1,ny
cc            ix = i + ilo - 1
cc            jx = j + jlo - 1
cc            f(ix,jx) = Q(ibuff + nx*(j-1) + (i-1)) 
cc         enddo
cc      enddo
cc
cc      call pg_dgop(1001, f, lsq, '+')
cc
cc      call gmem_free(ibuff)

      return
      end
c
c mult2_ga - version of q(dagger) h q where all matrices are
c            stored in gas
c            result is written back into h
c            scratch ga is required
c
      subroutine mult2_ga(h_handle, q_handle,
     &     scr_handle, ndim)
      implicit none
      integer h_handle, q_handle, scr_handle
      integer ndim
INCLUDE(common/timeperiods)
      character*1 xn,xt
      data xn,xt/'n','t'/

      call start_time_period(TP_GAMULT2)
      call pg_dgemm(xt,xn,ndim,ndim,ndim,1.0d0,
     &     q_handle, h_handle, 0.0d0, scr_handle)
      call pg_dgemm(xn,xn,ndim,ndim,ndim,1.0d0,
     &     scr_handle, q_handle, 0.0d0, h_handle)

      call end_time_period(TP_GAMULT2)
      end
c
c mult2t_ga - version of q h q(dagger) where all matrices are
c             stored in gas
c             result is written back into h
c             scratch ga is required
c
      subroutine mult2t_ga(h_handle, q_handle,
     &     scr_handle, ndim)
      implicit none
      integer h_handle, q_handle, scr_handle
      integer ndim
INCLUDE(common/timeperiods)
      character*1 xn,xt
      data xn,xt/'n','t'/

      call start_time_period(TP_GAMULT2)
      call pg_dgemm(xn,xn,ndim,ndim,ndim,1.0d0,
     &     q_handle, h_handle, 0.0d0, scr_handle)

      call pg_dgemm(xn,xt,ndim,ndim,ndim,1.0d0,
     &     scr_handle, q_handle, 0.0d0, h_handle)

      call end_time_period(TP_GAMULT2)
      end

      REAL function gax_absmax_patch(handle,ailo,aihi,ajlo,ajhi,
     +                               idomo,ivmo)
      implicit none
      integer handle,ailo,aihi,ajlo,ajhi
      integer idomo,ivmo
      integer ilo, ihi, jlo, jhi, ibuff, nx, ny, ix, jx, i, j
      integer ipg_nodeid, ipg_nnodes, igmem_alloc

INCLUDE(common/sizes)
INCLUDE(common/vcore)

      integer indexi(2,0:mxproc-1)
      REAL result, global_result
      integer me,nn,imax,jmax

      me = ipg_nodeid()
      nn = ipg_nnodes()

      call pg_distribution(handle, me, ilo, ihi, jlo, jhi)
c
      nx = ihi - ilo + 1
      ny = jhi - jlo + 1

      ibuff = igmem_alloc(max(1,nx*ny))

      call pg_get(handle,ilo,ihi,jlo,jhi,
     &              Q(ibuff),nx)

      result = 0.0d0
      imax = -1
      jmax = -1
      do i = max(ailo,ilo),min(aihi,ihi)
         do j = max(ajlo,jlo),min(ajhi,jhi)
            ix = i - ilo + 1
            jx = j - jlo + 1
            if(result .le. abs(
     &   Q(ibuff + nx*(jx-1) + (ix-1))))then
               result = abs(
     &   Q(ibuff + nx*(jx-1) + (ix-1)))
               imax = i
               jmax = j
            endif
         enddo
      enddo

      global_result = result
      call pg_dgop(1002,global_result, 1, 'absmax')

c
c The following is to ensure all node agree on the idomo/ivmo
c indices even allowing for rounding errors in the dgop
c

      do i = 0,nn-1
         if(i .eq. me .and. 
     &        abs(result - global_result) .lt. 1.0d-20)then
            indexi(1,i) = imax
            indexi(2,i) = jmax
         else
            indexi(1,i) = 0
            indexi(2,i) = 0
         endif
      enddo

      call pg_igop(1001,indexi, 2*nn, '+')

c Now all nodes take the indices with the highest node number

      do i = 0, nn-1
         if(indexi(1,i) .ne. 0)then
            idomo = indexi(1,i)
            ivmo = indexi(2,i)
         endif
      enddo

      call gmem_free(ibuff)

      gax_absmax_patch = global_result

      end
c
c parallel orthogonaliser, with apologies to rjh
c
c start with s matrix in the ao basis distributed in ga ih_ov
c porth will then take the input matrix of vectors v,
c construct the metric qdaggersq, and return orth'ed v
c

c     nu = no. of basis functions
c     ni = no. of vectors
c
c     overlap should be (nu,nu)
c
      subroutine porth(v,s,num,nvec)
      implicit none
      integer num, nvec
      REAL v(num,num), s(num,num)
      integer g_vecs, g_over, ulo, uhi, ld, uchunk
      integer k_tmp, k_s, k_over, k_w
      integer igmem_alloc_inf, ipg_nnodes, ipg_nodeid
      integer type, nu, ni
      logical opg_root
      character *4 fnm
      character *5 snm

INCLUDE(common/gmempara)
INCLUDE(common/gadiis)
INCLUDE(common/vcore)
INCLUDE(common/timeperiods)

      data fnm/'ga.m'/
      data snm/'porth'/

      g_vecs = ih_vec 
      g_over = ih_ovt

      if (supressed) call caserr('porth - no gas')

      call start_time_period(TP_GAORTHOG)
      call load_ga_from_square(ih_vec,v,num)
      call load_ga_from_triangle(ih_ovt,s,num)

c      write(6,*)'ga with input vecs'
c      call pg_print(g_vecs)

c      call mult2_ga(ih_ovt,ih_vec,ih_scr,num)
c      write(6,*)'overlap ga mo basis'
c      call pg_print(ih_ovt)


c
c     redistribute the input matrix ... block the leading
c     dimension, leave second dimension undistributed
c     ... each process has a(ulo:uhi,1:ni).  if
c

      call pg_inquire(g_vecs, type, nu, ni)
c
c     truncate ni to the actual number of vectors that need 
c     to be orthogonalised
c
      ni = min(nvec,ni)
      uchunk = max(8, (nu-1)/ipg_nnodes()+1)
      ulo  = ipg_nodeid()*uchunk + 1
      uhi  = min(ulo + uchunk - 1, nu)
      ld   = uhi - ulo + 1
      if (ulo .gt. uhi) then
         ulo = 0
         uhi = -1
         ld  = 1
      end if
c
c allocate local workspace
c
      k_tmp = igmem_alloc_inf(ld*ni,fnm,snm,'vecs',IGMEM_DEBUG)
      k_over = igmem_alloc_inf(ld*nu,fnm,snm,'metric',IGMEM_DEBUG)
      k_s = igmem_alloc_inf(ni,fnm,snm,'overlap',IGMEM_DEBUG)
      k_w = igmem_alloc_inf(nu,fnm,snm,'vec',IGMEM_DEBUG)
c
c  load by columns
c
      call pg_synch(1001)
      if (uhi .ge. ulo) then
         call pg_get(g_vecs, ulo, uhi, 1, ni,Q(k_tmp), ld)
	  endif
      if (uhi .ge. ulo) then
         call pg_get(g_over, 1, nu, ulo, uhi,Q(k_over), nu)
      end if

      call pg_synch(1002)

      call ga_orthog_vecs(Q(k_tmp), ld, ni, ulo, uhi, 
     &     Q(k_s), Q(k_over), nu, 
     &     Q(k_w), .true.)

      call pg_synch(1003)
c
c     put results back into to ga
c
      if (uhi .ge. ulo)call pg_put(g_vecs, ulo, uhi, 1, ni, 
     &     Q(k_tmp), ld)
c
c     tidy up memory
c
      call gmem_free_inf(k_w,fnm,snm,'vec')
      call gmem_free_inf(k_s,fnm,snm,'overlap')
      call gmem_free_inf(k_over,fnm,snm,'metric')
      call gmem_free_inf(k_tmp,fnm,snm,'vecs')

c this based on dgop
c maybe should compare with ga_get() of the whole array

      call load_square_from_ga(v,g_vecs,num)

c      if(opg_root())write(6,*)'result ga'
c      call pg_print(g_vecs)

      call end_time_period(TP_GAORTHOG)
      return
      end

      subroutine ga_orthog_vecs(vecs, ld, ni, ulo, uhi, s, o, nu, w, 
     $     ometric)
      implicit none
INCLUDE(common/sizes)
INCLUDE(common/harmon)
c     
      integer ld, ni, ulo, uhi, nu
      double precision 
     $     vecs(ulo:(ulo+ld-1),1:ni), 
     $     o(1:nu,ulo:(ulo+ld-1)), 
     $     s(ni),               
     $     w(nu)                
      logical ometric
c     
      integer i, j, u, npass
      double precision si, scale
      character*1 xn, xt
      integer count_zero
      data xn,xt/'n','t'/
c
      count_zero = 0
c     
c     orthogonalize columns of a matrix distributed so that
c     each process has vecs(ulo:uhi,1:ni) ... uses global sums only.
c     
      do i = 1, ni
         npass = 0
 10      npass = npass + 1
c     
c     if have a metric then first form overlap*vec(i)
c     
         if (ometric) then
            call dfill(nu, 0.0d0, w, 1)
 	    if ((uhi-ulo+1) .gt. 0) 
     $           call dgemv(xn, nu, (uhi-ulo+1) 
     $           ,1.0d0, o, nu
     $           ,vecs(ulo,i), 1, 0.0d0, w, 1)

            call pg_dgop(1001, w, nu, '+')
         else
            do u = ulo, uhi
               w(u) = vecs(u,i)
            end do
         end if
c     
c     now form overlap between vector i and vectors 1...i
c     
         call dfill(i, 0.0d0, s, 1)
         if ((uhi-ulo+1) .gt. 0) 
     $        call dgemv(xt, (uhi-ulo+1) 
     $        ,i, 1.0d0
     $        ,vecs(ulo,1), ld, w(ulo)
     $        ,1, 0.0d0, s(1), 1)

         call pg_dgop(1002, s, i, '+')

c
c     apply the rotation
c     
         if ((uhi-ulo+1).gt.0 .and. i.gt.1) then
            call dgemv(xn, (uhi-ulo+1), i-1, -1.0d0
     $           ,vecs(ulo,1), ld, s(1)
     $           ,1, 1.0d0, vecs(ulo,i), 1)
         end if
c     
c     renormalize vector i
c     
         si = s(i)
         do j = 1, i-1
            si = si - s(j)*s(j)
         end do
c     
c     if the vector norm changed a lot then repeat
c     
         if (i .gt. 1) then
            scale = si/s(i)
            if (scale .lt. 0.9d0) then
               if (npass .lt. 3) then
                  goto 10
               else
                  call caserr('ga_orthog: failed to orthog vector')
               end if
            end if
         end if
c
         if (si .lt. 0.0d0) call caserr('ga_orthog: negative')
         if (si .eq. 0.0d0) then
            count_zero = count_zero + 1
            if (count_zero.gt.newbas1-newbas0)
     1      call caserr('ga_orthog: hard zero')
            scale = 0.0d0
         else
            scale = 1.0d0/sqrt(si)
         end if
         do u = ulo, uhi
            vecs(u,i) = vecs(u,i) * scale
         end do
c     
      end do
c
      if (count_zero.ne.newbas1-newbas1) 
     1   call caserr('# zeros in ga_orthog < harmon-caused')
c
      end

**********************************************************************
*
      block data ga_file_initialize
*
*  block data to initialize the virtual file system which
* is implemented by global array (ga) tools.
*
      implicit none
*
*  all the relevant stuff is in an include file to avoid cock-ups.
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*

*  The default value for the amount of memory each node uses
* for the virtual file system. Best to make this a multiple
* of the block size. Default values are:
* 64 MByes nodes: just over 4 Mbytes per process shared over all 
* GAs (1000*511*8) 
* 128 MByte nodes: just over 12 Mbytes per process shared over all 
* GAs (3000*511*8) 
*
      Integer    ga_file_proc_mem_default
_IF(edgafs)
      Parameter( ga_file_proc_mem_default = 30000 * ga_file_block )
_ELSEIF(t3e,rs6000,r10000,ev6,ev5,linux)
      Parameter( ga_file_proc_mem_default = 3000 * ga_file_block )
_ELSE
      Parameter( ga_file_proc_mem_default = 1000 * ga_file_block )
_ENDIF
*
*  a silly parameter to help the data statements.
*
      integer    unit_range
      parameter( unit_range = max_unit_number - min_unit_number + 1 )
*
*  initialize to `sensible' values. in fact most of the values are
* silly to alllow the code to check for initialization errors. 
*
      data ga_file_initialized           / .false.                   /
      data ga_file_proc_mem              / ga_file_proc_mem_default  /
      data ga_file_total, ga_file_unused / not_in_use, not_in_use    /
      data ga_file_processes             / not_in_use                /
      data ga_file_proc_id               / not_in_use                /
      data ga_file_root                  / .false.                   /
      data ga_file_handle                / max_ga_files * not_in_use /
      data ga_file_open                  / max_ga_files * .false.    /
      data ga_file_size                  / max_ga_files * not_in_use /
      data ga_file_next                  / max_ga_files * not_in_use /
      data ga_file_mode                  / max_ga_files * .true.     /
      data ga_file_accessed              / max_ga_files * .false.    /
      data ga_file_leng1                 / max_ga_files * not_in_use /
      data ga_file_hist_next             / max_ga_files * 1          /
      data ga_file_unit_to_ga            / unit_range   * not_in_use /
      data ga_file_safety_mode           / ga_file_safe              /
*
      end
*
**********************************************************************
*
      subroutine wrt3_ga( array, words, block, unit )
*
*     subroutine to simulate wrt3 by using ga tools instead of disk.
*
      implicit none
*
      integer          words
      double precision array( 1:words )
      integer          block
      integer          unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/timeperiods)
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/iofile)
INCLUDE(common/errcodes)
*
*  functions
*
      logical  get_oswed3
      external get_oswed3
      logical  pg_locate_region
      external pg_locate_region
      logical  opg_root
      external opg_root
*
*     local variables:
*
      integer ga_to_use
      integer start, finish
      integer num_blocks
      integer hist_next
      integer map( 1:5, 1:max_processes ), n_owners
      MA_INTEGER work( 1:5, 1:max_processes )
      integer my_start, my_finish, my_words
      integer i
      logical worked
      logical do_put
      integer il,ih,jl,jh
      character*4 yga

*
      if (iwhere(unit).eq.6)then
       call start_time_period(TP_IOP6)
      else
       call start_time_period(TP_IOP7)
      endif
*
*     find out which ga is mapped to this unit number
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*     if the file is not open open it.
*
*
*
      if( .not. ga_file_open( ga_to_use ) ) then
         call rdtpnm_ga( unit )
      end if
*
*     check that the filing system has been initialized. it should be
*     but you never know.
*
      if( .not. ga_file_initialized ) then

         call gamerr('wrt3: ga not initalized - internal error',
     &        ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

      end if
*
*     if the main code changed the file io mode update the
*     local mode maintainer.
*
      ga_file_mode( ga_to_use ) = get_oswed3( unit )
*
*     if this is the first block (ie for dumpfiles the index block),
*     save the length. we may need to use the get_ga call to recover
*     the index block later (from secini).
*
      if( block .eq. 1 ) ga_file_leng1( ga_to_use ) = words
*
*
*     find out where this block starts
*
      start = ( block - 1 ) * ga_file_block + 1
*
*     and where it finishes. remember to catch zero word writes
*     which simply update the next block counter.
*
      finish = start + max( words, 1 ) - 1
c
c     special case for ed19 c/z vectors
c
      jl=1
      jh=1
      if (unit.eq.20) then
          if ( words .eq. 0 ) goto 999
          call get_ed19adres(start,finish,il,ih,jl,jh)
          start=il
          finish=ih
      endif
*
*     update the next block data on this file. again trap 
*     zero word writes.
*
      num_blocks = ( max( words, 1 ) - 1 ) / ga_file_block + 1
      ga_file_next( ga_to_use ) = block + num_blocks
      call set_ipos( ga_file_next( ga_to_use ), unit )
*
*     update the history array. 
*     ( this need only be done if the file is to be dumped at the 
*        end of the run. )
*     however we now always look after a history to allow the correct
*     value of the high water mark for the ga file to be found.
*     uncomment the if block below if you don't want the ( probably 
*     small ) overhead associated with this.
*
      if( ga_file_keep( ga_to_use ) ) then
*
         hist_next = ga_file_hist_next( ga_to_use )
         ga_file_history( 1, hist_next, ga_to_use ) = words
         ga_file_history( 2, hist_next, ga_to_use ) = block
         ga_file_history( 3, hist_next, ga_to_use ) = 
     +        block + num_blocks - 1
         ga_file_hist_next( ga_to_use ) = hist_next + 1
*
*        if about to fall off the end of the history compress it.
*
         if( hist_next + 1 .ge. max_ga_calls ) then
            call ga_file_history_compress( unit )
*
*           check have not truely run out of memory
*
            if( ga_file_hist_next( ga_to_use ) .gt. max_ga_calls ) then
               
               call gamerr('wrt3: run out of history space: '//
     &              yed(unit),
     &              ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

            end if
*
         end if
*
      end if
*
*     can skip the rest for zero length writes.
*
      if( words .gt. 0 ) then
*
*        check that we have enough memory to save this.
*
         if( finish .gt. ga_file_size( ga_to_use ) ) then
*
*           we do not. for a production version we should spill to disk,
*           but at present just flag an error.
*
            if (opg_root()) then
               write(iwr,*)'wrt3: current GA file         = '//
     &                     yed(unit)
               write(iwr,*)'wrt3: maximum size of GA file = ',
     &                     ga_file_size( ga_to_use ),' = ',
     &                     ga_file_size( ga_to_use )/512,' blocks '
               write(iwr,*)'wrt3: number of words to save = ',words
               write(iwr,*)'wrt3: final size of GA file   = ',finish,
     &                     ' = ',finish/512,' blocks '
               write(iwr,"(1x,120a)")
     &                   'wrt3: Final size of GA file exceeds the '//
     &                   'maximum size. Please increase the GA file '//
     &                   'in size !!!'
            endif
            call gamerr(
     &           'wrt3: insufficient space in global array: '//
     &           yed(unit),
     &           ERR_NO_CODE, ERR_USER_DIMENSION, ERR_ASYNC, ERR_NO_SYS)

         else
*
*           we do have enough memory, whack the array into it if we are
*           the root node. only one node ever does the `writing' to a
*           given area of memory to ensure `file' coherency.
*
            if (oadd(unit)) then
               call pg_acc( ga_file_handle( ga_to_use ),
     +                         start, finish,
     +                         jl, jh,
     +                         array,
     +                         words,1.0d0 )
            elseif( ga_file_mode( ga_to_use ) .and. COLLECTIVE ) then
*
*              working in colective io mode. as the data is replicated
*              the processor that has the relevant part of the ga in
*              local memory writes the corresponding part of the array.
*
               worked = pg_locate_region( ga_file_handle( ga_to_use ), 
     +                                    start, finish,
     +                                    jl, jh, map, n_owners, work)
               if( .not. worked ) then

                  call gamerr(
     &         'wrt3: failed to locate desired data',
     &         ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

               end if
*
*              search the map array to see if i own some, or possibly
*              all, of the ga.
*
               do_put = .false.
               do i = 1, n_owners
                  if( map( 5, i ) .eq. ga_file_proc_id ) then
                     my_start  = map( 1, i )
                     my_finish = map( 2, i )
                     my_words  = my_finish - my_start + 1
                     do_put = .true.
                     go to 10
                  end if
               end do
*
 10            continue
*
*              if i do own some put the corresponding bit of the array
*              into it.
*
               if( do_put ) then
                  call pg_put( ga_file_handle( ga_to_use ), 
     +                         my_start, my_finish,
     +                         jl, jh, 
     +                         array( my_start - start + 1 ), 
     +                         my_words )
               end if
*
*              if safe mode ensure all transactions have completed, 
*              to ensure that another processor does not try to read
*              what has just been put before the put completes.
*
               if( ga_file_safety_mode .eq. ga_file_safe ) then
                  call pg_synch( 160467 )
               end if
*
            else
*
*              working in single node io mode. only the root processor
*              should come through here.
*
               if( ga_file_root ) then
                  call pg_put( ga_file_handle( ga_to_use ), 
     +                         start, finish,
     +                         jl, jh, array, words )
               else
                  print*, 'wrt3: non root node attempting ' //
     +                    'single node io'
                  


               end if
*
            end if
*
         end if
*
      end if
*
*     finally indicate that something has happened to this file.
*
      ga_file_accessed( ga_to_use ) = .true.
*
 999  continue
      if (iwhere(unit).eq.6)then
       call end_time_period(TP_IOP6)
      else
       call end_time_period(TP_IOP7)
      endif

      end
*
***********************************************************************
*
      subroutine rdedx_ga( array, words, block, unit )
*
*     subroutine to simulate rdedx by using ga tools instead of disk.
*
      implicit none
*
      integer          words
      double precision array( 1:words )
      integer          block
      integer          unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/timeperiods)
INCLUDE(common/errcodes)
*
*     functions:
*
      logical  get_oswed3
      external get_oswed3
      logical  pg_locate_region
      external pg_locate_region
*
*     local variables:
*
      integer ga_to_use
      integer start, finish
      integer num_blocks
      logical present_mode
      logical worked
      integer map( 1:5, 1:max_processes ), n_owners, owner
      MA_INTEGER work( 1:5, 1:max_processes )
      integer left, out
      integer il,ih,jl,jh
*
      if (iwhere(unit).eq.6)then
       call start_time_period(TP_IOG6)
      else
       call start_time_period(TP_IOG7)
      endif
*
*     find out which ga is associated with this unit number.
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*     check that the filing system has been initialized
*
      if( .not. ga_file_initialized ) then

                  call gamerr(
     &         'rdedx: ga not initalized',
     &         ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

      end if
*
*     if this is the first access to this file it is assumed that
*     this is a restart file. open it if required, and then fill
*     it with data from the file(s) on the disk.
*
      if( .not. ga_file_accessed( ga_to_use ) ) then
         if( .not. ga_file_open( ga_to_use ) ) then
            call rdtpnm_ga( unit )
         end if
         call read_restart( unit, zedfil(unit) )
      end if
*
*     if this is not the first acces to this file check that the file 
*     is open. - unlike the write case the
*     file should not be opened `cos you can't read from an empty file.
*
      if( .not. ga_file_open( ga_to_use ) ) then

         print*, 'rdedx: unit ', unit, ' not opened'

                  call gamerr(
     &         'ga file internal error',
     &         ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

      end if
*
*     have to be a little careful here. if the main code has changed the
*     io mode from single node to collective there may still be
*     outstanding writes. hence a read of might result in an old or
*     corrupted result. hence the sync. also check that if doing single
*     node io only the root comes through the routine.
*
      present_mode = get_oswed3( unit )
      if( .not. present_mode .and. COLLECTIVE.and.unit.ne.20 ) then
         if( .not. ga_file_root ) then
            
            call gamerr(
     &           'rdedx: non root node attempting single node io',
     &           ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

         end if
      end if
      if( ( present_mode         .and. COLLECTIVE ) .and. .not.
     +    ( ga_file_mode( unit ) .and. COLLECTIVE ) ) then
         call pg_synch( 160467 )
      end if
*
      ga_file_mode( ga_to_use ) = present_mode
*
*     find out where this block starts
*
      start = ( block - 1 ) * ga_file_block + 1
*
*     and where it finishes, trapping zero word reads.
*
      finish = start + max( words, 1 ) - 1
c
c     special case ed19 c/z vectors ...
c
      jl=1
      jh=1
      if (unit.eq.20) then
         if ( words .eq. 0 ) goto 999
         call get_ed19adres(start,finish,il,ih,jl,jh)
         start=il
         finish=ih
      endif
c
*
*     update the next block data on this file, againg trapping
*     zero word reads.
*
      num_blocks = ( max( words, 1 ) - 1 ) / ga_file_block + 1
      ga_file_next( ga_to_use ) = block + num_blocks
      call set_ipos( ga_file_next( ga_to_use ), unit )
*
*     if this is a zero word read can skip the rest.
*
      if( words .gt. 0 ) then
*
*        check if this is all stored in memory.
*
         if( finish .gt. ga_file_size( ga_to_use ) ) then
*
*           we do not. for a production version it should have been
*           spilled to disk, but at present just flag an error.
*
            call gamerr(
     &          'rdedx: insufficient space in global arrays',
     &          ERR_NO_CODE, ERR_USER_DIMENSION, ERR_ASYNC, ERR_NO_SYS)

*
         else
*
*           it is all in memory. the method used to get it is 
*           a) if we are in collective io mode:
*           1) locate on which node(s) the data is located
*           2) if only one node owns all the data, that node broadcasts
*              it to all the other nodes
*           3) otherwise the root node gathers all the data together and
*              the broadcasts it.
*           this was found to be more efficient than all nodes getting,
*           contention probably being the reason.
*           b) if we are in single node io mode the root node just
*              gets the data
*
            if( ga_file_mode( ga_to_use ) .and. COLLECTIVE ) then
*
*              collective io
*              see next comment for thoughts on efficiency.
*
               worked = pg_locate_region( ga_file_handle( ga_to_use ), 
     +                                    start, finish,
     +                                    jl, jh, map, n_owners, work)
               if( .not. worked ) then

                  call gamerr('rdedx: failed to locate desired data',
     &                 ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

               end if
*
*              all of what now follows could probably be speeded up by
*              by having some sort of n_owners headed tree rather than
*              gather on one and then broadcast.
*
*
*              hack alert !
*
*              burnt in zero below represents the root node - may have 
*              to rework this.
*
               if( n_owners .eq. 1 ) then
                  owner = map( 5, 1 )
               else
                  owner = 0
               end if
*
               if( ga_file_proc_id .eq. owner ) then
                  call pg_get( ga_file_handle( ga_to_use ), 
     +                         start, finish, 
     +                         jl, jh, array, words )
               end if
*
*              hack alert !
*
*              burnt in 8 is the size in bytes of a double precision
*              real.
*
               if( ga_file_processes .ne. 1 ) then
                  left  = words
                  start = 1
                  do while( left .gt. 0 )
                     out = min( 16384, left )
                     call pg_brdcst( 160467, array( start ), 
     +                               out * 8, owner )
                     start = start + out
                     left  = left  - out
                  end do
               end if
*
            else
*
*              single node io mode
*
               call pg_get( ga_file_handle( ga_to_use ), 
     +                      start, finish, 
     +                      jl, jh, array, words )
*     
            end if
*
         end if
*
      end if
*
*     and finally indicate that something has happened to this file.
*
      ga_file_accessed( ga_to_use ) = .true.
*
 999  continue
      if (iwhere(unit).eq.6)then
       call end_time_period(TP_IOG6)
      else
       call end_time_period(TP_IOG7)
      endif

      end
*
***********************************************************************
*
      subroutine rdtpnm_ga( unit )
*
*  subroutine to simulate rdtpnm by using ga tools instead of disk.
*
      implicit none
*
      integer          unit
*
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
INCLUDE(common/ga_file)
INCLUDE(common/timeperiods)
INCLUDE(common/disc)
INCLUDE(common/discc)
INCLUDE(common/errcodes)
INCLUDE(common/iofile)
*
*  functions:
*
      logical  get_oswed3
      external get_oswed3
      logical  pg_create_inf
      external pg_create_inf
      logical  opg_root
      external opg_root
      integer idim,jdim,col_chunk

_IF(charmm,chemshell)
*
* Ensure loading of block data (cray requirement)
*
      external ga_file_initialize
_ENDIF
*
*  local variables:
*
      integer     ga_to_use
      integer     crap
      integer     row_chunk
      logical     worked

      if (iwhere(unit).eq.6)then
       call start_time_period(TP_IOO6)
      else
       call start_time_period(TP_IOO7)
      endif

*
*  opening a ga file must be done collectively. check that the file
* mode is correct.
*
      if( .not. get_oswed3( unit ) ) then

         print*, 'rdtpnm: attempt to open unit ', unit, ' in ' //
     +           'single node io mode.'

         call gamerr('ga internal error',
     &        ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

      end if
*
*  check that this is a valid unit number. it should have been mapped
* to a ga unit number by ga_file_set_up
*
      ga_to_use = ga_file_unit_to_ga( unit )
      if( ga_to_use .eq. not_in_use ) then
         print*, 'rdtpnm: invalid unit number ', unit

         call gamerr('ga internal error',
     &        ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

      end if
*
*  check that the file is not already open.
*
      if( ga_file_open( ga_to_use ) ) then
         print*, 'rdtpnm: unit ', unit, ' already open.'

         call gamerr('ga internal error',
     &        ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)
         
      end if
*
*  now open the file by creating a ga and setting up the various
* associated data things.
*
*
*  first of all update the memory counters.
*
      ga_file_total  = ga_file_total  + ga_file_size( ga_to_use )
      ga_file_unused = ga_file_unused - ga_file_size( ga_to_use )
*
*  next create the ga beast.
*
*  ga wants the minimum size for a chunk. hence the following is fine
* for the number of rows in a chunk on one processor.
*
      row_chunk = ga_file_size( ga_to_use ) / ga_file_processes
*
*  note that the gamess wrapper burns in the type parameter to
* ga_create - hence crap in the function reference below
* can be set to any old crap.
*
      if(opg_root())write(iwr,*)'allocating ga file memory: ',
     &   ga_file_size( ga_to_use ),' distributed words'
      jdim=1
      idim=ga_file_size(ga_to_use)
      col_chunk=1
      if (unit.eq.20) then
         jdim=2
         idim=ga_file_size(ga_to_use)/2
         row_chunk=ga_file_block
         col_chunk=2
      endif
      worked = pg_create_inf( crap, idim, jdim, 
     +                    yed(unit), row_chunk, col_chunk, 
     +                    ga_file_handle( ga_to_use ),'ga.m', 
     +                    'rdtpnm_ga',IGMEM_NORMAL)
c     worked = pg_create( crap, ga_file_size( ga_to_use ), 1, 
c    +                    yed(unit), row_chunk, 1, 
c    +                    ga_file_handle( ga_to_use ) )
*
      if( .not. worked ) then

         print*, 'rdtpnm: failed to create GA file ',yed(unit)
         print*, 'allocate more core memory or use more nodes'

         call gamerr('ga internal error',
     &        ERR_NO_CODE, ERR_USER_DIMENSION, ERR_ASYNC, ERR_NO_SYS)
         
      end if
      if (oadd(unit)) call pg_zero(ga_file_handle( ga_to_use ))
*
*  set flag indicating that the file is open
*
      ga_file_open( ga_to_use ) = .true.
*
*  `rewind' the file. for certain parts of gamess one of
* its internal data structures needs updating.
*
      ga_file_next( ga_to_use ) = 1
      call set_ipos( ga_file_next( ga_to_use ), unit )
*
*  and set the io mode.
*
      ga_file_mode( ga_to_use ) = COLLECTIVE
*
*  ensure all transactions have completed
*
      call pg_synch( 160467 )

      if (iwhere(unit).eq.6)then
       call end_time_period(TP_IOO6)
      else
       call end_time_period(TP_IOO7)
      endif

*
      end
*
**********************************************************************
*
      subroutine wrt3s_ga( array, words, unit )
*
*  subroutine to simulate wrt3s by using ga tools instead of disk.
*
      implicit none
*
      integer          words
      double precision array( 1:words )
      integer          unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local variables
*
      integer ga_to_use         
      integer block
*
*  find out which ga this unit is mapped to
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  find out which block to write to
*
      block = ga_file_next( ga_to_use )
*
*  and just use wrt3 to do the writing
*
      call wrt3_ga( array, words, block, unit )
*
      end 
*
***********************************************************************
*
      subroutine wrt3i_ga( array, words, block, unit )
*
*  subroutine to simulate wrt3i by using ga tools instead of disk.
*
      implicit none
*
      integer words
      integer array( 1:words )
      integer block
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  functions:
*
      integer  lenwrd
      external lenwrd
*
*  local static variables
*
      integer ints_per_real
      save    ints_per_real
      integer ints_in_buffer
      save    ints_in_buffer
      logical first
      save    first
*
*  local variables
*
      integer num_ints
      integer where, ints_to_write, reals_to_write
      integer old_safety
      integer i
      logical force_loop
*
*  local initialization
*
      data first / .true. /
*
*  the method used is to pack the integers into a double precision
* array and then use the double precision write. to avoid
* too much memory usage the writes are blocked. for this to
* work properly the parameter ga_file_buf_len in the include
* file must be an integer multiple of ga_file_block. check this !
*
*  first of all initialize on first call. need to know how the relative
* sizes of integers and double precisions, and can then set up the
* buffer length in terms of integers
*
      if( first ) then
         first          = .false.
         ints_per_real  = lenwrd()
         ints_in_buffer = ga_file_buf_len * ints_per_real
      end if
*
*  can save a lot of calls to the sychronization routine by
* changing the safety mode that the write runs under.
*
      old_safety = ga_file_safety_mode
      call set_file_safety( ga_file_risky, unit )
*
* num_ints looks after how many integers are yet to be written,
* where points to the next section of the array to be written.
*
      num_ints = words
      where    = 1
*
*  have to go through the loop at least once to cope with zero word 
* writes.
*
      force_loop = .true.
*
      do while( num_ints .gt. 0 .or. force_loop )
*
*  find out how much to write in terms of both reals
* and double precisions. we want to write as much as possible.
*
         ints_to_write  = min( num_ints, ints_in_buffer )
         reals_to_write = ( ints_to_write + ints_per_real - 1 ) / 
     +                    ints_per_real
*
*  pack the integer buffer.
*
         do i = 1, ints_to_write
            ga_file_integer_buf( i ) = array( i + where - 1 )
         end do
*
*  and write the double precision buffer that is equivalenced to
* the integer buffer. only need to specify a block number
* on the first call.
*
         if( where .eq. 1 ) then
            call wrt3_ga( ga_file_double_buf, reals_to_write, 
     +                  block, unit )
         else
            call wrt3s_ga( ga_file_double_buf, reals_to_write, unit )
         end if
*
*  and update the counters.
*
         num_ints = num_ints - ints_to_write
         where    = where    + ints_to_write
*
         force_loop = .false.
*
      end do
*
*  and reset the file safety mode
*
      call set_file_safety( old_safety, unit )
*
      end
*
***********************************************************************
*
      subroutine wrt3is_ga( array, words, unit )
*
*  subroutine to simulate wrt3is by using ga tools instead of disk.
*
      implicit none
*
      integer words
      integer array( 1:words )
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local variables
*
      integer ga_to_use         
      integer block
*
*  find out which ga this unit is mapped to
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  find out which block to write to
*
      block = ga_file_next( ga_to_use )
*
*  and just use wrt3i, and so eventually wrt3, to do the writing
*
      call wrt3i_ga( array, words, block, unit )
*
      end 
*
***********************************************************************
*
      subroutine wrtc_ga( array, words, block, unit )
*
*  subroutine to simulate wrtc by using ga tools instead of disk.
*
      implicit none
*
      integer     words
      character*8 array( 1:words )
      integer     block
      integer     unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local variables
*
      integer old_safety
      integer num_chars
      integer where
      integer chars_to_write
      integer i
      logical force_loop
*
*  the method used is to pack the integers into a double precision
* array and then use the double precision write. to avoid
* too much memory usage the writes are blocked. for this to
* work properly the parameter ga_file_buf_len in the include
* file must be an integer multiple of ga_file_block. check this !
*
*
*  can save a lot of calls to the sychronization routine by
* changing the safety mode that the write runs under.
*
      old_safety = ga_file_safety_mode
      call set_file_safety( ga_file_risky, unit )
*
* num_chars looks after how many integers are yet to be written,
* where points to the next section of the array to be written.
*
      num_chars = words
      where     = 1
*
*  have to go through loop once to cope with zero words writes.
*
      force_loop = .true.
*
      do while( num_chars .gt. 0 .or. force_loop )
*
*  find out how much to write. we want to write as much as possible.
*
         chars_to_write = min( num_chars, ga_file_buf_len )
*
*  pack the buffer.
*
         do i = 1, chars_to_write
            read( array( i + where - 1 ), '( a8 )' ) 
     +            ga_file_double_buf( i )
         end do
*
*  and write the double precsion buffer.
* only need to specify a block number on the first call.
*
         if( where .eq. 1 ) then
            call wrt3_ga( ga_file_double_buf, chars_to_write, 
     +                    block, unit )
         else
            call wrt3s_ga( ga_file_double_buf, chars_to_write, unit )
         end if
*
*  and update the counters.
*
         num_chars = num_chars - chars_to_write
         where     = where     + chars_to_write
*
         force_loop = .false.
*
      end do
*
*  and reset the file safety mode
*
      call set_file_safety( old_safety, unit )
*
      end
*
***********************************************************************
*
      subroutine wrtcs_ga( array, words, unit )
*
*  subroutine to simulate wrtcs by using ga tools instead of disk.
*
      implicit none
*
      integer     words
      character*8 array( 1:words )
      integer     unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local variables
*
      integer ga_to_use         
      integer block
*
*  find out which ga this unit is mapped to
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  find out which block to write to
*
      block = ga_file_next( ga_to_use )
*
*  and just use wrtc, and so eventually wrt3, to do the writing
*
      call wrtc_ga( array, words, block, unit )
*
      end 
*
***********************************************************************
*
      subroutine reads_ga( array, words, unit )
*
*  subroutine to simulate read3s by using ga tools instead of disk.
*
      implicit none
*
      integer          words
      double precision array( 1:words )
      integer          unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local variables
*
      integer ga_to_use         
      integer block
*
*  find out which ga this unit is mapped to
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  find out which block to read from
*
      block = ga_file_next( ga_to_use )
*
*  and just use rdedx to do the reading
*
      call rdedx_ga( array, words, block, unit )
*
      end 
*
***********************************************************************
*
      subroutine readi_ga( array, words, block, unit )
*
*  subroutine to simulate readi by using ga tools instead of disk.
*
      implicit none
*
      integer words
      integer array( 1:words )
      integer block
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  functions:
*
      integer  lenwrd
      external lenwrd
*
*  local static variables
*
      integer ints_per_real
      save    ints_per_real
      integer ints_in_buffer
      save    ints_in_buffer
      logical first
      save    first
*
*  local variables
*
      integer num_ints
      integer where, ints_to_read, reals_to_read
      integer i
      logical force_loop
*
*  local initialization
*
      data first / .true. /
*
*  the method used is to pack the integers into a double precision
* array and then use the double precision read. to avoid
* too much memory usage the reads are blocked. for this to
* work properly the parameter ga_file_buf_len in the include
* file must be an integer multiple of ga_file_block. check this !
*
*  first of all initialize on first call. need to know how the relative
* sizes of integers and double precisions, and can then set up the
* buffer length in terms of integers
*
      if( first ) then
         first          = .false.
         ints_per_real  = lenwrd()
         ints_in_buffer = ga_file_buf_len * ints_per_real
      end if
*
* num_ints looks after how many integers are yet to be read,
* where points to the next section of the array to be read.
*
      num_ints = words
      where    = 1
*
*  have to go through loop once to cope with zero word reads
*
      force_loop = .true.
*
      do while( num_ints .gt. 0 .or. force_loop )
*
*  find out how much to read in terms of both reals
* and double precsions. we want to read as much as possible.
*
         ints_to_read  = min( num_ints, ints_in_buffer )
         reals_to_read = ( ints_to_read + ints_per_real - 1 ) / 
     +                     ints_per_real
*
*  read into the double precision buffer that is equivalenced to
* the integer buffer. only need to specify a block number
* on the first call.
*
         if( where .eq. 1 ) then
            call rdedx_ga( ga_file_double_buf, reals_to_read, 
     +                     block, unit )
         else
            call reads_ga( ga_file_double_buf, reals_to_read, unit )
         end if
*
*  unpack into the integer buffer.
*
         do i = 1, ints_to_read
            array( i + where - 1 ) = ga_file_integer_buf( i ) 
         end do
*
*  and update the counters.
*
         num_ints = num_ints - ints_to_read
         where    = where    + ints_to_read
*
         force_loop = .false.
*
      end do
*
      end
*
***********************************************************************
*
      subroutine readis_ga( array, words, unit )
*
*  subroutine to simulate readis by using ga tools instead of disk.
*
      implicit none
*
      integer words
      integer array( 1:words )
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local variables
*
      integer ga_to_use         
      integer block
*
*  find out which ga this unit is mapped to
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  find out which block to read from
*
      block = ga_file_next( ga_to_use )
*
*  and just use readi, and so eventually rdedx, to do the reading
*
      call readi_ga( array, words, block, unit )
*
      end 
*
***********************************************************************
*
      subroutine rdchr_ga( array, words, block, unit )
*
*  subroutine to simulate rdchr by using ga tools instead of disk.
*
      implicit none
*
      integer     words
      character*8 array( 1:words )
      integer     block
      integer     unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local variables
*
      integer num_chars
      integer where
      integer chars_to_read
      integer i
      logical force_loop
*
*  the method used is to pack the integers into a double precision
* array and then use the double precision read. to avoid
* too much memory usage the reads are blocked. for this to
* work properly the parameter ga_file_buf_len in the include
* file must be an integer multiple of ga_file_block. check this !
*
*
* num_chars looks after how many integers are yet to be read,
* where points to the next section of the array to be read.
*
      num_chars = words
      where     = 1
*
*  have to go through loop once to cope with zero word reads
*
      force_loop = .true.
*
      do while( num_chars .gt. 0 .or. force_loop )
*
*  find out how much to write. we want to write as much as possible.
*
         chars_to_read = min( num_chars, ga_file_buf_len )
*
*  read into the double precsion buffer.
* only need to specify a block number on the first call.
*
         if( where .eq. 1 ) then
            call rdedx_ga( ga_file_double_buf, chars_to_read, 
     +                  block, unit )
         else
            call reads_ga( ga_file_double_buf, chars_to_read, unit )
         end if
*
*  unpack the buffer.
*
         do i = 1, chars_to_read
            write( array( i + where - 1 ), '( a8 )' ) 
     +             ga_file_double_buf( i )
         end do
*
*  and update the counters.
*
         num_chars = num_chars - chars_to_read
         where     = where     + chars_to_read
*
         force_loop = .false.
*
      end do
*
      end
*
***********************************************************************
*
      subroutine rdchrs_ga( array, words, unit )
*
*  subroutine to simulate rdchrs by using ga tools instead of disk.
*
      implicit none
*
      integer     words
      character*8 array( 1:words )
      integer     unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local variables
*
      integer ga_to_use         
      integer block
*
*  find out which ga this unit is mapped to
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  find out which block to read from
*
      block = ga_file_next( ga_to_use )
*
*  and just use rdchr, and so eventually rdedx, to do the writing
*
      call rdchr_ga( array, words, block, unit )
*
      end 
*
_IF(i8drct)
***********************************************************************
*
      subroutine readi8_ga( array, words, block, unit )
*
*  subroutine to simulate readi8 by using ga tools instead of disk.
*
      implicit none
*
      integer words
      integer *8 array( 1:words )
      integer block
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local static variables
*
      integer ints_per_real
      save    ints_per_real
      integer ints_in_buffer
      save    ints_in_buffer
      logical first
      save    first
*
*  local variables
*
      integer num_ints
      integer where, ints_to_read, reals_to_read
      integer i
      logical force_loop
*
*  local initialization
*
      data first / .true. /
*
*  the method used is to pack the integers into a double precision
* array and then use the double precision read. to avoid
* too much memory usage the reads are blocked. for this to
* work properly the parameter ga_file_buf_len in the include
* file must be an integer multiple of ga_file_block. check this !
*
*  first of all initialize on first call. need to know how the relative
* sizes of integers and double precisions, and can then set up the
* buffer length in terms of integers
*
      if( first ) then
         first          = .false.
         ints_per_real  = 1
         ints_in_buffer = ga_file_buf_len * ints_per_real
      end if
*
* num_ints looks after how many integers are yet to be read,
* where points to the next section of the array to be read.
*
      num_ints = words
      where    = 1
*
*  have to go through loop once to cope with zero word reads
*
      force_loop = .true.
*
      do while( num_ints .gt. 0 .or. force_loop )
*
*  find out how much to read in terms of both reals
* and double precsions. we want to read as much as possible.
*
         ints_to_read  = min( num_ints, ints_in_buffer )
         reals_to_read = ( ints_to_read + ints_per_real - 1 ) / 
     +                     ints_per_real
*
*  read into the double precision buffer that is equivalenced to
* the integer buffer. only need to specify a block number
* on the first call.
*
         if( where .eq. 1 ) then
            call rdedx_ga( ga_file_double_buf, reals_to_read, 
     +                     block, unit )
         else
            call reads_ga( ga_file_double_buf, reals_to_read, unit )
         end if
*
*  unpack into the integer buffer.
*
         do i = 1, ints_to_read
            array( i + where - 1 ) = ga_file_integer_buf( i ) 
         end do
*
*  and update the counters.
*
         num_ints = num_ints - ints_to_read
         where    = where    + ints_to_read
*
         force_loop = .false.
*
      end do
*
      end
*
***********************************************************************
*
      subroutine wrt3i8_ga( array, words, block, unit )
*
*  subroutine to simulate wrt3i8 by using ga tools instead of disk.
*
      implicit none
*
      integer words
      integer *8 array( 1:words )
      integer block
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  local static variables
*
      integer ints_per_real
      save    ints_per_real
      integer ints_in_buffer
      save    ints_in_buffer
      logical first
      save    first
*
*  local variables
*
      integer num_ints
      integer where, ints_to_write, reals_to_write
      integer old_safety
      integer i
      logical force_loop
*
*  local initialization
*
      data first / .true. /
*
*  the method used is to pack the integers into a double precision
* array and then use the double precision write. to avoid
* too much memory usage the writes are blocked. for this to
* work properly the parameter ga_file_buf_len in the include
* file must be an integer multiple of ga_file_block. check this !
*
*  first of all initialize on first call. need to know how the relative
* sizes of integers and double precisions, and can then set up the
* buffer length in terms of integers
*
      if( first ) then
         first          = .false.
         ints_per_real  = 1
         ints_in_buffer = ga_file_buf_len * ints_per_real
      end if
*
*  can save a lot of calls to the sychronization routine by
* changing the safety mode that the write runs under.
*
      old_safety = ga_file_safety_mode
      call set_file_safety( ga_file_risky, unit )
*
* num_ints looks after how many integers are yet to be written,
* where points to the next section of the array to be written.
*
      num_ints = words
      where    = 1
*
*  have to go through the loop at least once to cope with zero word 
* writes.
*
      force_loop = .true.
*
      do while( num_ints .gt. 0 .or. force_loop )
*
*  find out how much to write in terms of both reals
* and double precisions. we want to write as much as possible.
*
         ints_to_write  = min( num_ints, ints_in_buffer )
         reals_to_write = ( ints_to_write + ints_per_real - 1 ) / 
     +                    ints_per_real
*
*  pack the integer buffer.
*
         do i = 1, ints_to_write
            ga_file_integer_buf( i ) = array( i + where - 1 )
         end do
*
*  and write the double precision buffer that is equivalenced to
* the integer buffer. only need to specify a block number
* on the first call.
*
         if( where .eq. 1 ) then
            call wrt3_ga( ga_file_double_buf, reals_to_write, 
     +                  block, unit )
         else
            call wrt3s_ga( ga_file_double_buf, reals_to_write, unit )
         end if
*
*  and update the counters.
*
         num_ints = num_ints - ints_to_write
         where    = where    + ints_to_write
*
         force_loop = .false.
*
      end do
*
*  and reset the file safety mode
*
      call set_file_safety( old_safety, unit )
*
      end
_ENDIF
*
***********************************************************************
      subroutine clredx_ga
*
*  function to simulate clredx using ga tools. fairly 
* straight forward !
*
      end
*
**********************************************************************
*
      subroutine search_ga( block, unit )
*
*  subroutine to simulate search by using ga tools instead of disk.
*
*  this is only here for compatibility reasons. it is strongly
* recommened that the higher level routines ( wrt3 etc. ) are used.
*
      implicit none
*
      integer block
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/errcodes)
*
*  functions
*
      logical  get_oswed3
      external get_oswed3
*
*  local variables:
*
      integer ga_to_use
*
*  find out which ga is mapped to this unit number
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  if the file ain't open open it.
*
      if( .not. ga_file_open( ga_to_use ) ) then
         call rdtpnm_ga( unit )
      end if
*
*  check that the filing system has been initialized. it should be but
* you never know.
*
      if( .not. ga_file_initialized ) then
         print*, 'search: ga not initalized - internal error.'

         call gamerr('ga internal error',
     &        ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

      end if
*
*  if the main code changed the file io mode update the
* local mode maintainer.
*
      ga_file_mode( ga_to_use ) = get_oswed3( unit )
*
*  and here's the search
*
      ga_file_next( ga_to_use ) = block 
      call set_ipos( ga_file_next( ga_to_use ), unit )
*
*  and save the unit number
*
      ga_file_last_unit = unit
*
      end
*
***********************************************************************
*
      subroutine put_ga( array, words, unit )
*
* subroutine to simulate put by using ga tools instead of disk.
*
*  this is only here for compatibility reasons. it is strongly
* recommened that the higher level routines ( wrt3 etc. ) are used.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
      integer          words
      double precision array( 1:words )
      integer          unit
*
      call wrt3s_ga( array, words, unit )
*
      end
*
***********************************************************************
*
      subroutine get_ga( array, words )
*
* subroutine to simulate get by using ga tools instead of disk.
*
*  this is only here for compatibility reasons. it is strongly
* recommened that the higher level routines ( rdedx etc. ) are used.
*
*  in particular this routine is potentially dodgy. well no, not
* potentially, just plain is. this is because the ga files do not 
* store the number of records in a given block, and so this
* routine can not find how much has to be copied back into the
* array. hence it has to copy back the whole lot, leading to
* potential corruption of data after the end of the array.
*
*  there is one exception. restart runs go through secini which
* checks that the first block of the file is the right size. this
* one case, i.e. a get from block 1, is correctly handled.
*
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
      integer          words
      double precision array( 1:* )
*
*  local variables
*
      double precision tmp( 1:ga_file_block )
*
*  local variables
*
      integer          unit
*
*  pick up which unit was last searched.
*
      unit = ga_file_last_unit
  
ccc      write(6,*)'ga_get unit and block ',unit, 
ccc     &   ga_file_next(  ga_file_unit_to_ga( unit )  )

*
*  read the data and copy back.
*
      call reads_ga( tmp, ga_file_block, unit )
*
*  find out how many words were in block one on the disk.
*
      words = ga_file_leng1( ga_file_unit_to_ga( unit ) )

ccc      write(6,*)'ga_get words',words

*
*  and put the data back.
*
      call dcopy( words, tmp, 1, array, 1 )
*
      end
*     
***********************************************************************
*     
      subroutine delfil_ga( unit )
*     
*  subroutine to simulate delfil by using ga tools instead of disk.
*     
      implicit none
*     
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/discc)
INCLUDE(common/iofile)
*     
      integer unit
*     
* functions
*     
      logical  get_oswed3
      external get_oswed3

      logical pg_destroy_inf
      external pg_destroy_inf
*     
* local variables
*     
      integer ga_to_use
      integer last_block
      integer i
      logical mode
      logical first
      logical otmp
*
      save    first
*
      data first / .true. /

*
*  find out which ga file this is
*     
      ga_to_use = ga_file_unit_to_ga( unit )
*     
* if the file ain't open do nothing. this may seem a bit on the dodgy
* side but this is what the disk based delfil does.
*     
* also, if the ga file has not been written to, do nothing. needed when 
* an i/o error loading a ga from disk causes control to be passed
* to shut/delfil/delfil_ga
*
*
* Note for charmm/Chemshell work we need to close file properly
* even if it was not used yet
*
_IF(charmm,chemshell)
      if( ga_file_open( ga_to_use ) ) then
_ELSE
      if( ga_file_open( ga_to_use ) .and. 
     &    ga_file_hist_next( ga_to_use ) .gt. 1) then
_ENDIF
*     
*  check that all processors will go through this code. i.e.
* a ga file close is sychronizing.
*     
         mode = get_oswed3( unit )
*     
*  do the 'file' close and potential dump if oswed3 is hunky dory.
*     
         if( mode .and. COLLECTIVE ) then
*
*  now deal with the high water mark output. only the root node
* does this.
*

            if( ga_file_root ) then
*
*  if the first call write out a little header
*
               if( first ) then
                  write( iwr, * )
                  write( iwr, * ) ' gafs high water mark summary:'
                  write( iwr, * )
                  write( iwr, '( 1x, ''unit'', 5x, ''max blocks'' )' )
                  first = .false.
               end if
*
*  find the last block written to.
*
               last_block = -1
               do i = 1, ga_file_hist_next( ga_to_use ) - 1
                  last_block = max( last_block,
     +                         ga_file_history( 3, i, ga_to_use ) )
               end do
*
               write( iwr, '( 1x, i2, 6x, i7 )' ) unit - 1, last_block
*     
            end if
*     
*  ensure all transactions are complete. this is what replaces
* the call to waitedn in the disk based system.
*     
            call pg_synch( 160467 )
*     
*  if this file needs to be kept dump it to disk.
*     
            if( ga_file_keep( ga_to_use ) ) then
               call dump_ga_file( unit, zedfil(unit) )
            end if
*     
*  destroy the ga file
*     
            otmp = pg_destroy_inf( ga_file_handle( ga_to_use ),
     &                             yed(unit),'ga.m','delfil_ga' )
*     
*  restore the memory used by this file to the pool
* available for other files.
*     
            ga_file_total  = ga_file_total  - ga_file_size( ga_to_use )
            ga_file_unused = ga_file_unused + ga_file_size( ga_to_use )
*     
*  set the ga file descriptor back to its closed state values.
*     
            ga_file_handle   ( ga_to_use ) = not_in_use
            ga_file_open     ( ga_to_use ) = .false.
            ga_file_next     ( ga_to_use ) = not_in_use
            ga_file_hist_next( ga_to_use ) = 1
            ga_file_mode     ( ga_to_use ) = .true.
*     
*  need to set a couple of gamess's internal values.
*     
            call set_ipos( 0, unit )
            if( ga_file_keep( ga_to_use ) ) then
               call set_istat( 0, unit )
            else
               call set_istat( 2, unit )
            end if
*
*  and ensure the close is complete on all nodes.
*
            call pg_synch( 160467 )
*     
         else
*     
*  attempt to close a file while not doing collective
* io. flag a warning. can't flag an error because
* of deadlock.
*     
            print*, 'delfil: attempt to close a ga file while ' //
     +              'doing single node io'
            print*, 'request ignored'
         end if
*     
      end if
*
      end
*
***********************************************************************
*
      subroutine set_file_safety( safety_mode, unit )
*
*  subroutine to look after the safety mode that the file system is 
* running under.
*
      implicit none
*
      integer safety_mode
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  functions
*
      logical get_oswed3
*
*  local variables
*
      logical io_mode
*
*  find out how the file is currently being accessed
*
      io_mode = get_oswed3( unit )
*
      if( .not. ( io_mode .and. COLLECTIVE ) ) then
*
*  if only one process is reading from/writing to the file the is
* not a lot we can do except set the new safety mode.
*
         ga_file_safety_mode = safety_mode
*
      else
*
*  in this case all files are participating. if we are changing to the 
* safe mode do a synchronization first to ensure everything is cleared
* up
*
         if( safety_mode .eq. ga_file_safe ) then
            call pg_synch( 160467 )
         end if
*
         ga_file_safety_mode = safety_mode
*
      end if
*
      end
*
***********************************************************************
*
      subroutine ga_file_summary( unit, where_write )
*
*  noddy routine to write a summary of the history of a ga file. 
*
      implicit none
*
      integer unit
      integer where_write
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  functions
*
      logical  opg_root
      external opg_root
*
*  local variables
*
      integer ga_to_use
      integer i
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
      if( ga_file_open( ga_to_use ) ) then
*
         if( opg_root() ) then
*
            write( where_write, '( ''the following writes '',
     +                   ''occured to unit:'', i2, '':'' )' ) unit
*
            write( where_write, * )
*
            do i = 1, ga_file_hist_next( ga_to_use ) - 1
*
               write( where_write, 
     +             '( 1x, i7, 1x, i7, '' words to blocks '', i7,
     +                                '' to '', i7 )' ) i,
     +              ga_file_history( 1, i, ga_to_use ),
     +              ga_file_history( 2, i, ga_to_use ),
     +              ga_file_history( 3, i, ga_to_use )
*
            end do
*
         end if
*
      end if
*
      end 
*
***********************************************************************
*
      subroutine ga_file_history_compress( unit )
*
*  subroutine to remove duplicate entries from the history of
* writes to this ga file.
*
      implicit none
*
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
*  functions
*
      logical  opg_root
      external opg_root
*
*  local variables
*
      integer ga_to_use
      integer last_entry, next_entry
      integer words, start, finish
      integer j_start, j_finish
      integer i, j
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
      last_entry = ga_file_hist_next( ga_to_use ) - 1
*
*  loop over the entries backwards marking earlier entries
* that overlap with the current one under consideration.
*
      do i = last_entry, 1, -1
*
*  if words is -1 it means that this entry has been overwritten,
* so we can skip it.
*
         words  = ga_file_history( 1, i, ga_to_use )
         if( words .ne. -1 ) then
*
            start  = ga_file_history( 2, i, ga_to_use )
            finish = ga_file_history( 3, i, ga_to_use )
*
*  look for overlaps.
*
            do j = i - 1, 1, -1
*
               j_start  = ga_file_history( 2, j, ga_to_use )
               j_finish = ga_file_history( 3, j, ga_to_use )
*
               if( ( j_start  .ge. start    .and. 
     +               j_start  .le. finish ) .or.
     +             ( j_finish .ge. start    .and.
     +               j_finish .le. finish ) ) then
                  ga_file_history( 1, j, ga_to_use ) = -1
               end if
*
            end do
*
         end if
*
      end do
*
*  now go through the history cutting out the flagged entries.
*
      next_entry = 1
*
      do i = 1, last_entry
*
         if( ga_file_history( 1, i, ga_to_use ) .ne. -1 ) then
*
            ga_file_history( 1, next_entry, ga_to_use ) =
     +           ga_file_history( 1, i, ga_to_use ) 
            ga_file_history( 2, next_entry, ga_to_use ) =
     +           ga_file_history( 2, i, ga_to_use ) 
            ga_file_history( 3, next_entry, ga_to_use ) =
     +           ga_file_history( 3, i, ga_to_use ) 
*
            next_entry = next_entry + 1
*
         end if
*
      end do
*
*  and finally update the length indicator
*
      ga_file_hist_next( ga_to_use ) = next_entry
*
      end
*      
***********************************************************************
*     
      subroutine ga_file_history_sort( n, history )
*     
*  sort ga history file by heapsort. this is simply a slight tidy up of
* the numerical recipies routine. sorry there are no comments, i don't
* really understand the algorithm. it was chosen as it needs no
* workspace.
*
      implicit none
*
      integer n
      integer history( 1:3, 1:n )
*
*  local variables
*
      integer temp1, temp2, temp3
      integer p, q, i, j
*     
      p = n / 2 + 1
      q = n
*     
      do while( q .ne. 1 .or. p .ne. 1 )
*     
         if( p .gt. 1 ) then
*     
            p     = p - 1
            temp1 = history( 1, p )
            temp2 = history( 2, p )
            temp3 = history( 3, p )
*     
         else
*     
            temp1 = history( 1, q )
            temp2 = history( 2, q )
            temp3 = history( 3, q )
*     
            history( 1, q ) = history( 1, 1 )
            history( 2, q ) = history( 2, 1 )
            history( 3, q ) = history( 3, 1 )
*     
            q = q - 1
*     
         end if
*     
         if( q .eq. 1 .and. p .eq. 1 ) then
            history( 1, 1 ) = temp1
            history( 2, 1 ) = temp2
            history( 3, 1 ) = temp3
*     
         else
*     
            i = p
            j = p + p
*     
            do while( j .le. q )
               if( j .lt. q ) then
                  if( history( 2, j     ) .lt.
     +                history( 2, j + 1 ) ) then
                     j = j + 1
                  end if
               end if
               if( temp2 .lt. history( 2, j ) ) then
                  history( 1, i ) = history( 1, j )
                  history( 2, i ) = history( 2, j )
                  history( 3, i ) = history( 3, j )
                  i = j
                  j = j + j
               else
                  j = q + 1
               end if
            end do
*     
            history( 1, i ) = temp1
            history( 2, i ) = temp2
            history( 3, i ) = temp3
*     
         end if
*     
      end do
*     
      end
*
*********************************************************************
*
      subroutine dump_ga_file( unit, file_name )
*
*  subroutine to dump a ga file to disk in gamess format. to help
* reasonable transfer rates, both out of the ga file and also to
* the disk the internal buffer that is also used for data format
* conversion is used.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/timeperiods)
*
      integer unit
      character file_name*(*)
*
*  functions
*
      integer  lenwrd
      external lenwrd
*
*  local variables
*
      integer sizes( 1:blocks_wanted )
      integer ga_to_use
      integer error
      integer block
      integer blocks_write
      integer length
      integer word_size
      logical i_write

      call start_time_period(TP_IO_GAFILE_DUMP)

*
*  find out which ga file this unit is mapped to.
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  find out if this process will actually access the disk. the root
*  node will always do, the other depending on the dump mode.
*
      i_write = ga_file_root .or. 
     +          ga_file_dump_mode( ga_to_use ) .eq. 6
*
*  if appropriate open the file. 
*
      if( i_write ) then
         call strtrm( file_name, length )
         call opencc( file_name, length, unit, error )
      end if
*
*  sort out the history. compress removes duplicate entries in it
* while sort ermmm sorts it.
*
            call ga_file_history_compress( unit )
            call ga_file_history_sort( 
     +           ga_file_hist_next( ga_to_use ) - 1,
     +           ga_file_history( 1, 1, ga_to_use ) )
*
*  i personally would like to make this an option driven
* by the input file - it can really help in debugging !
*
*           call ga_file_summary( unit, 98 )
*
*  only process 0 is sure to have a complete version of the history,
* send it to the other nodes. this is required so that the all
* node participating version of rdedx_ga will, 'cos then all nodes
* will make consistent requests.
*
*
*  have to know how long an integer is to send these message.
*
      word_size = lenwrd()
*
* tell all other nodes how big the root node's history is for this 
* file.
*
      call pg_brdcst( 160467, ga_file_hist_next( ga_to_use ), 
     +                8 / word_size, 0 )
*
*  now send out the history.
*      
      length = 3 * ga_file_hist_next( ga_to_use ) * 8 / word_size
      call pg_brdcst( 160467, ga_file_history( 1, 1, ga_to_use ), 
     +                length, 0 )
*
*  block is the first block that we want to try to write to file.
*
      block = 1
*
*  look at the history and try and find upto blocks_wanted contiguous
* blocks that contain data. blocks_write returns how many blocks
* we actually found. for these blocks in sizes reurn how
* much data is in each of the blocks.
*
      call ga_file_history_interpret( block, blocks_wanted, 
     +                                blocks_write, sizes, ga_to_use )
*
*  while we have some blocks ....
*
      do while( blocks_write .ne. 0 )
*
*  read the blocks out out of the ga file.
*
         call rdedx_ga( ga_file_double_buf, 
     +                  blocks_write * ga_file_block,
     +                  block, unit )
*
*  change the data to the format required by the disk i.e. at the
* end of each block insert the amount of data present information.
*
         call reformat_ga_to_disk( block, blocks_write, sizes,
     +                             ga_file_double_buf, 
     +                             ga_file_double_buf )
*
*  write out the blocks.
*
         if( i_write ) then
            call write_blocks( block, blocks_write, 
     +                         ga_file_double_buf, unit )
         end if
*
*  update the block counter.
*
         block = block + blocks_write
*
*  reinvestigate the history looking for so more blocks to write.
*
         call ga_file_history_interpret( block, blocks_wanted, 
     +                                   blocks_write, sizes, 
     +                                   ga_to_use )
*
      end do
*

      if( i_write ) then
          call closecc( unit )
      endif

      call end_time_period(TP_IO_GAFILE_DUMP)

      end
*
***********************************************************************
*
      subroutine ga_file_set_up( n, dump_mode, length, keep, name)
*
*  subroutine that initializes the ga files and sets up the 
* characteristics of each file. n.b. it does not actually open
* the files.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/errcodes)
*
      integer      n
      integer      dump_mode( 1:n )
      integer      length   ( 1:n )
      integer      keep     ( 1:n )
      character*(*) name     ( 1:n )
*
*  functions
*
      integer  ipg_nnodes
      external ipg_nnodes
      integer  ipg_nodeid
      external ipg_nodeid
      logical  opg_root
      external opg_root
*
*  local variables
*
      integer total_gas
      integer default_size
      integer ga_to_use
      integer i
*
*  set up the ga file system.
*
      ga_file_initialized = .true.
      ga_file_processes   = ipg_nnodes()
      ga_file_proc_id     = ipg_nodeid()
      ga_file_root        = opg_root()
*
*  count how many ga files we are going to use. the burnt
* in 6 and 7 are inherited from gamess, would be better parameterized.
*
      total_gas = 0
      do i = 1, n
         if( dump_mode( i ) .eq. 6 .or. dump_mode( i ) .eq. 7 ) then
            total_gas = total_gas + 1
         end if
      end do
*
*  check that we can use this amount of ga files.
*
      if( total_gas .gt. max_ga_files ) then
         print*, 'too many ga files requested'
         print*, 'maximum is ', max_ga_files
         print*, 'requested  ', total_gas

         call gamerr('ga dimensioning problem',
     &        ERR_NO_CODE, err_dimension, ERR_ASYNC, ERR_NO_SYS)

      end if
*
      if( total_gas .ne. 0 ) then
*
*  we have some ga files - set up their characteristics.
* first of all the default amount of memory allocated to a file.
*
         default_size = ga_file_proc_mem * ga_file_processes / 
     +                  total_gas
*
*  now set up the characteristics.
*
         ga_to_use = 0
         do i = 1, n
*
            if( dump_mode( i ) .eq. 6 .or. 
     +          dump_mode( i ) .eq. 7 ) then
*
               ga_to_use = ga_to_use + 1
*
*  set up the unit to ga mapping array
*
               ga_file_unit_to_ga( i ) = ga_to_use
*
*  a name, a keep status and a dump mode are always provided.
*
               ga_file_keep     ( ga_to_use ) = keep( i ) .ne. 0
               ga_file_dump_mode( ga_to_use ) = dump_mode( i )
*
*  on the other hand a length need not be. use the default if
* it ain't. also non-default values are in blocks.
*
               if( length( i ) .gt. 0 ) then
                  ga_file_size( ga_to_use ) = length( i ) * 
     +                                        ga_file_block
               else
                  ga_file_size( ga_to_use ) = default_size
               end if
*
            end if
*
         end do
*
*  now work out how much memory we are going to use
*
         ga_file_total_mem = 0
         do i = 1, total_gas
            ga_file_total_mem = ga_file_total_mem + ga_file_size( i )
         end do
*
*  and finally set up the counters on this memory
*
         ga_file_total       = 0
         ga_file_unused      = ga_file_total_mem
*
      end if
*
      end
*
***********************************************************************
*
      logical function get_oswed3( unit )
*
*  function that returns the current mode of access to a file
* as controlled by the oswed3 variable.
*
      implicit none
*
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/nodeio)
INCLUDE(common/parcntl)
*
      if (unit.eq.20) then
         get_oswed3 = oswed3( unit )
      else
         get_oswed3 = oswed3( unit ) .or. ( ipiomode .eq. IO_NZ) 
      endif
*
      end
*
***********************************************************************
*
      subroutine set_ipos( position, unit )
*
*  subroutine to set gamess's internal file position value.
*
      implicit none
*
      integer position
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/disc)
*
      ipos( unit ) = position
*
      end
*
***********************************************************************
*
      subroutine set_istat( status, unit )
*
*  subroutine to set gamess's internal file status value.
*
      implicit none
*
      integer status
      integer unit
*
INCLUDE(common/sizes)
INCLUDE(common/disc)
*
      istat( unit ) = status
*
      end
*
***********************************************************************
*
      subroutine reformat_ga_to_disk( base_block, blocks, sizes, 
     +                                in_array, out_array )
*
*  subroutine that takes an array in the format as stored in
* ga and blocks it in the form suitable to write to the disk.
* we have to be a little careful because in_array and out_array 
* may be the same.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*     
      integer          base_block
      integer          blocks
      integer          sizes( 1:blocks )
      double precision in_array ( 1:blocks * ga_file_block )
      double precision out_array( 1:blocks * ( ga_file_block + 1 ) )
*
*
*  functions
*
_IFN(bits8)
      double precision pack2
      external         pack2
_ENDIF
*
*  local variables
*
      integer this_block
      integer unpack
      integer pack
      integer i, j
*
      this_block = base_block + blocks - 1
*
      pack   = ( blocks - 1 ) * ( ga_file_block + 1 ) +
     +         ga_file_block + 1
*
*  pack into the output array going backwards to avoid corrupting 
* data if the arrays are the same.
*
      do i = blocks, 1, -1
*
_IFN(bits8)
         out_array( pack ) = pack2( this_block, sizes( i ) )
_ELSE
         call pack2( this_block, sizes( i ), out_array( pack ) )
_ENDIF
         this_block        = this_block - 1
*
         unpack = ( i - 1 ) * ga_file_block         + sizes( i )
         pack   = ( i - 1 ) * ( ga_file_block + 1 ) + sizes( i )
*
         do j = sizes( i ), 1, -1
*
            out_array( pack ) = in_array( unpack )
*
            unpack = unpack - 1
            pack   = pack   - 1
*
         end do
*
      end do
*
      end
*
***********************************************************************
*
      subroutine ga_file_history_interpret( base_block, 
     +                                      blocks_required,
     +                                      blocks_got, sizes,
     +                                      ga_to_use )
*
*  subroutine that looks at the history an returns the sizes of the
* blocks requested. basically you ask it what are the sizes of
* the next blocks_required from base_block. if this takes you
* over the end of the history blocks_got will be less than
* blocks_required
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
      integer base_block
      integer blocks_required
      integer blocks_got
      integer sizes( 1:blocks_required )
      integer ga_to_use
*
*  local variables
*
      integer this_block
      integer history_position
      integer words_from_entry
      integer entry_start
      integer entry_end
      integer last_start
*
*  find out where in the history this block is.
*
      call locate_block( base_block, ga_to_use, history_position )
*
      blocks_got = 0
*
      this_block = base_block
      last_start = this_block
*
*  while we have not got the requesite blocks and there are still
* blocks to get investigate the history.
*
      do while( blocks_got .lt. blocks_required .and.
     +          history_position .lt. ga_file_hist_next( ga_to_use ) )
*
         blocks_got = blocks_got + 1
*
*  find out the actual details of this history entry.
*
         words_from_entry = ga_file_history( 1, history_position, 
     +                                       ga_to_use )
         entry_start      = ga_file_history( 2, history_position, 
     +                                       ga_to_use )
         entry_end        = ga_file_history( 3, history_position, 
     +                                       ga_to_use )
*
*  if this entry is for a latter block it means that there must
* be some blank blocks. skip them one at a time. this could probably
* be improved.
*
         if( entry_start .gt. this_block ) then
            sizes( blocks_got ) = 0
*
         else
*
*  o.k. this history entry covers the block of interest,
* but the block of interest may be right in the middle of
* a transfer. e.g. we may be looking at block 41, but the
* history entry may cover blocks 39 to 45. this if condition
* accounts for this.
*
            if( entry_start .lt. this_block ) then
               words_from_entry = words_from_entry - ga_file_block *
     +                ( this_block - entry_start )
*
*  if words_from_entry is less than zero then there is a gap
* in the dump file.
*
               words_from_entry = max( words_from_entry, 0 )
*
            end if
*
*  work out the size of this block.
*
            sizes( blocks_got ) = min( ga_file_block, 
     +                                 words_from_entry )
*
*  if this corresponds to the last block in this history
* entry mover onto the next entry.
*
            if( this_block .eq. entry_end ) then
               history_position = history_position + 1
            end if
*
         end if
*
         this_block = this_block + 1
*
      end do
*
      end
*
**********************************************************************
*
      subroutine locate_block( block, ga_to_use, history_position )
*
*  subroutine to locate the history entry corresponding to this block.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
      integer block
      integer ga_to_use
      integer history_position
*
*  local variables
*
      integer last_history
      integer last_block
      integer wanted
      integer lower, middle, upper
*
      last_history = ga_file_hist_next( ga_to_use ) - 1
      last_block   = ga_file_history( 3, last_history, ga_to_use )
*
*  first check that we are not over the end of the history. if
* so we can return immediately.
*
      if( block .gt. last_block ) then
         history_position = last_history + 1
      else
*
*  otherwise use a binary search to locate the data. remember that
* the history has been sorted by this point.
*
         upper = last_history + 1
         lower = 0
*
         do while( upper - lower .gt. 1 )
            middle = ( upper + lower ) / 2
            if( ga_file_history( 2, middle, ga_to_use ) .gt.
     +          block ) then
               upper = middle
            else
               lower = middle
            end if
         end do
*
         history_position = lower
*
      end if
*
      end
*
***********************************************************************
*
      subroutine write_blocks( base_block, blocks, array, unit )
*
*  dump blocks blocks of array to disk starting at base_block
* in unit unit.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/errcodes)

*
      integer          base_block
      integer          blocks
      double precision array( 1:blocks * ( ga_file_block + 1 ) )
      integer          unit
*
*  local variables
*
      integer error
*
      call srchcc( unit, base_block, error )
      if( error .ne. 0 ) then

         call gamerr('i/o error (search) when writing ga to disk',
     &        ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_SYS)

      end if
*
      call putccn( unit, array, blocks, error )
      if( error .ne. 0 ) then

         call gamerr('i/o error (put) when writing ga to disk',
     &        ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_SYS)

      end if
*
      end
*
***********************************************************************
*
      subroutine read_restart( unit, file_name )
*
*  subroutine to read a disk file into the gas. to help
* reasonable transfer rates, both into the ga file and also from
* the disk the internal buffer that is also used for data format
* conversion is used.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/timeperiods)
INCLUDE(common/errcodes)
*
      character file_name*(*)
      integer unit
*
*  functions
*
      integer  lenwrd
      external lenwrd
*
*  local variables
*
      integer sizes( 1:blocks_wanted )
      integer ga_to_use
      integer file_length
      integer error
      integer dummy
      integer length
      integer this_block
      integer blocks_read
      integer next_block
      integer i
      logical exist
*
      call start_time_period(TP_IO_GAFILE_READ)
*
*  map the unit number to the appropriate ga
*
      ga_to_use = ga_file_unit_to_ga( unit )
*
*  only node 0 will ever do any disk accesses.
*
      if( ga_file_root ) then
*
*  first things first. does the file exist ?
* this is bodged - what if in indi iomode.
*
         inquire( file = file_name,  exist = exist )
         if( .not. exist ) then

            write(6,*)'problem reading from file',file_name
            call gamerr(
     &   'ga file error - attempt to read from non-existent file',
     &           ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_NO_SYS)

         end if
*
*  find out how big the file is in blocks.
*
         call length_file( file_name, file_length,
     +                     error )
         if( error .ne. 0 ) then

            call gamerr('read_restart: error finding file length',
     &           ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

         end if
*
*  trap empty files.
*
         if( file_length .eq. 0 ) then

            call gamerr('read_restart: empty file ',
     &           ERR_NO_CODE, ERR_INTERNAL, ERR_ASYNC, ERR_NO_SYS)

         end if
*
*  open the file. 
*
         call strtrm( file_name, length )
         call opencc( file_name, length, unit, error )
*
*  to fool secini in gamess we need to look at
* the first block on the disk, find out how long it is and
* store it. this is used in get_ga.
*
         call read_blocks( 1, 1, ga_file_double_buf, unit )
         call upack2( ga_file_double_buf( 512 ), dummy,
     +                ga_file_leng1( ga_to_use ) )
*
      end if
*
*  make the length of the file general knowledge.
*
      call pg_brdcst( 160467, file_length, 8 / lenwrd(), 0 )
*
*  make the length of block 1 general knowledge.
*
      call pg_brdcst( 160467, ga_file_leng1( ga_to_use ), 
     +                8 / lenwrd(), 0 )
*
*  now just bring in the file of the disk in as big a chunk as 
* possible, reformat the chunk and put it in the ga. note that
* as this is the first time the ga is accessed the file is already
* `rewound'
*
      this_block = 1
*
      do while( this_block .lt. file_length )
*
         blocks_read = min( blocks_wanted, 
     +                      file_length - this_block + 1 )
*
         if( ga_file_root ) then
            call read_blocks( this_block, blocks_read, 
     +                        ga_file_double_buf, unit )
         end if
*
*  make the data read in general knowledge 'cos the ga write will
* be done in collective io mode.
*
         call pg_brdcst( 160467, ga_file_double_buf, 
     +                   blocks_read * ( ga_file_block + 1 ) * 8, 0 )
*
         call reformat_disk_to_ga( blocks_read, sizes, 
     +                             ga_file_double_buf,
     +                             ga_file_double_buf )
*
*  have to be careful to ensure that the history of the writes is
* correct so that when the file gets dumped at the end of the
* job it is in the correct format.
*
         do i = 1, blocks_read
            next_block = ( i - 1 ) * ga_file_block + 1
            call wrt3s_ga( ga_file_double_buf( next_block ),
     +                     sizes( i ), unit )
         end do
*
         this_block = this_block + blocks_read
*
      end do
*
*  close the file so that the unit can be reused latter.
*
      call closecc( unit )
*
*  what if the file needs deleting at the end of the job ?
*
*
*  and rewind the ga file.
*
      call search_ga( 1, unit )
*
      call end_time_period(TP_IO_GAFILE_READ)

      end
*
***********************************************************************
*
      subroutine read_blocks( base_block, blocks, array, unit )
*
*  read blocks blocks into array from disk starting at base_block
* in unit unit.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
INCLUDE(common/errcodes)
*
      integer          base_block
      integer          blocks
      double precision array( 1:blocks * ( ga_file_block + 1 ) )
      integer          unit
*
*  local variables
*
      integer error
*
      call srchcc( unit, base_block, error )
      if( error .ne. 0 ) then

         call gamerr('i/o error (search) when reading ga from disk',
     &        ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_SYS)

      end if
*
      call getccn( unit, array, blocks, error )
      if( error .ne. 0 ) then

         call gamerr('i/o error (get) when reading ga from disk',
     &        ERR_NO_CODE, ERR_NO_CLASS, ERR_ASYNC, ERR_SYS)

      end if
*
      end
*
***********************************************************************
*
      subroutine reformat_disk_to_ga( blocks, sizes,
     +                                in_array, out_array )
*
*  subroutine to change an array read infrom the disk to the form
* suitable for the ga tools. sizes returns the size of the blocks.
* have to be a little careful because in_array and out_array may
* be the same.
*
      implicit none
*
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
*
      integer          blocks
      integer          sizes( 1:blocks )
      double precision in_array ( 1:blocks * ( ga_file_block + 1 ) )
      double precision out_array( 1:blocks * ga_file_block )
*
*  local variables
*
      integer disk, ga
      integer dummy
      integer i, j
*
      disk = 1
      ga   = 1
*
      do i = 1, blocks
         do j = 1, ga_file_block
            out_array( ga ) = in_array( disk )
            disk = disk + 1
            ga   = ga + 1
         end do
         call upack2( in_array( disk ), dummy, sizes( i ) )
         disk = disk + 1
      end do
*
      end
*
***********************************************************************

      block data block_anal_storage
INCLUDE(common/gaanal)
      data declared /.false./
      data supressed /.false./
      end

      subroutine declare_anal_storage(ndim)
      implicit none
      integer ndim, igmem_alloc

INCLUDE(common/gaanal)
INCLUDE(common/parcntl)
INCLUDE(common/vcore)
INCLUDE(common/iofile)

      logical pg_create

      integer i, ltri

_IF(charmm,chemshell)
*
* Ensure loading of block data (cray requirement)
*
      external block_anal_storage
_ENDIF
      if (declared) return

      if(idpdiis .eq. 99999999 .and.
     &   idporth .eq. 99999999 .and.
     &   idpmult2 .eq. 99999999)then
c
         write(iwr,*)'ANAL GA allocation supressed'

         supressed = .true.

      endif	

c
c one square buffer for scratch sum/mult2
c
      if (.not. pg_create(0,ndim,ndim,'sqbuf',0,0,
     &     ih_scr)) then
         call pg_error('failed to create buffer GA ',i)
      endif
c
c square buffer for scratch sum/mult2
c
      if (.not. pg_create(0,ndim,ndim,'sqbuf2',0,0,
     &     ih_scr2)) then
         call pg_error('failed to create buffer GA ',i)
      endif
c
c one square buffer for vectors
c
      if (.not. pg_create(0,ndim,ndim,'vec',0,0,
     &     ih_vec)) then
         call pg_error('failed to create vec GA ',i)
      endif
c
c buffer for overlap matrix  (store now)
c
      ltri = (ndim+1)*ndim/2

c
      if (.not. pg_create(0,ndim,ndim,'ov',0,0,
     &     ih_ov)) then
         call pg_error('failed to create overlap GA ',i)
      endif

      declared = .true.

      end

      subroutine destroy_anal_storage
INCLUDE(common/gaanal)
      logical pg_destroy,o
      if (supressed) return
      o = .true.
      o = o .and. pg_destroy(ih_ov)
      o = o .and. pg_destroy(ih_scr)
      o = o .and. pg_destroy(ih_scr2)
      o = o .and. pg_destroy(ih_vec)
      if(.not. o)call caserr('error destrying ANAL GAs')
      declared = .false.
      end

      subroutine get_ed19adres(istart,iend,il,ih,jl,jh)
      implicit REAL (a-h,o-z)
INCLUDE(common/sizes)
INCLUDE(common/ga_file)
      common/remc_ed19/iblok
      if (iblok.eq.0) call caserr('error 1 in get_ed19adres')
      ik=iblok*ga_file_block
c default c-vector
      il=istart
      ih=iend
      jl=1
      jh=1
      if (il.gt.ik) then
c z-vector ...
         il=il-ik
         jl=2
         jh=2
         if (ih.lt.ik)call caserr('what do you want c or z??')
         ih=ih-ik
      endif
      if (il.gt.ik)call caserr('internal c/z-file error')
      return
      end
_IFN(scalapack)
c
c     GA tools were build without the ScaLAPACK interfaces.
c     Therefore we need to trap attempts to use the missing routines.
c
      subroutine ga_pdsyev(g_a,g_q,e,nb)
      implicit none
      integer g_a, g_q, nb
      REAL e(*)
      call caserr("ga_pdsyev: no ScaLAPACK interfaces in this build")
      end
c
      subroutine ga_pdsyevx(g_a,g_q,e,nb)
      implicit none
      integer g_a, g_q, nb
      REAL e(*)
      call caserr("ga_pdsyevx: no ScaLAPACK interfaces in this build")
      end
c
      subroutine ga_pdsyevd(g_a,g_q,e,nb)
      implicit none
      integer g_a, g_q, nb
      REAL e(*)
      call caserr("ga_pdsyevd: no ScaLAPACK interfaces in this build")
      end
c
      subroutine ga_pdsyevr(g_a,g_q,e,nb)
      implicit none
      integer g_a, g_q, nb
      REAL e(*)
      call caserr("ga_pdsyevr: no ScaLAPACK interfaces in this build")
      end
c
_ENDIF
_ELSE
c 
c ga tools were not selected when the code was generated
c
_ENDIF(ga)
      subroutine ver_ga(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/ga.m,v $
     +     "/
      data revision /"$Revision: 6115 $"/
      data date /"$Date: 2009-12-17 11:37:45 +0100 (Thu, 17 Dec 2009) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
