*********************************************************************
c
c     ZORA
c     Zeroth Order Regular Approximation 
c     Scalar Relativistic Correction
c     S. Faas 1997-2000
c
c     the following approximations/choices are available
c     Integral approximations:
c     - Full Coulomb (icoul_z=1) - no approximations
c     - Atomic (icoul_z=2) both 1- and 2-electron ints for ZORA
c              are 1-center. All other matrices complete.
c     - No Coulomb (icoul_z=3) Only 1-electron coulomb ints, 
c              both 1- and 2-center
c     - Small Coulomb (icoul_z=4) as full, but the 2-electron 
c              coulomb operator is generated in external basis.
c     - The exchange may be added to the coulomb to mimic the
c              DFT ZORA  (for icoul_z <=2)
c     Density matrix used to calculate zora corrections:
c     - Molecular - may be picked up (get) or selfconsistently
c                   calculated.
c     - Atomic  (use get atom / get atoi ..)
c       This implies 1-center 1- and 2-electron integrals
c       (icoul_z >1 still possible) and all transformations
c       and zora corrections to be 1-center
c     Unscaled, scaled and iora ZORA 
c       The iora not for Atomic density
c     Choices for internal basis 
c     - may be generated automatically or by hand
c     - restriction to s,p,d or f functions
c     - contraction of external basis may be retained
c     Print options : print (some) and pri2 (a lot)
c     The speed of light may be set (:))
c     
c***********************************************************************
c 
      subroutine trans_z(sd1,sx,s2,num)
c
      implicit REAL (a-h,p-z),integer (i-n),logical (o)
      dimension sd1(num,num),sx(num,num),s2(num,num)
c   
INCLUDE(common/zorac)
INCLUDE(common/iofile)
c
c all matrices triangular, except s2;  anti-symmetrize and transform
c to orthogonal basis
c
c transform (a' = s2 * a * s2^+ ) s2: sym.
c
      l2 = num*(num+1)/2
      l3 = num*num
c
c...  s**-1/2 (internal)
c
      call rdedx(s2,l3,ibsmin_z,num8)
c
c...   x
c
      call rdedx(sx,l2,ibsx_z,num8)
      call antisquare(sd1,sx,num)
      call mxmaa(s2,num,1,sd1,num,1,sx,num,1,num,num,num)
      call mxmaa(sx,num,1,s2,num,1,sd1,num,1,num,num,num)
      call wrt3(sd1,l3,ibsx_z,num8)
c
c...   y
c
      call rdedx(sx,l2,ibsy_z,num8)
      call antisquare(sd1,sx,num)
      call mxmaa(s2,num,1,sd1,num,1,sx,num,1,num,num,num)
      call mxmaa(sx,num,1,s2,num,1,sd1,num,1,num,num,num)
      call wrt3(sd1,l3,ibsy_z,num8)
c
c...   z
c
      call rdedx(sx,l2,ibsz_z,num8)
      call antisquare(sd1,sx,num)
      call mxmaa(s2,num,1,sd1,num,1,sx,num,1,num,num,num)
      call mxmaa(sx,num,1,s2,num,1,sd1,num,1,num,num,num)
      call wrt3(sd1,l3,ibsz_z,num8)
c
      return
      end
c
c****************************************************************
c
      subroutine scale_z(q,eng,engo,eigv,t2,ehf,etot,ne,num,imode)
c
c subroutine calculates scaling of orbital energies
c
c eng contains fock operator in basis of old orbitals (previous
c iteration). so diagonal with orbital energies.
c
c ... imode = 0 : perform original scaling on orbital energies after
c ...             scf convergence
c
c ... imode = 1 : perform scaling on MO's during SCF iterations
c ...             (development version) save scaling factors here
c ...             used (for scaling MO's) in rhfclm

      implicit REAL    (a-h,p-z),integer   (i-n),logical    (o)
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
c
INCLUDE(common/zorac)
INCLUDE(common/iofile)
INCLUDE(common/scfwfn)
c
      dimension eng(*),engo(*),eigv(num,*),t2(*),q(*)
c
      if (nopen+npair.eq.0) nco = ne/2
c...    do not know what to scale in GVB
      if (npair.ne.0) return
c
      write(iwr,605)
605   format(2x,30('='),' orbital energies scaled ',30('='))
c
      call rdedx(t2,nwcor_z,ibscale_z,num8)
c
c calculate total scaled energy 
c
      if (imode .eq. 0 .and. .not. oatscf_z) then
c
          sum=0.0d0
c         ne2 = ne/2
c only occupied orbitals; no open shells!
          do i=1,nco
             sum=sum+eng(iky(i+1))*2.0d0
          end do
          do i=nco+1,nco+nopen
             sum=sum+eng(iky(i+1))*1.0d0
          end do
c
          if (opzora) write(iwr,603) sum
          ehf = ehf - sum
          etot = etot - sum
c
      endif 
c
c...  allocate space for scale factors 
c
      iscal1 = igmem_alloc(num_ext)
c
      do k=1,num
         ikk=iky(k+1) 
         sum=0.0d0
         do i=1,num
            do j=1,num
               ij = min(i,j) + iky(max(i,j))
               sum = sum - eigv(j,k)*eigv(i,k)*t2(ij)
            enddo
         enddo 
         engo(ikk)  = eng(ikk)
         if (imode .ne. 1) eng(ikk) = eng(ikk)/(1.0d0+sum) 
         q(iscal1+k-1) = 1.0d0/(1.0d0+sum)
      enddo
c
c ... write scale factors to file sf
c
      call wrt3(q(iscal1),num_ext,ibscalf_z,num8)
c     
      call gmem_free(iscal1)

      if (imode .ne. 0) return
c
      if (.not. oatscf_z) then
         sum=0.0d0
         do i=1,nco
            sum=sum+eng(iky(i+1))*2.0d0
         enddo
         do i=nco+1,nco+nopen
            sum=sum+eng(iky(i+1))*1.0d0
         end do
         ehf = ehf + sum
         etot = etot + sum
         if (opzora) then
            write(iwr,602) sum
            write(iwr,601) ehf
            write(iwr,604) etot
         end if
      end if
c
c...  now compute the scaling corrections to the one-electron
c...  integrals and replace the unscaled corrections on disk by them
c
      do i=1,num
         engo(i) = -engo(iky(i+1)) + eng(iky(i+1))
      end do
c
      isd6 = igmem_alloc(4*num*num+1)
      isd7 = isd6 + num*num
      isd8 = isd7 + num*num 
      isd9 = isd8 +  num*num
c invert matrix of eigenvectors
      call dcopy(num*num,eigv,1,q(isd8),1)
      d = 0.0d0
      call minvrt(q(isd8),num,d,q(isd9),q(isd9+num))
c store scaling corrections in square matrix for transformation
      call vclr(q(isd6),1,2*num*num)
      do i=1,num
          q(isd6+(i-1)*num+i-1) = engo(i)
      enddo
c transform scaling corrections (MO-basis) to AO basis
      call mxmaa(q(isd6),1,num,q(isd8),1,num,q(isd6+num*num),
     + 1,num,num,num,num)
      call mxmaa(q(isd8),num,1,q(isd6+num*num),1,num,q(isd6),
     + 1,num,num,num,num)
c store triangular
      call triangle(q(isd6),q(isd6),num)
      call dcopy(nwcor_z,q(isd6),1,t2,1)
c
      call gmem_free(isd6)
c     saveascaling corrections in ao-basis
      call wrt3(t2,nwcor_z,ibscalao_z,num8)
c
601   format ('electronic scaled zora energy:',8x,f12.4)
602   format ('sum of scaled orbital energies:',2x,f12.4)
603   format ('sum of unscaled orbital energies:',f12.4)
604   format ('total scaled zora energy:',13x,f12.4)
c
      return
      end
c
      subroutine scale_uhf(q,eng,eigv,t2,descal,ne,
     +                     num,ne_z)
c
c subroutine calculates scaling of orbital energies
c
c eng contains fock operator in basis of old orbitals (previous
c iteration). so diagonal with orbital energies.
c
      implicit REAL    (a-h,p-z),integer   (i-n),logical    (o)
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
c
INCLUDE(common/zorac)
INCLUDE(common/iofile)
c
      dimension eng(*),eigv(*),t2(*),q(*)
c
      write(iwr,605)
605   format(2x,30('='),' orbital energies scaled ',30('='))
c
      l2 = num*(num+1)/2
      isd1 = igmem_alloc(2*num*num)
      isd2 = isd1 + num*num
c
      call rdedx(t2,nwcor_z,ibscale_z,num8)
      do i=1,l2
         q(isd2+i-1)=-1.0d0*t2(i)
      end do
      call mult2(eigv,q(isd1),q(isd2),num,num,num)
c
c     do i=1,num
c        t2(i) = q(isd1+iky(i))
c     end do
c
      call dcopy(l2,q(isd1),1,t2,1)
      do i=1,num
         t2(iky(i+1))=1.0d0*t2(iky(i+1))+1.0d0
      enddo
c
c calculate total scaled energy 
c
      sum1=0.0d0
c
c uhf scaling
c 
      do 10 i=1,ne_z
         sum1=sum1+eng(i)
10    continue
c
76    format (i2,f14.4)
c
      do i=1,num
         eng(i) = eng(i)/t2(iky(i+1))
      end do
      sum=0.0d0
      do i=1,ne_z
         sum=sum+eng(i)
      enddo
c difference between sum of scaled and unscaled energies
      descal = descal + sum1 -sum
      if (opzora) write(iwr,604)descal 
c
604   format ('diff between sum of scaled and unscaled:',13x,f12.4)
c
      call gmem_free(isd1)
      return
      end
c
c***********************************************************************
c
      subroutine pre_zora(sd)
c
      implicit REAL    (a-h,p-z),integer   (i-n),logical    (o)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/mapper)
INCLUDE(common/scra7)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/zorac)
c
      integer onormf1,onormp1
      common/nshel_z/expz(mxprim),csz(mxprim),cpz(mxprim),cdz(mxprim),
     +               cfz(mxprim),cgz(mxprim),
     +       kstartz(mxshel),katomz(mxshel),ktypez(mxshel),
     +       kngz(mxshel),klocz(mxshel),kminz(mxshel),kmaxz(mxshel),
     +               nshellz,nonz,numorbz,ndummz
     +      ,onormf1,onormp1,c2(3,maxat)
c
      logical ostf(6)
      dimension sd(*),array(10)
c
      data ostf/.true.,5*.false./
c
      o1e_zora = .false.
c compute one electron integrals again ab
c sf    compute derivative and transformation matrices needed for
c       zora calculation.
c
      non_z=non
      numorbz=numorb
c
c generate internal basis
c
      num2 = nbas_zora
c
      lneed = 3*num2*num2  + 2*num2 + num2*(num2+1)/2
c
      lavail = igmem_max_memory()
c
      if (lavail .lt. lneed) call caserr ('ran out of memory ')
c
      iscr = igmem_alloc(lneed)
      isq = iscr + num2*num2
      iout = isq + num2*num2
      itr = iout + num2*num2
      id1 = itr + num2*(num2+1)/2
c
c construct overlap between external and internal basis and calculate internal 1-el. ints
c
      call scrosz(sd,sd(iscr),nbas_zora)
c
c     sd(iscr_ = scros = overlap between external and internal basis
c
      if (op2zora) then
         write(iwr,*) 'overlap between internal and external basis'
         call prsq (sd(iscr),num,num2,num2)
      endif
c
c
      call getmat(sd(itr),
     1            sd(isq),sd(isq),sd(isq),sd(isq),sd(isq),array,
     2            num,ostf,ionsec)
c
c...  block diagonalise over atoms for atomic startup 
c... (s-matrix was complete)
c
      if (nat_z .gt. 0)  call block_z(sd(itr),1,'tri')
c
      d=0.0d0 
      call square(sd(isq),sd(itr),num,num)
c invert
      call minvrt(sd(isq),num,d,sd(id1),sd(id1+num)) 
c
c  matrix multiplication to generate projection external to internal
c  is s(external)**-1 * scros 
c
      call mxmaa(sd(isq),1,num,sd(iscr),num2,1,sd(iout),1,num,
     1           num,num,num2)
      call wrt3(sd(iout),nwtrat_z,ibtrin_z,num8)
c
c constructing s^-1 internal basis
c
      call rdedx(sd(itr),num2*(num2+1)/2,ibls_z,num8)
      if (op2zora) then 
         write(iwr,*) 'overlap in internal basis'
         call prtri(sd(itr),num2)
      endif
c
      call square(sd(isq),sd(itr),num2,num2)
c
      d=0.0d0
      call minvrt(sd(isq),num2,d,sd(id1),sd(id1+num2))
      call mxmaa(sd(isq),1,num2,sd(iscr),1,num2,sd(iout),1,num2,
     1           num2,num2,num)
      call wrt3(sd(iout),nwtrat_z,ibtrout_z,num8)
c
      if (op2zora) then
         write(iwr,*) 'resulting matrix 2:'
         call prsq(sd(iout),num,num2,num2)
      endif
c
c....  make s**0.5 and s**-0.5 in internal basis
c....  iscr is scratch now
c     
      call rdedx(sd(itr),num2*(num2+1)/2,ibls_z,num8)
c generate s**1/2
      call symtrn(sd(itr),sd(iscr),sd(isq),
     +            iky,sd(id1),num2,num2*(num2+1)/2,num2,' pre_zora ')
c
      if (op2zora) then
         write(iwr,*) ' S**0.5 '
         call prtri(sd(itr),num2)
      end if
c
      call square(sd(iout),sd(itr),num2,num2)
c
      call wrt3(sd(iout),nwmat_z,ibsplus_z,num8)
c
c inverting s**1/2
c
      d=0.0d0
      call minvrt(sd(iout),num2,d,sd(id1),sd(id1+num2))
c
      call wrt3(sd(iout),nwmat_z,ibsmin_z,num8)
c
      if (op2zora) then
         write(iwr,*) 's2 new',d
         call prsq (sd(iout),num2,num2,num2)
      endif
c
c now sx matrices, transform and antisym; minus sign!
c transform and write  derivative matrices to orthogonal basis
c
      call trans_z(sd(isq),sd(iscr),sd(iout),num2)
c  
      call gmem_free(iscr)
c
      return
      end
c
c
c***********************************************************************
c
c
      subroutine scrosz(q,sout,mdim)
c
c...   call scros after taking care that the external basis is in junk
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension q(*),sout(mdim,mdim)
INCLUDE(common/sizes)
INCLUDE(common/nshel)
INCLUDE(common/infoa)
INCLUDE(common/zorac)
INCLUDE(common/restrl)
c
      REAL ex2, cs2, cp2, cd2, cf2, cg2
      integer kstrt2, ktom2, ktype2, kng2, kloc2, kmin2, kmax2
      integer nshel2, nspace
c
      common/junk/ilifc3(maxorb),ntran3(maxorb),itran3(mxorb3),
     * ctran3(mxorb3),itrn3(2),
     *ex2(mxprim),cs2(mxprim),cp2(mxprim),cd2(mxprim),cf2(mxprim),
     *cg2(mxprim),cspace(maxat),
     +kstrt2(mxshel),ktom2(mxshel),ktype2(mxshel),kng2(mxshel),
     +kloc2(mxshel),kmin2(mxshel),kmax2(mxshel),nshel2,nspace(3),
     *nprin1,itol1,icut1,normf1,normp1
     *,nprini(695),iovmat,iosvec,ioproj,iorthg,c2(3,maxat)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)
INCLUDE(common/symtry)
INCLUDE(common/iofile)
c
c.... copy nshel to junk
c
      call dcopy(6*mxprim,ex,1,ex2,1)
_IF1(c)      call dcopy(7*mxshel+4,kstart,1,kstrt2,1)
_IFN1(c)      call icopy(7*mxshel+4,kstart,1,kstrt2,1)
      normp1 = 0
      normf1 = 0
c
c.... and coordinates from infoa
c
      call dcopy(3*maxat,c,1,c2,1)
c
      call swap_zora('internal')
c
      opbas = .false.
      odbas = .false.
      ofbas = .false.
      ogbas = .false.
      do 100 loop = 1,nshell
      if (ktype(loop) .eq. 2) opbas = .true.
      if (ktype(loop) .eq. 3) odbas = .true.
      if (ktype(loop) .eq. 4) ofbas = .true.
      if (ktype(loop) .eq. 5) ogbas = .true.
  100 continue 
c
c.... copy coordinates to nshel (internal)
c
      call dcopy(maxat*3,c2,1,c,1)
c
      oatint_z= nat_z.gt.0
c...    atomic scf requires atomic transformation 
c...    (not defined with full coulomb) (but who am i ....)
c
      call scross(sout,mdim)
c
c
c----- read in transformation matrices for s,p,d,f and g basis functions.
c
      call rdedx(ptr,nw196(1),ibl196(1),idaf)
      if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
      if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
      if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
c
c .... read seems to be necessary. scross writes on bufb
c
      call standv(1,q)
c
c
      call sder('zora',q)
c
      oatint_z = .false.
      call swap_zora('external')
c
      return
      end
c
c
c
c*********************************************************************
c
      subroutine zora(q,dens,t1,mode)
c
c this version is called in case of oprojz      
c projected zora version (only remaining)
c S. Faas (Utrecht,1998)
c
c     q : blank common scratch
c     dens : current density matrix (only input) (external triangular)
c     t1: external 1-electron integrals to be corrected
c         if oscale t1 is the scalingfactor (external triangular)
c     mode  : 'start' : startup mode (no density matrix)
c           : 'force' : forced recalculaion of zora corrections
c           : 'calc'  : calculate zora corrections if not available 
c           : 'scale' : calculate scaling
c           : 'read'  : requires zora corrections be available on disk
c           : 'read1' : requires zora corrections be available on disk 
c                       they are read in dens in case q is not available
c           : 'dump'  : dump the zora corrctions on ed7 to section 425 
c                       (for atomic corrections)
c
      implicit REAL    (a-h,p-z), integer (i-n), logical (o)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
c
INCLUDE(common/zorac)
INCLUDE(common/restri)
c 
      dimension q(*),dens(*),t1(*)
      character*(*) mode
c
c     parameter (c1  = 137.0359895d0)
      c1 = cspeed
      call secinf(idum,numdu,idum,idum)
c
      if (mode.eq.'dump') then
         izzz = igmem_alloc(nwcor_z)
         i425 = 425
         ll = lensec(nwcor_z)
         call secput(isect(425),i425,ll,iblock_z)
         call rdedx(q(izzz),nwcor_z,ibcor_z,num8)
         call wrt3(q(izzz),nwcor_z,iblock_z,numdu)
         irest_z = 1
         call gmem_free(izzz)
         return
      end if
c
      if (mode.eq.'start'.or.irest_z.eq.2) then
         if (mode.eq.'scale') call caserr('jvl confused')
         if (irest_z.eq.2) then
            i425 = 425
            izzz = igmem_alloc(nwcor_z)
            if (numdu_z.ne.0) then
               call revise
               call secini(ibldu_z,numdu_z)
               call sectst(isect(425),iretrn)
               if (iretrn.eq.0) 
     1         call caserr('foreign dumpfile contains no zora')
               call secget(isect(425),i425,k)
               call rdedx(q(izzz),nwcor_z,k,numdu_z)
               call secini(ibl3d,idaf)
               ll = lensec(nwcor_z)
               call secput(isect(425),i425,ll,iblock_z)
               call wrt3(q(izzz),nwcor_z,iblock_z,numdu)
            end if
            call sectst(isect(425),iretrn)
            if (iretrn.eq.0) call caserr('dumpfile has no zora data')
            call secget(isect(425),i425,iblock_z)
            call rdedx(q(izzz),nwcor_z,iblock_z,numdu)
            call vadd(q(izzz),1,t1,1,t1,1,nwcor_z)
            call gmem_free(izzz)
            opre_zora = .true.
            irest_z = 1
            return
         else if (irest_z.eq.1) then
            izzz = igmem_alloc(nwcor_z)
            call rdedx(q(izzz),nwcor_z,iblock_z,numdu)
            call vadd(q(izzz),1,t1,1,t1,1,nwcor_z)
            call gmem_free(izzz)
            return
         else
            opre_zora = .false.
         end if
      end if
c
      tcoul = -1.0d0
      if ((mode.eq.'calc'.or.mode(1:4).eq.'read').and.opre_zora) then
         if (mode.eq.'read1') then
            if (irest_z.eq.0) then
               call rdedx(dens,nwcor_z,ibcor_z,num8)
            else
               call rdedx(dens,nwcor_z,iblock_z,numdu)
            end if
c...  add correction to 1-electron matrix
            call vadd(dens,1,t1,1,t1,1,nwcor_z)
         else
            isd1 = igmem_alloc(nwcor_z)
            if (irest_z.eq.0) then
               call rdedx(q(isd1),nwcor_z,ibcor_z,num8)
            else
               call rdedx(q(isd1),nwcor_z,iblock_z,numdu)
            end if
c...  add correction to 1-electron matrix
            call vadd(q(isd1),1,t1,1,t1,1,nwcor_z)
            call gmem_free(isd1)
         end if
         return
      end if
c
      if (mode.eq.'read') 
     1 call caserr('zora corrections not on disk for read')
      if (irest_z.ne.0) call caserr('zora restore inhibits calc/force')
c
      if (.not.opre_zora) then
         call pre_zora(q)
         opre_zora = .true.
      end if
c
c...  get "external" density
c
      if (mode.ne.'scale'.and.mode.ne.'start') then
         idext = igmem_alloc(num*(num+1)/2)
         if (is_z.ne.0) then
c...        get alien density from abroad, if requested
            call dens_z(q,q(idext))
         else
c...        copy external density
            call dcopy(num*(num+1)/2,dens,1,q(idext),1)
c
         end if
      end if
c
      numold = num
      call swap_zora('internal')
c
      l1 = num
      l2 = num*(num+1)/2
      l3 = num*num
      if (opzora) then
      if (mode.eq.'scale') write(iwr,*) ' adding ZORA scaling '
      end if
c
c set pointers
c
      idens = igmem_alloc(l3)
      isd1  = igmem_alloc(l3) 
c
      if (mode.eq.'scale') then
c
c...     the coulomb matrix was just calculated so read it 
c
         call rdedx(q(isd1),nwv_z,ibcoul_z,num8)
c
      else if ((icoul_z.eq.1.or.icoul_z.eq.2)
     1         .and.mode.ne.'start') then
c
c...    save the density matrix (external basis)
c
         call wrt3(q(idext),nwdens_z,ibdens_z,num8)
c 
         iout2 = igmem_alloc(l3)
         call square(q(idens),q(idext),numold,numold)
c
c...   read transformation large to small
c
         call rdedx(q(iout2),nwtrat_z,ibtrout_z,num8)
c
c...  transform density matrix to internal (large) basis
c
         call mxmaa(q(iout2),1,num,q(idens),1,numold,q(isd1),1,num,
     *              num,numold,numold)
c
         call mxmaa(q(isd1),1,num,q(iout2),num,1,q(idens),1,num,
     *              num,numold,num)
c
c store density matrix in triangular form
c
         call triangle(q(idens),q(idens),num)
c
         call gmem_free(iout2)
c
c...   generate internal coulomb matrix
c...   density matrix at idens, resulting coulomb matrix at isd1
c
         tcoul = cpulft(1)
         call coul_z(q,idens,isd1)
         tcoul = cpulft(1)-tcoul
c
      else if (icoul_z.eq.3.or.mode.eq.'start') then
c...     **no** coulomb matrix
        call vclr(q(isd1),1,l2)
c
      else if (icoul_z.eq.4) then
c...     ** coulomb matrix in external basis (**debug**))
         call swap_zora('external')
c
         call dcopy(num*(num+1)/2,q(idext),1,q(idens),1)
         call swap_zora('small')
         tcoul = cpulft(1)
         call coul_z(q,idens,isd1)
         tcoul = cpulft(1)-tcoul
         call swap_zora('nonsmall')
c
         call swap_zora('internal')
c...     project coulomb matrix to internal basis
         iout2 = igmem_alloc(l3)
c
         call square(q(idens),q(isd1),numold,numold)
c
c...   read transformation small to large (exactly opposite to dens)
c
         call rdedx(q(iout2),nwtrat_z,ibtrin_z,num8)
c
c...  transform coulomb matrix to internal (large) basis
c
         call mxmaa(q(iout2),numold,1,q(idens),1,numold,q(isd1),1,num,
     *              num,numold,numold)
c
         call mxmaa(q(isd1),1,num,q(iout2),1,numold,q(idens),1,num,
     *              num,numold,num)
c
c store coulomb matrix in triagular form
c
         call triangle(q(idens),q(isd1),num)
c
         call gmem_free(iout2)
      else
         call caserr('invalid icoul_z in zora')
      end if
c
      if (tcoul.gt.0.0d0) write(iwr,9089) tcoul
9089  format('  ZORA coulomb matrix : ',f14.2,' seconds')
c
      call wrt3(q(isd1),nwv_z,ibcoul_z,num8)
c
c set pointers for remainder
c

      iout2 = igmem_alloc(l3)
      isd2  = igmem_alloc(l3)
      it1   = igmem_alloc(l3)
      ivc   = igmem_alloc(l3)
      is2   = igmem_alloc(l3)
      is4   = is2
      isx   = igmem_alloc(l3)
      isy   = igmem_alloc(l3)
      isz   = igmem_alloc(l3)
      is    = igmem_alloc(l3)
c
c compute vc
c     using a coulomb matrix at q(isd1)
c
      cc2=2.0d0*c1*c1     
c
c compute v
c       
      call rdedx(q(ivc),nwv_z,iblv_z,num8)
      call rdedx(q(is),nwv_z,ibls_z,num8)
c
      if (.not. (mode.eq.'scale')) then
          do 198 i=1,l2
             q(ivc+i-1)=q(is+i-1)-(q(ivc+i-1)+q(isd1+i-1))/cc2
198       continue
          if (op2zora) then
            write(iwr,*)'vc matrix in zora, no scaling' 
            call prtri(q(ivc),num)
          end if 
      else
          do 298 i=1,l2
             q(ivc+i-1)=(-q(ivc+i-1)-q(isd1 + i-1) + q(is+i-1)*cc2)/c1
298       continue
          if (op2zora) then
              write(iwr,*)'vc matrix in zora, scaling'
              call prtri(q(ivc),num)
          end if
      endif
c
c transform to orthogonal basis:
c transform matrix
c 
      call rdedx(q(is2),nwmat_z,ibsmin_z,num8)
      call mult2(q(is2),q(isd2),q(ivc),num,num,num)
c
c compute zora corrections
c
      n2zora=num*num
c     if (.not.(mode.eq.'scale')) then
c
         call rdedx(q(isx),nwmat_z,ibsx_z,num8)
         call rdedx(q(isy),nwmat_z,ibsy_z,num8)
         call rdedx(q(isz),nwmat_z,ibsz_z,num8)
c
      if (.not.(mode.eq.'scale')) then
         call mxmaa(q(isx),num,1,q(isx),num,1,q(it1),num,1,
     +              num,num,num)
c
         call mxmb(q(isy),num,1,q(isy),num,1,q(it1),num,1,
     +             num,num,num)
c
         call mxmb(q(isz),num,1,q(isz),num,1,q(it1),num,1,
     +             num,num,num)
c
         call dscal(l3,-1.0d0,q(it1),1)
c
c...  it1 contains the  -sum(ig) <f|p_i|g><g|p_i|f>
c
      else
        call vclr(q(it1),1,l3)
      end if
c make vc square for inversion
c
      call square(q(ivc),q(isd2),num,num)
c           
      d=-1.0d0
      call minvrt(q(ivc),num,d,q(isd2),q(isd1))
c
      if (mode.eq.'scale') then
         call mxmaa(q(ivc),num,1,q(ivc),num,1,q(isd1),num,1,
     +              num,num,num)
         call dcopy(num*num,q(isd1),1,q(ivc),1)
      endif
c
c compute term 2
c
      call mxmaa(q(isx),num,1,q(ivc),num,1,q(isd1),num,1,
     +           num,num,num)
      call mxmb(q(isd1),num,1,q(isx),num,1,q(it1),num,1,
     +           num,num,num)
c
      call mxmaa(q(isy),num,1,q(ivc),num,1,q(isd1),num,1,
     +           num,num,num)
      call mxmb(q(isd1),num,1,q(isy),num,1,q(it1),num,1,
     +             num,num,num)
c
      call mxmaa(q(isz),num,1,q(ivc),num,1,q(isd1),num,1,
     +           num,num,num)
      call mxmb(q(isd1),num,1,q(isz),num,1,q(it1),num,1,
     +          num,num,num)
c
      if (.not.(mode.eq.'scale')) then
         call dscal(l3,-0.5d0,q(it1),1)
c
c    q(it1) :  0.5*( -sum(ig) <f|p_i|g><g|p_i|f> +
c                    +sum(igk) <f|p_i|g><g|1-v/2c**2|k>**-1<k|p_i|f> )
      endif
c
c transform back to original basis
c
      call rdedx(q(is4),num*num,ibsplus_z,num8)
      call mxmaa(q(it1),1,num,q(is4),1,num,q(isd1),1,num,
     1           num,num,num)
      call mxmaa(q(is4),num,1,q(isd1),1,num,q(it1),1,num,
     1           num,num,num)
c
c transform back to original external basis
c
      call rdedx(q(iout2),nwtrat_z,ibtrout_z,num8)
      if (op2zora) then
         write(iwr,*) 'kinetic corrections before trafo:'
         call prsq(q(it1),num,num,num)
         write(iwr,*) 'sout2 in zora;'
         call prsq(q(iout2),num,numold,numold)
      endif
c
      call mxmaa(q(it1),1,num,q(iout2),1,num,q(isd1),1,num,
     1           num,num,numold)
      call mxmaa(q(iout2),num,1,q(isd1),1,num,q(it1),1,numold,
     1           numold,num,numold)
c     
c...  now we are back in external basis
c...  => num = numold
c
      call swap_zora('external')
c
c make triangular again
c
      call triangle(q(it1),q(isd1),num)
c
c     **NOW** we have in q(isd1) the correction to the kinetic energy
c 
      if (op2zora) then
         if (.not.(mode.eq.'scale')) then
            write(iwr,*) 'zora corrections in zora'
         else
            write(iwr,*) 'back transformed scaling factors:'
         end if
         call prtril(q(isd1),num)
      endif
c
c...   write either the kinetic energy correction or the ao-basis scaling
c 
      if (.not.(mode.eq.'scale')) then
         call nonrel_atom(q(isd1))
         call wrt3(q(isd1),nwcor_z,ibcor_z,num8)
c...  add correction to 1-electron matrix
         call vadd(q(isd1),1,t1,1,t1,1,nwcor_z)
      else
         call wrt3(q(isd1),nwcor_z,ibscale_z,num8)
c...  move scaling to t1
         call dcopy(nwcor_z,q(isd1),1,t1,1)
      end if
c 
      call gmem_free(is)
      call gmem_free(isz)
      call gmem_free(isy)
      call gmem_free(isx)
      call gmem_free(is2)
      call gmem_free(ivc)
      call gmem_free(it1)
      call gmem_free(isd2)
      call gmem_free(iout2)
      call gmem_free(isd1)
      call gmem_free(idens)
      if (mode.ne.'scale'.and.mode.ne.'start') call gmem_free(idext)
c
      return
9088  format (/20x,23('*')/20x,'symmetrized matrix'/20x,23(
     +     '-'))
5030  format (2i4,3f16.9)
      end
c
c **********************************************************************
c
      subroutine ver_zora(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/zora.m,v $
     +     "/
      data revision /"$Revision: 6289 $"/
      data date /"$Date: 2013-10-15 13:24:06 +0200 (Tue, 15 Oct 2013) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
      subroutine basis_zora(q,csinp,cpinp,cdinp,cfinp,cginp,
     * ztita,czana,zatagg,
     *csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     *ict,newsh,iso,isoc,natmax,nshmax,ngsmax,onel)
c
c...  fill in extended basis for zora and calculate iso etc arrays
c...  ** note various commons get messed up; swap zora should save them
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/zorac)
INCLUDE(common/infoa)
c
INCLUDE(common/nshel)
c     dimension newsh(*),iso(*)
      dimension ict(natmax,*),newsh(nshmax,*),iso(*),isoc(*)
      dimension ztita(*),czana(*),zatagg(*),q(*)
      dimension csinpi(ngsmax),cpinpi(ngsmax),cdinpi(ngsmax)
      dimension cfinpi(ngsmax),cginpi(ngsmax)
      dimension csinp(ngsmax),cpinp(ngsmax),cdinp(ngsmax)
      dimension cfinp(ngsmax),cginp(ngsmax)
c
      num_ext = num
      iexpm = igmem_alloc(natmax*5)
      call basexp_zora(q,q,q(iexpm),
     +                 csinpi,cpinpi,cdinpi,cfinpi,cginpi,ngsmax)
      call gmem_free(iexpm)
      call swap_zora('internal')
      num = 0
      call atoms2(csinp,cpinp,cdinp,cfinp,cginp,
     * ztita,czana,zatagg,
     *csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     *ict,newsh,iso,isoc,natmax,nshmax,ngsmax,onel)
      nbas_zora = num
      call swap_zora('external')
c
      if (nbas_zora.lt.num) call caserr
     +   ('Internal basis smaller than external')
c
      return
      end

      subroutine basexp_zora(q,iq,expm,
     +                       csinpi,cpinpi,cdinpi,cfinpi,cginpi,
     +                       ngsmax)
c
c this routine expands original basis set to obtain a suitable internal
c basis set (scalar zora calculation). 
c
c Variable int_zora defines type of addition:
c int_zora = 0: No extra functions added (internal basis is external basis)
c          = 1: Exponents from L-1 block larger than largest exponent
c               in L-block are added automatically to L-block
c          = 2: Extra exponents are read from input
c
c Variables onlys, onlyp, onlyd, onlyf, onlyg define which types of
c functions are added
c onlys = true: No functions are added
c onlyp = true: p-functions
c onlyd = true: p- and d-functions
c onlyf = true: p-, d- and f-functions 
c onlyg = true: p-, d-, f- and g-functions 
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/czmat)
INCLUDE(common/nshel)
INCLUDE(common/zorac)
      common/junk/trmat(432),
     +   ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +   intyp(mxshel),ns(maxat),ks(maxat),
     +   intypi(mxshel),nsi(maxat),ksi(maxat),
     +   expz(mxprim),
     +   csz(mxprim),cpz(mxprim),cdz(mxprim),cfz(mxprim),cgz(mxprim),
     +   kstartz(mxshel),katomz(mxshel),ktypez(mxshel),kngz(mxshel),
     +   klocz(mxshel),kminz(mxshel),kmaxz(mxshel),nshellz
c
      parameter (iz_atom=100,iz_exp=1000)
      common/zora_char/ zoratom(iz_atom),ztyp(iz_atom)
      common/zora_ibas/ exp_zora(iz_exp),
     +                  ib_z(iz_atom),ie_z(iz_atom),izm,ib
c
      logical opresent
      dimension map_zora(iz_atom)
c
      dimension expm(6,maxat)
      dimension q(*),iq(*)
      dimension csinpi(ngsmax),cpinpi(ngsmax),cdinpi(ngsmax),
     +          cfinpi(ngsmax),cginpi(ngsmax)
c
      parameter (pi32=5.56832799683170d0)
c
c search for maximum exponent on each block. and each atom
c
      do i = 1,iz_atom
         map_zora(i) = 0
      enddo
c
      if (op2zora) write(iwr,1) 'Number of Atoms:',nat
 1    format (a16,i4)
c
      do iat = 1,nat
         do i=1,6
            expm(i,iat) = 0.0d0
         end do
      end do
c
      if (op2zora) then
c data from original basis set
         do i=1,nshell
              write(iwr,2) kstart(i),katom(i),ktype(i),kng(i),
     +                     kloc(i),kmin(i),kmax(i)
         enddo
2        format (7(2x,i2))
         write(iwr,*) 'the other integers:'
         write(iwr,*) nshell,non,numorb,ndumm ,num
      endif 
c
      do iat = 1,nat
         istart = 1
         do i=1,nshell
          if (katom(i).eq.iat) then
            iend = istart + kng(i) -1
            do j = istart,iend
             if (ex(j).gt.expm(ktype(i),iat)) expm(ktype(i),iat) = ex(j)
            enddo
          endif
          istart = istart + kng(i)
         enddo 
      enddo
c           
      if (op2zora) then
         do iat = 1,nat
            write (iwr,5)'Maximum exponents for atom:',iat 
            write (iwr,4)'Maximum exponent on s-block is:',expm(1,iat) 
            write (iwr,4)'Maximum exponent on p-block is:',expm(2,iat)
            write (iwr,4)'Maximum exponent on d-block is:',expm(3,iat)
            write (iwr,4)'Maximum exponent on f-block is:',expm(4,iat)
            write (iwr,4)'Maximum exponent on g-block is:',expm(5,iat)
            write (iwr,*)' '
         enddo
      endif
4     format (a31,2x,f14.4)
5     format (a27,2x,i4)
c
c set maxima according to input of limit
c
      if (onlys .or. (int_zora.eq.0)) then
         istart = 2
         iend = 1
         write (iwr,*) 'No functions added'
      elseif (onlyp) then
         istart = 3
         iend = 2
         write (iwr,*) 'No d-, f- and/or g-functions added'
      elseif (onlyd) then
         istart = 4
         iend = 3
         write (iwr,*) 'No f- and/or g-functions added'
      elseif (onlyf) then
         istart = 5
         iend = 4
         write (iwr,*) 'No g-functions added'
      elseif (onlyg) then
c ... h functions not available
         istart = 6
         iend = 5
      endif
c
      do iat = 1,nat
         do itype = istart,6
            expm(itype,iat) = 1.0d99
         enddo
      enddo
c
      if (onoext_z) then
c do not copy external onto internal
         do iat = 1,nat
            do itype = 1,iend
               expm(itype,iat) = 0.0d0
            enddo
         enddo
      endif
c
c augment external basis with high exponent function  
c
      nonz = non
      numorbz = numorb
      ndummz = ndumm
      num2 = 0
c
c count number of external functions
c
      do i=1,nshellz
         num2=num2+kngz(i)
      enddo
c
c initialize memory for temporary saving of additional functions
c
      nipw = lenwrd()
      if (int_zora.ne.0) then
         if (int_zora.eq.1) then
            length = num2
         else if (int_zora.eq.2) then
            length = ib
         endif
         if (nipw.eq.1) then
            i10 = igmem_alloc(length*3)
            i20 = i10 + length
            i30 = i20 + length
         else if (nipw.eq.2) then
            i10 = igmem_alloc(length*2)
            i20 = 2*(i10+length)-1
            i30 = i20 + length
         endif
      endif
c
c generating additional functions 
c
      if (int_zora.eq.1) then
c automatic internal basisset generation
         istart = 1
         nbfz = 0
c       
         do i=1,nshellz
            do iat=1,nat
               if (katomz(i).eq.iat) then
                  iend = istart+kngz(i)-1
                  do j = istart,iend
                     if (expz(j) .gt. expm(ktypez(i)+1,iat)) then 
                        opresent = .false.
                        do ib = 1, nbfz
                           opresent = opresent .or. 
     +                       ((q(i10+ib-1)  .eq. expz(j)).and.
     +                        (iq(i20+ib-1) .eq. ktypez(i)+1).and.
     +                        (iq(i30+ib-1) .eq. iat))
                        enddo
                        if (.not.opresent) then
                           nbfz = nbfz + 1
                           q(i10+nbfz-1) = expz(j)
                           iq(i20+nbfz-1) = ktypez(i)+1
                           iq(i30+nbfz-1) = iat
                        endif
                     endif
                  enddo
               endif
            enddo
            istart = istart + kngz(i)
         enddo
c
      else if (int_zora.eq.2) then
c manual basis set generation
         nbfz = 0
         do i=1,izm
            do iat=1,nat
               if ((zoratom(i).eq.zaname(iat)) .and.
     +             (map_zora(i).eq.0)) then
                  map_zora(i) = 1
                  do j = ib_z(i),ie_z(i)
                     if (ztyp(i).eq.'s') then
                        nbfz = nbfz + 1
                        q(i10+nbfz-1) = exp_zora(j)
                        iq(i20+nbfz-1) = 1
                        iq(i30+nbfz-1) = iat
                     else if (ztyp(i).eq.'p') then
                        nbfz = nbfz + 1
                        q(i10+nbfz-1) = exp_zora(j)
                        iq(i20+nbfz-1) = 2
                        iq(i30+nbfz-1) = iat
                     else if ((ztyp(i).eq.'d').and.(.not.onlyp)) then
                        nbfz = nbfz + 1
                        q(i10+nbfz-1) = exp_zora(j)
                        iq(i20+nbfz-1) = 3
                        iq(i30+nbfz-1) = iat
                     else if ((ztyp(i).eq.'f').and.(.not.onlyp)
     +                         .and.(.not.onlyd)) then
                        nbfz = nbfz + 1
                        q(i10+nbfz-1) = exp_zora(j)
                        iq(i20+nbfz-1) = 4
                        iq(i30+nbfz-1) = iat
                     else if ((ztyp(i).eq.'g').and.(.not.onlyp)
     +                         .and.(.not.onlyd).and.(.not.onlyf)) then
                        nbfz = nbfz + 1
                        q(i10+nbfz-1) = exp_zora(j)
                        iq(i20+nbfz-1) = 5
                        iq(i30+nbfz-1) = iat
                     endif
                  enddo 
               endif
            enddo
         enddo
c
      endif
c
      if (onoext_z) then
         nshellz = 0
         num2 = 0
      endif
c
      if (int_zora.ne.0) then
         call add_zora(q(i10),iq(i20),iq(i30),nbfz,num2,
     +                 csinpi,cpinpi,cdinpi,cfinpi,cginpi,ngsmax)
         call gmem_free(i10)
      endif
c
      if (op2zora) then
         istart = 1
         do j=1,nshellz
            iend = istart + kngz(j) - 1 
            do i=istart,iend
               write(iwr,7) expz(i),csz(i),cpz(i),cdz(i),cfz(i),cgz(i)
            enddo 
            istart = istart + kngz(j)
            write(iwr,*) istart
         enddo
7        format (6f12.4)
         do i=1,nshellz
            write(iwr,8) ktypez(i),klocz(i),kstartz(i),kminz(i),
     +                   kmaxz(i),kngz(i)
         enddo
8        format (6i4)
      endif
c
      return
      end
c
      subroutine swap_zora(aim)
c
c swaps basis set information internal and external basis 
c in common blocks nshel (external) and nshel_z (external)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/zorac)
INCLUDE(common/nshel)
      integer onormf1,onormp1
      common/nshel_z/exp_z(mxprim),
     +               cs_z(mxprim),cp_z(mxprim),cd_z(mxprim),
     +               cf_z(mxprim),cg_z(mxprim),
     +               kstart_z(mxshel),katom_z(mxshel),ktype_z(mxshel),
     +               kng_z(mxshel),kloc_z(mxshel),kmin_z(mxshel),
     +               kmax_z(mxshel),nshell_z,non_z,numorb_z,ndumm_z
     +              ,onormf1,onormp1,c2(3,maxat)
INCLUDE(common/infob)
c
c from March 97, amass is used only for user-input
c mass values
c
      integer nonsym_z, map80_z
      common/infob_z/nonsym_z,map80_z(maxat)
c
INCLUDE(common/symtry)
      integer iliso_z,ilisoc_z
      common/symtry_z/iliso_z(48),ilisoc_z(48)
c
INCLUDE(common/infoa)
      integer int_z
c ... int_z contains num and nx (changed in atoms2)
      REAL czan_z,c_z
      common/infoa_z/int_z(2),czan_z(maxat),c_z(3,maxat)
c
INCLUDE(common/cslosc)
      common/cslosc_z/odscf_z,ifock_z,idmat_z,irdmat_z,iprefa_z,iexch_z
c
INCLUDE(common/runlab)
      common/runlab_z/zscftp_z
c
INCLUDE(common/restar)
      common/restar_z/nprint_z,nopk_z,intg76_z,m2file_z(61)
c
INCLUDE(common/transf)
      common/transf_z/psmal_z,qsmal_z,rsmal_z,
     +                pnew_z,qnew_z,rnew_z,pp_z,qp_z,rp_z
c
INCLUDE(common/restrl)
c      (just these 4 jvl/sf 1999)
      logical opbas_z,odbas_z,ofbas_z,ogbas_z
      common /restrl_z/ opbas_z,odbas_z,ofbas_z,ogbas_z
c
c
      character*(*) aim
      character*8 state
      save state
      data state/'external'/
c
      if (state.eq.aim) then
         write(iwr,'(a,a)') ' current swap state is ',aim
         call caserr(' swap_zora error')
      end if
c
c...  for small coulomb swapping (in the end debug)
c
      if (aim.eq.'small'.or.aim.eq.'nonsmall') then
         if (state.ne.'external') call caserr('aim small for no extern')
         if (osmall_zora.and.aim.eq.'small') call caserr('swap  small')
         if (.not.osmall_zora.and.aim.eq.'nonsmall') call caserr('swap')
         if (aim.eq.'small') osmall_zora = .true.
         if (aim.eq.'nonsmall') osmall_zora = .false.
c...  runlab
         zzzz = zscftp_z
         zscftp_z = zscftp
         zscftp = zzzz
c...  restar
         call iswap(1,nprint,1,nprint_z,1)
         call iswap(1,nopk,1,nopk_z,1)
         call iswap(1,intg76,1,intg76_z,1)
         call iswap(61,m2file,1,m2file_z,1)
c...  cslosc
         call iswap(5,odscf,1,odscf_z,1)
         call iswap(1,iexch,1,iexch_z,1)
         return
      end if
c
      state = aim
c
      if (state.eq.'internal') then
        oint_zora = .true.
      else
        oint_zora = .false.
      end if
c    
c...  nshell
c
      call dswap(mxprim*6,exp_z,1,ex,1)
      call iswap(mxshel*7+4,kstart_z,1,kstart,1)
c
c...  infob
c
c     call iswap(maxat+1,nonsym,1,nonsym_z,1)
      call iswap(maxat,map80,1,map80_z,1)
c 
c...  symtry
c
      call iswap(96,iliso,1,iliso_z,1)
c
c...  infoa
c
      call iswap(2,num,1,int_z,1)
      call dswap(maxat*4,czan,1,czan_z,1)
c
c...  cslosc
c
      call iswap(5,odscf,1,odscf_z,1)
      call iswap(1,iexch,1,iexch_z,1)
c
c...  runlab
c
      zflop = zscftp
      zscftp = zscftp_z
      zscftp_z = zflop
c
c...  restar
c
      call iswap(1,nprint,1,nprint_z,1)
      call iswap(1,nopk,1,nopk_z,1)
      call iswap(1,intg76,1,intg76_z,1)
      call iswap(61,m2file,1,m2file_z,1)
c
c...  transf
c
      call dswap(9,psmal,1,psmal_z,1)
c
c...  restrl
c
      call iswap(4,opbas,1,opbas_z,1)
c
      return
      end
c
      subroutine coul_z(q,idens,icoul)
c
c     direct coulomb matrix builder for ZORA
c     density matrix at q(idens)
c     resulting coulomb matrix at q(icoul)
c     restarting disabled
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension q(*)
c
INCLUDE(common/sizes)
INCLUDE(common/cslosc)
INCLUDE(common/infoa)
INCLUDE(common/restri)
INCLUDE(common/restrl)
      common/scra/iso(mxshel,48)
INCLUDE(common/nshel)
INCLUDE(common/iofile)
INCLUDE(common/runlab)
INCLUDE(common/zorac)
INCLUDE(common/sortp)
INCLUDE(common/symtry)
INCLUDE(common/restar)
INCLUDE(common/mapper)
      common/bufb/ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720)  
c
c...  read iso ...
c
      nav = lenwrd()
      if (oint_zora) then
         call readi(iso,nwiso_z*nav,ibiso_z,num8)
      else
         call readi(iso,nw196(5)*nav,ibl196(5),idaf)
      end if
c
c...  and dtr etc ..... to be sure 
c
      if (nt.gt.1) then
        if (odbas) call rdedx(dtr,nw196(2),ibl196(2),idaf)
        if (ofbas) call rdedx(ftr,nw196(3),ibl196(3),idaf)
        if (ogbas) call rdedx(gtr,nw196(4),ibl196(4),idaf)
      endif
c ...
c     we must reset "zscftp" in /runlab/ to ensure
c     the closed shell fock builder is called from intega etc
c
      zscftp = 'coul'
      odscf = .true.
      nopk = 1
      oschw_save = oschw
      oschw = .false.
      intg76 = 0
      oatint_z = icoul_z.eq.2
      if (.not.opzora) nprint = -5
c
c...  oschw+intg76 => only jkinta
c..
c..   core for hstar separate (i10/iexch/idmat)
c
      out = .false.
c
      l1 = num
      l2 = num*(num+1)/2
c
      ifock = icoul
      idmat = idens
c
      irdmat = igmem_alloc(nshell*(nshell+1)/2)
      iprefa = igmem_alloc(nshell*(nshell+1)/2)
      i20    = igmem_alloc(l2)
      iexch  = i20
c
c...  nothing is written because no restart possible
c...  build prefactor and reduced density matrix 
c
      call rdmake(q(iprefa))
      call mkrdmt('rhf',q(irdmat),q(idmat),l2,nprint)
c
c     ----- construct a skeleton coulomb matrix -----
c     h       at q(i10)
c     density at q(idmat)
c     exch    at q(iexch) (scratch)
c
      call dhstar(q,q(ifock),q(idmat),q(iprefa),q(irdmat),0)
c
      if (out) then
         write (iwr,6060)
         call prtri(q(ifock),l1)
      end if
c
c     ----- symmetrize skeleton fock matrix -----
c
c     -h- at q(ifock)
c     scratch area at q(i20)
c
      call symh(q(ifock),q(i20),iky,0,0)
c
      if (out) then
         write (iwr,6070)
         call prtri(q(ifock),l1)
      end if
c...    reset
      oschw = oschw_save
      call readi(iso,nw196(5)*nav,ibl196(5),idaf)
      oatint_z = .false.
c
c...   reset core
c
      call gmem_free(i20)
      call gmem_free(iprefa)
      call gmem_free(irdmat)
c
      return
 6060 format (/20x,20('-')/20x,'skeleton coulomb matrix'/20x,20('-'))
 6070 format (/20x,23('-')/20x,'symmetrized coulomb matrix'/20x,
     +        23('-'))
      end
      subroutine antisquare(r,a,n)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c..   convert a triangular matrix (with diagonal)
c..   into an antisymmetric square (with 0.0 diagonal)
c
      dimension r(n,n),a(*)
c
      k = 0
      do i=1,n
       do j = 1,i-1
          k = k + 1
          r(i,j) = a(k)
          r(j,i) = -a(k)
       end do
       k = k + 1
       r(i,i) = 0.0d0
      end do
c
      return
      end
c
      subroutine zora_int
c
c...  read the specification for the internal basis
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
c
      parameter (iz_atom=100,iz_exp=1000)
      common/zora_char/ zoratom(iz_atom),ztyp(iz_atom)
      common/zora_ibas/ exp_zora(iz_exp),
     +                  ib_z(iz_atom),ie_z(iz_atom),izm,ib
c
      izm = 0
      ib = 0
c
10    call inpa4(ystring)
      if (ystring.ne.'end') then
         izm = izm + 1
         if (izm.gt.iz_atom) call caserr('zoratom overflow')
         zoratom(izm) = ystring
         call inpa(ztyp(izm))
         ib_z(izm) = ib + 1
         call input
20       call inpf(real)
         if (real.ne.0.0d0) then
            ib = ib + 1
            if (ib.gt.iz_exp) call caserr('exp_zora overflow')
            exp_zora(ib) = real
            go to 20
         else
            ie_z(izm) = ib
            call input
            go to 10
         end if
c
      else
         call input
      end if
c
      return
      end
c
      subroutine add_zora(expdummy,itypdummy,iatdummy,n,num2,
     +                    csinpi,cpinpi,cdinpi,cfinpi,cginpi,ngsmax)
c
c...  generate internal basis by adding functions to external
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/zorac)
      common/junk/trmat(432),
     +   ptr(3,144),dtr(6,288),ftr(10,480),gtr(15,720),
     +   intyp(mxshel),ns(maxat),ks(maxat),
     +   intypi(mxshel),nsi(maxat),ksi(maxat),
     +   expz(mxprim),
     +   csz(mxprim),cpz(mxprim),cdz(mxprim),cfz(mxprim),cgz(mxprim),
     +   kstartz(mxshel),katomz(mxshel),ktypez(mxshel),kngz(mxshel),
     +   klocz(mxshel),kminz(mxshel),kmaxz(mxshel),nshellz
c
      dimension expdummy(*),itypdummy(*),iatdummy(*)
      dimension csinpi(ngsmax),cpinpi(ngsmax),cdinpi(ngsmax),
     +          cfinpi(ngsmax),cginpi(ngsmax)
c
      parameter (pi32=5.56832799683170d0)
      parameter (d43=4.0d0/3.0d0, dinv18=1.0d0/1.875d0,
     1           dinv65=1.0d0/6.5625d0)
c
      do i=1,n
         num2 = num2 + 1
         if (num2.gt.mxprim)
     +      call caserr('Too many primitives added in internal basis')
         nshellz = nshellz + 1
         if (nshellz.gt.mxshel)
     +      call caserr('Too many shells in internal basis')
c
         expz(num2)   = expdummy(i)
         ee           = expz(num2) + expz(num2)
c...     dummys now used inverted ...
c...     avoid squaring and root-taking
         dummys       = (ee*dsqrt(ee))/pi32
         dummyt       = (dsqrt(ee))/pi32
         if (itypdummy(i).eq.1) then
            csz(num2)    = dsqrt(dummys)
            csinpi(num2) = 1.0d0
         else if (itypdummy(i).eq.2) then
            dummy        = 2.0d0*dummyt
            cpz(num2)    = dsqrt(dummy)*ee
            cpinpi(num2) = 1.0d0
         else if (itypdummy(i).eq.3) then
            dummy        = d43*dummys
            cdz(num2)    = dsqrt(dummy)*ee
            cdinpi(num2) = 1.0d0
         else if (itypdummy(i).eq.4) then
            dummy        = dinv18*dummyt
            cfz(num2)    = dsqrt(dummy)*ee*ee
            cfinpi(num2) = 1.0d0
         else if (itypdummy(i).eq.5) then
            dummy        = dinv65*dummys
            cgz(num2)    = dsqrt(dummy)*ee*ee
            cginpi(num2) = 1.0d0
         endif   
c
         if (onoext_z.and.(nshellz.eq.1)) then
            kstartz(nshellz) = 1
         else
            kstartz(nshellz) = kstartz(nshellz-1) + kngz(nshellz-1)
         endif
         katomz(nshellz)  = iatdummy(i)
         ktypez(nshellz)  = itypdummy(i)
         kngz(nshellz)    = 1
         if (ktypez(nshellz).eq.1) then
            intypi(nshellz) = 16
         else if (ktypez(nshellz).eq.2) then
            intypi(nshellz) = 17
         else if (ktypez(nshellz).eq.3) then
            intypi(nshellz) = 18
         else if (ktypez(nshellz).eq.4) then
            intypi(nshellz) = 19
         else if (ktypez(nshellz).eq.5) then
            intypi(nshellz) = 20
         endif
c
         if (onoext_z.and.(nshellz.eq.1)) then
            itypez = 1
         else
            itypez = ktypez(nshellz-1)
         endif
         if (itypez .eq. 1) then
              if (nshellz.eq.1) then
                 klocz(nshellz) = 1
              else
                 klocz(nshellz) = klocz(nshellz-1) + 1
              endif
         else if (itypez .eq. 2) then
              ikmin = kminz(nshellz-1)
              if (ikmin .eq. 1) then
                 klocz(nshellz) = klocz(nshellz-1) + 4
              else if (ikmin .eq. 2) then
                 klocz(nshellz) = klocz(nshellz-1) + 3
              else
                 call caserr ('unknown function type 2')
              endif
         else if (itypez .eq. 3) then
              klocz(nshellz) = klocz(nshellz-1) + 6
         else if (itypez .eq. 4) then
              klocz(nshellz) = klocz(nshellz-1) +10
         else if (itypez .eq. 5) then
              klocz(nshellz) = klocz(nshellz-1) +15
         else
              call caserr ('unknown function type')
         endif
         if (ktypez(nshellz).eq.1) then
            kminz(nshellz) = 1
            kmaxz(nshellz) = 1
         else if (ktypez(nshellz).eq.2) then
            kminz(nshellz) = 2
            kmaxz(nshellz) = 4
         else if (ktypez(nshellz).eq.3) then
            kminz(nshellz) = 5
            kmaxz(nshellz) = 10
         else if (ktypez(nshellz).eq.4) then
            kminz(nshellz) = 11
            kmaxz(nshellz) = 20
         else if (ktypez(nshellz).eq.5) then
            kminz(nshellz) = 21
            kmaxz(nshellz) = 35
         endif
      enddo
c
      return
      end
c
c
c**********************************************************************
c
c
      subroutine scf_so(q)
c
c driver for spin orbit zora calculation. Simon Faas may 1998
c Henk Eshuis 2004
c level shifting enabled 
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      logical ocrit,hm
      COMPLEX alpha,beta
      double precision ode,dzero,done
c 20/6/07      parameter (m100 =400000)
      parameter (m100 =900000)
c
INCLUDE(common/sizes)
INCLUDE(common/mapper)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/scra7)
INCLUDE(common/dump3)
INCLUDE(common/scfopt)
INCLUDE(common/zorac)
INCLUDE(common/occnr_z)
INCLUDE(common/harmon)
INCLUDE(common/atmol3)
c
      COMPLEX mat1,mat2,sum_eng
      character*6 type_z
      common/somatrix1/mat1(m100),mat2(m100)
      common/somatrix2/rmat1(m100),rmat2(m100)
c
      data oconvd,ospinbaan /.false.,.true. / 
      data dzero,done /0.0d0,1.0d0/
c
      dimension q(*)
c
      write(iwr,200)
c      type_z = 'closed'
      type_z = 'open'
      if (no_z.ne.0) type_z= 'open' 
c
      l1 = num
      l0 = newbas0
      l2 = l1*(l1+1)/2
      l3 = l1*l1
      if (l1*l1 .gt. m100) call caserr('spinorbit memory problem')
      alpha = (1.0d0,0.0d0)
      beta = (0.0d0,0.0d0)
c

      ifock = igmem_alloc(8*l3)
      ivect = igmem_alloc(8*l3)
      idnew = igmem_alloc(8*l3)
      idreal= igmem_alloc(l3)
      isd1  = igmem_alloc(l3)
      ieaai = igmem_alloc(l3)
      ieba  = igmem_alloc(l3)
      iebai = igmem_alloc(l3)
      idens = igmem_alloc(l2)
      ijandk= igmem_alloc(l2)
      jblkf = igmem_alloc(l2)
      iw1   = igmem_alloc(2*l1)
c
      if (type_z.eq.'open') then
         iebbi = igmem_alloc(l3)
         idrealb = igmem_alloc(l3)
         ijandkb = igmem_alloc(l2)
         ida = igmem_alloc(l2)
         idb = igmem_alloc(l2)
         ivt = igmem_alloc(l3)
         i80 = igmem_alloc(l1)
         i81 = igmem_alloc(l1)
         i10 = igmem_alloc(l3)
         i50 = igmem_alloc(l3)
         i20 = igmem_alloc(l3)
         iha = igmem_alloc(l2)
         ihb = igmem_alloc(l2)
         i30 = igmem_alloc(l3)
      end if
c
c ... set up occupation numbers for kramers degenerate eigv
c not available yet
c
       call occ_nr
c
c....  make orthonormalising transformation 
c....  iscr is scratch now
c
      is2 = igmem_alloc(l3)
      isq = igmem_alloc(l3)
      itr = igmem_alloc(l2)
      iscr= igmem_alloc(l2)
c
      call rdedx(q(iscr),num*(num+1)/2,ibl7s,num8)
c
      if(op2zora) then
       write(iwr,9008)
       write(iwr,*) 'before tranp: q(iscr)'
       call prtri(q(iscr),l1)
      endif

c ....tranp creates symmetry adapted orbitals

      call tranp(q(iscr),q(itr))

      if(op2zora) then
       write(iwr,9008)
       write(iwr,*) 'q(itr) after tranp'
       call prtri(q(itr),l1)
      endif
c
c generate orthonormalising transformation
c
c .... disable section isecqm to force calculation in qmat
c
      call secdrp(isecqm)
      call qmat(q,q(itr),q(is2),q(isq),q(isd1),iky,l0,l1,l3,l1,
     +          op2zora)
      call dcopy(l3,q(is2),1,q(i20),1)
c      
      if (op2zora) then
         write(iwr,9008) 
         write(iwr,*) 'q(i20) after qmat'
         call prsq(q(i20),l1,l1,l1)
      end if
c
      call secdrp(isecqm)
      call tdown(q(is2),ilifq,q(is2),ilifq,l0)
      if (op2zora) then
         write(iwr,*) 'q(is2) after tdown'
         call prsq(q(is2),l1,l1,l1)
      end if
c
c build total orthonormalising transformation in ivect
c ivect = complex (2*l1 x 2*l1)
c
      call totals2(q(is2),q(ivect),l0,l1)
      
      if(op2zora) then
        write(iwr,*)'orthonormalising transformation q(ivect)'
        call prevz(q,q(ivect),2*l1)
      end if
c
      call gmem_free_set(is2,iscr)
c
c     check for orthonormality
      if (op2zora)then
       write(iwr,*) 'Checking orthonormality of vectors'
       call ortho_check_z(q,q(ivect),l1,l0)      
      endif 
c
c start scf
c
      iter = 1
      iter2=0
      nshift = 0
100   continue 
c
      call vclr(q(ifock),1,8*l3) 
c
c read in 1e part of fock matrix
c      
      call rdedx(q(jblkf),l2,ibl7f,num8)
c
      if (op2zora) then
	  write(iwr,*) '1st part fock-matrix in cs'
	  call prtri(q(jblkf),num)
      endif
c
c build density out of previous iter.
      if (iter .gt. 1) then
	 call zgemm('n','c',2*l1,2*l1,ne,(1.0d0,0.0d0),q(ivect),2*l1
     +           ,q(ivect),2*l1,(0.0d0,0.0d0),q(idnew),2*l1)

        if(op2zora) then
         itt  = igmem_alloc(l3)
         call densdecomplex(q(idnew),q(itt),l1)
         write(iwr,*) 'Electron density, begin: ',iter
         call prsq(q(itt),l1,l1,l1)
         call gmem_free(itt)
         
         write(iwr,*) 'eigenvectors q(ivect), begin',iter
         call prevz(q,q(ivect),2*l1)
         write(iwr,*) 'q(idnew), begin',iter
         call prevz(q,q(idnew),2*l1)
       endif

      else

c ....create initial density
c ....the  old way for closed shell 
c ....build from guess-vectors for open shell (uhf-like)
c ....closed

      if (op2zora) then
          write(iwr,9008)
          write(iwr,*)'Creating initial density, type_z is: ',type_z
      end if

      if (type_z .eq. 'closed') then
        call rdedx(q(idens),l2,ibl3pa,idaf)
        call dscal(l2,0.5d0,q(idens),1)
        if (op2zora) then
           write(iwr,*) 'Initial electron density'
           call prtri(q(idens),num)
        end if
        call scalarf_cs(q,q(idens),q(idens),q(idnew),l1)
c
c.... open

      else
c     ----- occupation numbers -----
c
c     -oa- at q(i80)
c     -ob- at q(i81)
c
      do i = 1,l1
         popa=dzero
         popb=dzero
         if(i.le.na)popa=done
         q(i-1+i80)=popa
         if(i.le.nb)popb=done
         q(i-1+i81) = popb
      end do
     
c
cc     ----- compute initial guess density matrices -----
c     -da- in q(ida)
c     -db- in q(idb)
c
c ----- ensure vector orthogonality
      ipass=0
      nspin=na
      i=ibl3qa
      iocc=i80
      iblkp=ibl3pa
      if (ibl3pb.eq.0) then
         ibl3pb = ibl3pa
         write(iwr,*) ' CLOSED UHF simulation JVL2007'
      end if
 801  call rdedx(q(ivt),l3,i,idaf)
      call rdedx(q(i50),l2,ibl7st,num8)
c
      if (op2zora) then
         write(iwr,*) 'q(ivt) as read',i
         call prsq(q(ivt),l1,l1,l1)
         write(iwr,*) 'q(i50) as read'
         call prtri(q(i50),num)
      end if
c
      call mult2(q(ivt),q(i10),q(i50),l0,l0,l1)
c
      if (op2zora) then
         write(iwr,9008)
         write(iwr,*) 'q(i10), before orfog'
         call prtri(q(i10),num)
         write(iwr,*) 'q(i50)'
         call prtri(q(i50),num)
         write(iwr,*) 'q(ivt)'
         call prsq(q(ivt),l1,l1,l1)
      end if
      if(op2zora) call prsq(q(ivt),l1,l1,l1)

      if (nspin.ne.0) then
         call orfog(q(ivt),q(ivt),q(i10),q(i20),iky,ilifq,l0,l1,1)
         if (op2zora) then
            write(iwr,*) 'q(ivt),after orfog',i
            call prsq(q(ivt),l1,l1,l1) !may be filled only to l0 -XJ
         endif
      end if
c
      call wrt3(q(ivt),l3,i,idaf)
      call tdown(q(ivt),ilifq,q(ivt),ilifq,l0)
c
      if (l1.ge.maxorb) call caserr('l1 > maxorb, dmtx, scf_so')
c
      call dmtx(q(i10),q(ivt),q(iocc),iky,nspin,l1,l1)
c
      if (op2zora) then
         write(iwr,*) 'q(i10),dmtx'
         call prtri(q(i10),num)
      end if
c
      if(ipass.eq.1)go to 802
      i=ibl3qb
      iocc=i81
      ipass=ipass+1
      call wrt3(q(i10),l2,iblkp,idaf)
      iblkp=ibl3pb
      nspin=nb
      go to 801
 802  continue
      call dcopy(l2,q(i10),1,q(idb),1)
      call wrt3(q(i10),l2,iblkp,idaf)
      call rdedx(q(ida),l2,ibl3pa,idaf)
c
      call scalarf_cs(q,q(ida),q(idb),q(idnew),l1)
c
      if (op2zora) then
         write(iwr,*) 'density alpha'
         call prtri(q(ida),num)
         write(iwr,*) 'density beta'
         call prtri(q(idb),num)
      end if
c
      end if 
c
      end if
c
c ... idnew contains complex total density matrix 
c
      if(op2zora)then
        write(iwr,*)'q(idnew) at the start of iteration: ',iter
        call prevz(q,q(idnew),2*l1)
      endif

      if (type_z .eq. 'closed') then
         call extractd_z(q(idnew),q(idreal),l1,'aa')
c closed shells only one double coulomb
         do i=1,l2
            q(idens+i-1) = 2.0d0*q(idreal+i-1)
         enddo
         call hstar(q(idens),q(ijandk),q(isd1),1)
c
         if (op2zora) then
            write(iwr,*) 'fock matrix Jaa-Kaa q(ijandk)'
            call prtri(q(ijandk),num)
         end if
c
      else if (type_z .eq. 'open') then
         call extractd_z(q(idnew),q(idreal),l1,'aa')
         call extractd_z(q(idnew),q(idrealb),l1,'bb')
         call hstaru(q(idreal),q(idrealb),q(ijandk),q(ijandkb),1)
         call vadd (q(idreal),1,q(idrealb),1,q(idens),1,l2)
c
         if(op2zora) then
            write(iwr,*) 'density alpha real'
            call prtri(q(idreal),l1)
            write(iwr,*) 'density beta  real'
            call prtri(q(idrealb),l1)
            write(iwr,*) 'fock matrix Jaa-Kaa jandk'
            call prtri(q(ijandk),l1)
            write(iwr,*) 'fock matrix Jbb-Kbb jandkb'
            call prtri(q(ijandkb),l1)
            write(iwr,*) 'idens'
            call prtri(q(idens),l1)
         endif
c
         call symh(q(ijandk),q(i30),iky,0,0)
         call symh(q(ijandkb),q(i30),iky,0,0)
c
         if (op2zora) then
            write(iwr,*) 'jandk, after symh'
            call prtri(q(ijandk),l1)
            write(iwr,*) 'jandkb, after symh'
            call prtri(q(ijandkb),l1)
         end if
c
      else 
         call caserr ('unknown type_z type in scf_so')
      end if   
c
c build exchange contributions to spin orbit fock matrix
c
      call vclr(q(ieaai),1,l3)
      call vclr(q(ieba),1,l3)
      call vclr(q(iebai),1,l3)
      call extractd_z(q(idnew),q(idreal),l1,'aai')
      call hstar_z2(q(idreal),q(ieaai),1)
      call extractd_z(q(idnew),q(idreal),l1,'ba')
      call hstar_z2(q(idreal),q(ieba),1)
      call extractd_z(q(idnew),q(idreal),l1,'bai')
      call hstar_z2(q(idreal),q(iebai),1)
      if (type_z .eq. 'open') then
         call vclr(q(iebbi),1,l3)
         call extractd_z(q(idnew),q(idreal),l1,'bbi')
         call hstar_z2(q(idreal),q(iebbi),1)
      end if
c
      if (op2zora) then
         write(iwr,*) 'fock matrix eaai'
         call prsq(q(ieaai),num,num,num)
         write(iwr,*) 'fock matrix ebai'
         call prsq(q(iebai),num,num,num)
         write(iwr,*) 'fock matrix eba'
         call prsq(q(ieba),num,num,num)
         if (type_z .eq. 'open') then
           write(iwr,*) 'fock matrix ebbi'
           call prsq(q(iebbi),num,num,num)
         end if
      endif
c
      do i=1,num
         i1 = (i-1)*num
         do j=1,num
            i2 = i1 + j
            i3 = (j-1)*num + i
            eij = (q(ieaai+i2-1) - q(ieaai+i3-1))/2.0d0
            q(ieaai+i2-1) = eij
            q(ieaai+i3-1) = -eij
            eij = (q(ieba+i2-1) - q(ieba+i3-1))/2.0d0
            q(ieba+i2-1) = eij
            q(ieba+i3-1) = -eij
            eij = (q(iebai+i2-1) - q(iebai+i3-1))/2.0d0
            q(iebai+i2-1) = eij
            q(iebai+i3-1) = -eij
            if (type_z .eq. 'open') then
                eij = (q(iebbi+i2-1) - q(iebbi+i3-1))/2.0d0
                q(iebbi+i2-1) = eij
                q(iebbi+i3-1) = -eij
            end if 
         enddo
      enddo
c
c.... compute scalar zora corrections
c
      if ((mod(iter,niter_z).ne.0.and..not.ocz)
     1     .or. nat_z.ne.0) then
         call zora(q,q(idens),q(jblkf),'calc')
      else
         call zora(q,q(idens),q(jblkf),'force')
         write(iwr,*) 'forced recalculation of zora corrections'
      end if
c
c ... save scalar corrected 1e fock for total energy 
c
      if (op2zora) then
         write(iwr,*) ' zora corrected fock 1e'
         call prtri(q(jblkf),num)
         write(iwr,*) 'density: q(idens)'
         call prtri(q(idens),l1) 
      endif
c
c ... build scalar zora fock
c
      if (type_z .eq. 'closed') then
         call vadd(q(ijandk),1,q(jblkf),1,q(ijandk),1,l2)
      else if (type_z .eq. 'open') then
c
         if (op2zora) then
            write(iwr,*) 'ijandk,(fock matrix) voor vadd'
            call prtri(q(ijandk),l1)
            write(iwr,*) 'jblkf,(fock matrix) voor vadd'
            call prtri(q(jblkf),l1)
         end if
c
         call vadd(q(ijandk),1,q(jblkf),1,q(ijandk),1,l2)
c
         if (op2zora) then
            write(iwr,*) 'ijandk, na vadd'
            call prtri(q(ijandk),l1)
         end if
c
         call vadd(q(ijandkb),1,q(jblkf),1,q(ijandkb),1,l2)
      else
         call caserr ('invalid type_z type in scf_so')
      end if 
c
c...  de fock matrices voor beide spins worden in een 2*l1 bij 2*l1 
c...  complexe matrix q(ifock) gedaan. Using only 2*l0*2*l1!
c
      if (type_z .eq. 'closed') then
         call scalarf_cs(q,q(ijandk),q(ijandk),q(ifock),l1)
      else  if (type_z .eq. 'open') then
         call scalarf_cs(q,q(ijandk),q(ijandkb),q(ifock),l1)
      else
         call caserr('unknown type_z type')
      end if
c
      if(op2zora)then
         write(iwr,*) 'hermi check scalar fock'
         call hermi_check(q(ifock),2*l1)
      endif
c
c... compute spin orbit zora corrections
c
c...  add spin orbit parts to scalar zora fock matrix
c
      if (ospinbaan) then
         itert = iter
         if (ocz) itert = 999999999
         call so_zora(q,q(idens),.false.,itert)
         if (type_z .eq. 'closed') then
            call spin_orbitf_cs2(q,q(ifock),q(ieaai),q(ieaai),q(ieba),
     1                           q(iebai),mat2,rmat2,l1)
         else  if (type_z .eq. 'open') then
            call dscal(l3,-1.0d0,q(iebbi),1)
            call spin_orbitf_cs2(q,q(ifock),q(ieaai),q(iebbi),q(ieba),
     1                           q(iebai),mat2,rmat2,l1)
         else 
            call caserr ('unknown type_z type')
         end if
      else
         write(iwr,*) 'no spin orbit corrections'
      endif
c
c
c.... check if fock matrix is hermitian
c...  vies vies vies
c      call herm(q(ifock),2*l0,hm,itol)
c      call hermi_check(q(ifock),2*l1)
c     if (op2zora) then
c        print *,'hermiticity:',hm,itol
c     end if
c
      ihF = igmem_alloc(8*l1*l1)
c    
      ehf1 = etot_cs2(q(idnew),q(ifock),q(jblkf),
     +              rmat2,mat2,q(ihF),q,l1)
c
      call gmem_free(ihF)
c
c.... check if orthonormalising transformation is orthonormal
      if(op2zora) then
       write(iwr,*)'check orthonormalising transformation'
       call ortho_check_z(q,q(ivect),l1,l0)
      endif
c
c .... transform fock matrix onto basis of old orbitals
c
      iiscr = igmem_alloc(8*l3)
      call zgemm('c','n',2*l0,2*l1,2*l1,alpha,q(ivect),
     +        2*l1,q(ifock),2*l1,beta,q(iiscr),2*l0)
      call zgemm('n','n',2*l0,2*l0,2*l1,alpha,q(iiscr),
     +        2*l0,q(ivect),2*l1,beta,q(ifock),2*l0)
      call gmem_free(iiscr)
c
c
c	find largest off diagonal element
c       set criterium
c
	call lode(q(ifock),2*l0,ode)
	ocrit=ode.lt.acurcy
c
c
c .... level shifting
c
      if (lshift_z.eq.1 .and. iter .gt. 1) then
c
         if (dabs(de) .gt. 1.0d0) then
            vshift = 5.0d0*dsqrt(dabs(de))
         else if (dabs(de) .gt. 0.01) then
            vshift = 1.0d0*dsqrt(dabs(de))
         else if (dabs(de) .gt. 0.0000001) then
             vshift =0.1d0*dsqrt(dabs(de))
         else 
             vshift = 0.0d0
         end if 
         if (iter .gt. 50) then
            if (dabs(de) .gt. 0.02) then
               vshift = 2.0d0
            else 
               vshift = 0.1d0
            end if
         end if
         call zshift_cs(q(ifock),vshift,2*l0)
c
       else if (lshift_z.eq.0.and.iter.gt.1) then
c...      just one shift for now
           vshift = gapa1
           if (iter.gt.ibrk) vshift = gapa2
           call zshift_cs(q(ifock),vshift,2*l0)
       else if (lshift_z.eq.-1.and.iter.gt.1) then
c...       do nothing
       end if
c
c .... allocate space for diagonalisation routine
c
      lscr1 = 2*(4*l0-1 + l3)
      lscr2 = 6*l0 - 2
      lscratch = lscr1 + lscr2
c
      iscr1 = igmem_alloc (lscratch)
      iscr2 = iscr1 + lscr1
c
      if(op2zora)then
        write(iwr,*) 'hermicity of q(ifock) tested before zheev:'
        call hermi_check(q(ifock),2*l1)
      endif
      call zheev('v','l',2*l0,q(ifock),2*l0,q(iw1),q(iscr1),lscr1/2,
     +            q(iscr2),info)
c...  NOTE:q(ifock) now is a matrix with eigenvectors of the fock matrix
      if (op2zora) then
        if (info.eq.0) write(iwr,*) 'diagonalization fock succesfull'
        if (info.ne.0) write(iwr,*)'diagonalization fock failed'
      endif
c
      if(op2zora)then
         write(iwr,*) 'orthogonality of eigenvectors in q(ifock) tested'
         itemper = igmem_alloc(8*l3)
         call ortho_check_kolom(q(ifock),q(itemper),2*l1)
         call gmem_free(itemper)
      endif
      call gmem_free(iscr1)
c
      if (info .ne. 0) then
         write(iwr,*)'info:',info
         call caserr('error in scf_so_zora, zheev')
      endif 
c
c...  write eigenvalues
c
      if (op2zora) then
         write(iwr,665) iter
         do i=1,2*num
            write(iwr,670)i,q(iw1+i-1)
         enddo
      endif
c 
c .... transform eigenvectors back to original basis
c
      iiscr = igmem_alloc(8*l3)
c
      if (op2zora)then
         write(iwr,*)'q(ivect) orthonormalising transformation'
         call prevz(q,q(ivect),2*l1) 
      endif
c
      call zgemm('n','n',2*l1,2*l0,2*l0,alpha,q(ivect),
     +          2*l1,q(ifock),2*l0,beta,q(iiscr),2*l1)
      call zcopy(4*l3,q(iiscr),1,q(ivect),1)
c...  NOTE q(ivect) is now a matrix with the eigenvectors of the fock matrix
      if(op2zora)then
       write(iwr,*)'orthogonality of eigenvectors after transformation'
       call ortho_check_z(q,q(ivect),l1,l0)
      endif
      call gmem_free(iiscr)
c
c.... perform scaling of orbital energies
c
c.... compute scalar scaling factors
c
      if (oscalz) then
         call zora(q,q(idens),q(jblkf),'scale')
         if (ospinbaan) then
            call so_zora(q,q(idens),.true.,iter)
         else
              write(iwr,*) '** skip spinbaanscaling **'
         endif
c
c... allocate space for scaling factor matrix
c
         iscal = igmem_alloc(8*l3)
         call vclr(q(iscal),1,8*l3)
         call scale_so_cs(q,q(iscal),q(jblkf),q(ivect),q(iw1),
     +                 ospinbaan,oconvd,ocz) 
         call gmem_free(iscal) 
      end if 
c
      call etot_cs(q,ehf1,ode,vshift,iter,oconvd,ocz,ocrit)
c
c      print values  of
c      parameters for stopcriteria 
c
c     test for convergence
      iter = iter + 1
      iter2 =iter2+1
      nshift = nshift + 1
      if (.not. oconvd .and. iter .lt. maxcyc) then
          ocz = .false.
          go to 100
      else if (.not. ocz .and. iter .lt. maxcyc) then
          nshift = 0
          iter2= 0
          ocz = .true. 
          ospinbaan=.true.
          go to 100
      else if (ocz .and. iter .lt. maxcyc) then
      	  write(iwr,667) iter-1
      else if (iter .ge. maxcyc) then
          write(iwr,664)
      endif
c
c.... writing orbital energies and eigenvectors to file
c
      norb = 0
      j = 0
      write (iwr,668)
      do 90 i=1,2*num,1
        if (q(iw1+i-1)-q(iw1+i).lt.-0.0000001) then
            norb = norb + 1
            ndeg = i - j
            write(iwr,666)norb,q(iw1+i-1),ndeg
            j = i
        endif
90    continue
      write(iwr,669)
c
      write(iwr,663)
      call prevz(q,q(ivect),2*l1)
c
c.... in case there is only one iteration...
      if (iter.eq.2) then 
      write(iwr,*) ' ** there was only one it. - call zgemm again **'
          call zgemm('n','c',2*l1,2*l1,ne,(1.0d0,0.0d0),q(ivect),2*l1
     +               ,q(ivect),2*l1,(0.0d0,0.0d0),q(idnew),2*l1)
      end if
  
      write(iwr,672)
      itt  = igmem_alloc(l3)
      call densdecomplex(q(idnew),q(itt),l1)
      call prsq(q(itt),l1,l1,l1)
      call gmem_free(itt)
c
c....  analyse
c
      write(iwr,671)  
      call anal_so(q,q(ivect),q(idnew),q(iw1),l1)
c
      if (type_z .eq. 'open') call gmem_free_set(iebbi,i30)
      call gmem_free_set(ifock,iw1)   
c
      return
c
200   format(40x,'*****************'
     +      /40x,' ZORA SPIN-ORBIT'
     +      /40x,'*****************'///)
c
663   format(/40x,'  ------------'
     +       /40x,'<-eigenvectors->'
     +       /40x,'  ------------')
c
664   format('maximum number of iterations exceeded')
665   format('eigenvalues on',i3,' iteration')
666   format(i3,11x,f14.6,10x,i3)
667   format(/' energy converged after',i3,' iterations')
c
668   format(/
     +       ' ======================================================='/
     +       '   orbital              energy      #deg states        '/
     +       ' =======================================================')
669   format(' =======================================================')
670   format(i3,f22.8)
c
671   format(/'============================'
     +       ,'============================'/
     +   /40x,'*********************'/
     +    40x,' spin-orbit analyses '/
     +    40x,'*********************')
c
672   format(/40x,'  ----------------'
     +       /40x,'<-Electron density->'
     +       /40x,'  ----------------')
9008  format(1x)
      end
cc
c
      subroutine anal_so(q,vect,dens,ener,l1)
c
c...  analyse Spin-orbit 
c...  first version - jvl 2006
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/zorac)
INCLUDE(common/atmol3)
INCLUDE(common/dump3)
INCLUDE(common/scra7)
INCLUDE(common/common)
INCLUDE(common/iofile)
c
      COMPLEX vect(2*l1,2*l1),dens(2*l1,2*l1)
      REAL ener(2*l1)
c
      dimension q(*)
c
      l2 = l1*(l1+1)/2
      l3 = l1*l1
c   
      call secget(isect(497),19,ibl3pa)
      ibl3ea   = ibl3pa + lensec(l2)
      ibl3pb   = ibl3ea + lensec(l1)
      ibl3eb   = ibl3pb + lensec(l2)
      ibl3pc   = ibl3eb + lensec(l1)
      ibl3ec   = ibl3pc + lensec(l2)
      ibl7ew_z = ibl7la

c...  write MO eigenvalues
      itemp = igmem_alloc(l1)
       call dcopy(l1,ener(1),2,q(itemp),1)
       call wrt3(q(itemp),l1,ibl3ea,idaf)
       call dcopy(l1,ener(2),2,q(itemp),1)
       call wrt3(q(itemp),l1,ibl3eb,idaf)
      call gmem_free(itemp)
c
c...  generate density matrix
c
      call zgemm('no transpose','conjugate transpose',2*l1,2*l1,ne,
     +           (1.0d0,0.0d0),vect,2*l1,vect,2*l1,(0.0d0,0.0d0),
     +           dens,2*l1)
c
      if (op2zora) then
         write(iwr,*) 'eigenvector, end'
         call prevz(q,vect,2*l1)
         write(iwr,*)'check orthogonality eigenvectors'
         call ortho_check_z(q,vect,l1,l1)
         
         write(iwr,*) 'density, end'
         idensity = igmem_alloc(2*l3)
         call densdecomplex(dens,q(idensity),l1)
         call prevz(q,q(idensity),l1)
         call gmem_free(idensity)
      end if

c...  test, write energy weighted density
        ittest = igmem_alloc(2*2*l1*2*l1)
        iterst = igmem_alloc(l3)

        call dens_eweight(vect,ener,q(ittest),2*l1,ne)
        call densdecomplex(q(ittest),q(iterst),l1)

c...   triangularise and save
        call triangle(q(iterst),q(iterst),l1)
        call wrt3(q(iterst),l2,ibl7ew_z,num8)
        ibl7la = iposun(num8)
        call gmem_free(iterst)
        call gmem_free(ittest)
      kvo = igmem_alloc(l3)
      call rdedx(q(kvo),l3,ibl3qa,idaf)
c
      kd = igmem_alloc(l3)
      iprinv = 1
      if (opzora) iprinv = 10
c
c.....De dichtheden geschreven
c
      call extractd_z(dens,q(kd),l1,'aa') 
      call wrt3(q(kd),l2,ibl3pa,idaf)
      call sonat(q,q(kd),q(kvo),mouta,iprinv,'alpha')
      call extractd_z(dens,q(kd),l1,'bb')   
      call wrt3(q(kd),l2,ibl3pb,idaf)
      call sonat(q,q(kd),q(kvo),moutb,iprinv,'beta')
c
c       check imaginary
c
      dimax = 0.0d0
      call extractd_z(dens,q(kd),l1,'aai')   
      do i=1,l3
        dimax = dmax1(dimax,q(kd+i-1))
      end do
      call extractd_z(dens,q(kd),l1,'bbi')   
      do i=1,l3
        dimax = dmax1(dimax,q(kd+i-1))
      end do
      call extractd_z(dens,q(kd),l1,'bai')   
      do i=1,l3
        dimax = dmax1(dimax,q(kd+i-1))
      end do
c
c       check hermicity
c
      dherm = 0.0
      call extractd_z(dens,q(kd),l1,'ba')   
      do i=1,l1
         do j=1,i-1
            dherm = dmax1(dherm,dabs(q(kd-1+(i-1)*l1+j)-
     1                               q(kd-1+(j-1)*l1+i)))
         end do
      end do
c
c...  triangularise
c
      kk = kd - 1
      do i=1,l1
         do j=1,i
            kk = kk + 1
            q(kk) = q(kd-1+(i-1)*l1+j)
         end do
      end do
      call wrt3(q(kd),l2,ibl3pc,idaf)
      call sonat(q,q(kd),q(kvo),moutb+1,iprinv,'alfabeta')
c
      write(iwr,10) dimax,dherm
10    format(1x,/' $$$ max imag en non-herm ',1p,2e14.5)
c
      call gmem_free(kd)
      call gmem_free(kvo)
c
c...  analyse the occupied orbitals
c
      kd = igmem_alloc(l2)
      ks = igmem_alloc(l2)
      call rdedx(q(ks),l2,ibl7s,num8)
cjvl      call tranp(q(ks),q(ks))
c...    scale off-diagonal by 2.0 (cf QC-Utrecht reader)
      kk = ks -1
      do i=1,l1
         do j=1,i-1
            kk = kk + 1
            q(kk) = q(kk) * 2.0d0
         end do
         kk = kk + 1
      end do
c     print *,' Warning .... adapt not supported yet ; a bit iffy '
      write(iwr,*) ' '
      do imo = 1,ne
         call zgemm('n','c',2*l1,2*l1,1,(1.0d0,0.0d0),vect(1,imo),2*l1
     +             ,vect(1,imo),2*l1,(0.0d0,0.0d0),dens,2*l1)
         call extractd_z(dens,q(kd),l1,'aa') 
         ssa = ddot(l2,q(kd),1,q(ks),1)
         call extractd_z(dens,q(kd),l1,'bb')   
         ssb = ddot(l2,q(kd),1,q(ks),1)
         write(iwr,11) imo,ssa*100.0d0,ssb*100.0d0
11       format('  MO ',i5,4x,f6.2,' % alpha  ',f6.2,' % beta ')
      end do
c
      call gmem_free(ks)
      call gmem_free(kd)
c
      return
      end
**==sonat.f
      subroutine sonat(q,dens,vo,nsav,iprinv,text)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
c..    get natural orbitals from so densities 
c..    copied from uhfnat (anala.m)
c      density matrices need to be transformed inverse
c..
c..     ** called from anal_so **
c..      jvl Utrecht 2008
c
      dimension q(*),dens(*),vo(*)
      character*(*) text
c
INCLUDE(common/sizes)
c
INCLUDE(common/infoa)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
INCLUDE(common/dump3)
INCLUDE(common/runlab)
INCLUDE(common/unocas)
INCLUDE(common/zorac)
c
       zcom(5) = text
c..
c..    core partitioning
c..
       l1 = num
       l2 = num*(num+1)/2
       l3 = num*num
c
c      first compute core requirements
c
       i10 = igmem_alloc(2*l2)
       i20 = igmem_alloc(2*l2)
       i30 = igmem_alloc(2*l2)
       i50 = igmem_alloc(2*l1)
       i60 = igmem_alloc(l1)
c
c       i50 is scratch for minv (crayxmp) : 2*l1
c       i50,i60 scratch for minvrt (others) l1 each
c
c..   transform density matrix to orthonormal basis (UHF A-vectors)
c
      call dcopy(l3,vo,1,q(i20),1)
      call tdown(q(i20),ilifq,q(i20),ilifq,l1)
      call vclr(q(i50),1,2*l1)
c
c      invert transformation matrix t -> t(inv)
c
      d = 0.0d0
_IFN1(c)      call minvrt(q(i20),l1,d,q(i50),q(i60))
_IF1(c)      call minv(q(i20),l1,l1,q(i50),d,1.0d-30,0,1)
c
c      transpose it  t(inv) -> t(inv(dagger))
c
      call dagger(l1,l1,q(i20),l1,q(i30),l1)
c
c      then make t(inv).d.t(inv(dagger))
c      i30: t / i20: output d / i10: input d
c
      call mult2(q(i30),q(i20),dens,l1,l1,l1)
c
c..    diagonalise orthonormalised d-matrix (sort in decreasing order)
c..    result : occupations in i50 ; vectors in i30
c
      call gldiag(l1,l1,l1,q(i20),q(i60),q(i50),q(i30),iky,3)
c
c...  back-transform the eigenvectors (with symmetry adapted vectors)
c
      call dcopy(l3,vo,1,q(i10),1)
      call tfsqc(q(i10),q(i30),q(i20),l1,l1,l1)
c
c...   save ** mo's and eigenvalues = occupations ** on section nsav
c
      if (nsav.ne.0) then
c        if (iprinv.gt.0) then
            write(iwr,9028) text,nsav
c        endif
c         write(iwr,*)'spin-orbit val',oso
c saving i10 instead of i30; symmetry adapted non-backtransformed expected in scf -XJ
c         if(oso) !no longer needed
          call putq(zcom,ztitle,q(i50),q(i50),l1,l1,l1,1,1,
     *            q(i10),nsav,ibl3qa)
      end if
c
c..     transform back to original basis and print if requested
c
      if (iprinv.gt.1) then
         call tdown(q(i10),ilifq,q(i10),ilifq,l1)
         write(iwr,9030) text
         call prev(q(i10),q(i50),l1,l1,l1)
      endif
      if (iprinv.eq.1) then
         write(iwr,9032) text
         write(iwr,9040) (q(i50+i-1),i=1,l1)
         write(iwr,9040)
      end if
c
      call gmem_free_set(i10,i60)
c
      return
9028  format(/' ** ',a8,' SO natural orbitals saved in section'
     +       ,i4,' **')
9030  format(//30x,29('-')/30x,a8,' SO natural orbitals'/
     +       30x,29('-')/)
9032  format(10x,'----- ',a8,' SO natural orbital occupations -----')
9040  format(10x,7f14.7)
      end
      subroutine lode(fock,n,h)
c     finds largest off diagonal element of fock matrix
      COMPLEX fock(n,n)
      double precision a,h
      h=0.0
      a=0.0
      im = 0
      jm = 0
      do i=1,n
        do j=1,n
c          print *,'comm:info',i,j,fock(i,j)
           if (i.ne.j) then
           a=abs(fock(i,j))
           if (a.gt.h) then
               h=a
	       im = i
	       jm = j
           end if
          end if
         end do
      end do
c
      return
      end
c
      subroutine herm(fock,n,tf,itol)
c     checks if fock matrix is hermitian
      COMPLEX fock(n,n),df
      logical tf
      double precision dg
      integer itol
INCLUDE(common/iofile)
      df=(0.0,0.0)
      dg=0.0
      tf = .true.
      do i=1,n
         do j=i+1,n
            df=conjg(fock(i,j))-fock(j,i)
            dg=abs(df)
            if (dg.gt.1.0d-14 ) tf = .false. 
            if (dg.gt.1.0d-8 ) then
               itol = 1
               write(iwr,100) dg
            else if (dg.gt.1.0d-9 ) then
               itol = 2
            else if (dg.gt.1.0d-10) then
               itol = 3
            else itol = 0
            end if
          end do
       end do
       return
100    format(1x,'herm diff gt 1.0d-8, diff is ',d18.10)
       end
c**************************************************************
      subroutine hermi_check(fock,l)
c...  this subroutine checks the hermicity of a matrix
c...  a hermitian matrix should be equal to its conjugate
c...  transpose, so xresidu should be zero
      COMPLEX     fock,xresidu
      DIMENSION   fock(l,l)
      real        wresidu,iresidu,tresidu
      integer     i,j,l
      rresidu     =0.0
      iresidu     =0.0
      do i=1,l
       do j=i,l
       xresidu = fock(j,i)-conjg(fock(i,j))
       wresidu = wresidu + abs(real(xresidu))
       iresidu = iresidu + abs(imag(xresidu))
       enddo
      enddo
      tresidu = wresidu+iresidu
      if(tresidu.gt.1.0d-8)write(iwr,*)'hermicity in danger',wresidu,iresidu
      return 
      end
c****************************************************************
      subroutine ortho_check_z(q,vect,l1,l0)
c
c     checks for orthogonality by doing c*sc in which c is de 
c     eigenvector and s the overlap matrix. in case of 
c     orthogonality the resulting matrix is the identity.
c     the eigenvectors of a hermitian matrix are orthogonal
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/iofile)
INCLUDE(common/scra7)
      COMPLEX vect,sb,scratch,alpha,beta
      REAL s,q,tol,amreal,amimag,dia
c
      parameter(alpha = (1.0d0,0.0d0),beta = (0.0d0,0.0d0),
     +          tol1=1.0d-14,tol2=1.0d-9)
c
      dimension vect(2*l1,2*l1),q(*)
c      dimension vect(2*l1,2*l1),sb(2*l1,2*l1),s(*),scratch(*),q(*)
c
      is       = igmem_alloc(l1*(l1+1)/2)
      isb      = igmem_alloc(8*l1*l1)
      iscratch = igmem_alloc(8*l1*l1)
c
c     read triangular part overlap matrix
      call rdedx(q(is),l1*(l1+1)/2,ibl7s,num8)
c     fill overlap matrix with alpha and beta spin parts
      call scalarf_cs(q,q(is),q(is),q(isb),l1)
c     conjg(vect)*fock*vect should be identity if vect is orthogonal
      call zgemm('c','n',2*l1,2*l1,2*l1,alpha,vect,
     +        2*l1,q(isb),2*l1,beta,q(iscratch),2*l1)
      call zgemm('n','n',2*l1,2*l1,2*l1,alpha,q(iscratch),
     +        2*l1,vect,2*l1,beta,q(isb),2*l1)
      
c
      amreal = 0.0d0
      amimag = 0.0d0
      dia    = 0.0d0
      do i=1,2*l0
         do j=1,2*l0
            if (i .ne. j) then
               amreal = dmax1(amreal,dabs(complex2real(q(isb),l1,j,i,
     1                        'dreal')))
            else
               dia = dmax1(dia,dabs(complex2real(q(isb),l1,i,i,'dreal')
     1                     -1.0d0))
            end if
            amimag = dmax1(amimag,dabs(complex2real(q(isb),l1,j,i,
     1                     'dimag')))
         end do
      end do
c
      write(iwr,900)' orthocheck: amreal, dia, animag',amreal,dia
     +                  ,animag
c
      if (amreal+dia+amimag .gt. tol1) then
           write(iwr,*) 'orthonormalisation in trouble'
           write(iwr,*) 'max deviation of diagonal:',dia
           write(iwr,*) 'max deviation of off diagonal:',amreal
           write(iwr,*) 'max deviation of imaginary part:',amimag
          if (amreal+dia+amimag .gt. tol2) then
              write(iwr,*) 'max deviation gt tol2!!!!'
          end if
      end if
c
      call gmem_free(iscratch) 
      call gmem_free(isb)      
      call gmem_free(is)  
c
  900 format(A,X,3(E10.4,X),/)
      return
      end
      REAL function complex2real(sb,l1,j,i,mode)
c
      COMPLEX sb(l1,l1)
      character*(*) mode
      integer l1,j,i
c
      if (mode.eq.'dreal') complex2real = dreal(sb(j,i))
      if (mode.eq.'dimag') complex2real = dreal(sb(j,i))
c
      return 
      end
c**********************************************************************
c
      subroutine spin_orbitf_cs2(q,fock,eaai,ebbi,eba,ebai,
     +                           mat2,rmat2,l1)
c
c adds spin orbit matrices to original (scalar zora) fock matrix
c
c     s.f. 12/99
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/zorac)
c
c
      COMPLEX fock,mat2
      dimension q(*),fock(2*l1,2*l1),eaai(l1,l1),ebbi(l1,l1)
     +     ,eba(l1,l1),ebai(l1,l1),mat2(l1,l1),rmat2(l1,l1)
c
c ... now add spin orbit corrections to scalar zora fock matrix
c
      do i=1,l1
       do j=1,l1
         fock(j,i) = fock(j,i) + 
     +               dcmplx(0.0d0,rmat2(j,i)-eaai(j,i))
         fock(j+l1,i) = fock(j+l1,i) + dconjg(mat2(j,i)) +
     +                     dcmplx(eba(j,i),-ebai(j,i))
         fock(j,i+l1) = fock(j,i+l1) + mat2(i,j) + 
     +                  dcmplx(eba(i,j),ebai(i,j))
         fock(j+l1,i+l1) = fock(j+l1,i+l1)+
     +                     dcmplx(0.0d0,-rmat2(j,i)+ebbi(j,i))
       end do
      end do
c
      return
      end
c
c**********************************************************************
c 
        subroutine trace_z(q,tra,l1)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
c     implicit character *8 (z),character *1 (x)
c     implicit character *4 (y)
INCLUDE(common/zorac)
      dimension q(*) 
c calculate trace
         tra = 0.0d0
         do i=1,l1
           do j=1,l1
           if (i.eq.j) then
           itri = max(i,j)*(max(i,j)-1)/2 + min(i,j) 
           tra = tra + q(itri)
           else
           end if
           enddo
         enddo
      return
      end
c
c
      subroutine scalarf_cs(q,fa,fb,fock,l1)
c
c fock = (1.0d0,0.0d0)*f (scalar part of fock matrix)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/zorac)
INCLUDE(common/iofile)
c
      COMPLEX fock
c
      dimension fa(*),fb(*),fock(2*l1,2*l1),q(*)
c
      itri(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)  
c
       if (op2zora) then
          write(iwr,'(a)') ' == A Fock matrix '
          call prtri(fa,l1) 
          write(iwr,'(a)') ' == B Fock matrix '
          call prtri(fb,l1) 
       end if
c
      call vclr(fock,1,8*l1*l1)
      do i=1,l1
         do j=1,l1
           fock(j,i) = dcmplx(fa(itri(j,i)),0.0d0)
           fock(j+l1,i+l1) = dcmplx(fb(itri(j,i)),0.0d0)
         end do
      end do
c
       if (opzora) then
          write(iwr,'(a)') ' == scalarf_cs: total Fock matrix =='
          call prevz(q,fock,2*l1)
       end if
c
      return
      end 
c
c**********************************************************************
c
      subroutine extractd_z(dens,bai,l1,mode)
c
c extract real l1 x l1 block from total complex density matrix
c result in bai 
c
      implicit REAL   (a-h,p-z), integer(i-n), logical (o)
      COMPLEX dens 
      character*(*) mode
      dimension  bai(l1,l1),dens(2*l1,2*l1)
c
c
      if (mode .eq. 'aa' .or. mode .eq. 'bb') then
         kk = 0
         do i=1,l1
            do j=1,i
              kk = kk+1
              if (mode .eq. 'aa') bai(kk,1) = dreal(dens(j,i))
              if (mode .eq. 'bb') bai(kk,1) = dreal(dens(j+l1,i+l1))
            end do 
         end do
      else 
         do i=1,l1
            do j=1,l1
              if (mode .eq. 'aai')  bai(j,i) = dimag(dens(j,i))
              if (mode .eq. 'bbi')  bai(j,i) = dimag(dens(j+l1,i+l1))
              if (mode .eq. 'ba')  bai(j,i) = dreal(dens(i+l1,j))
c..           if (mode .eq. 'ba')  bai(j,i) = dreal(dens(j+l1,i))
c.. maar waarom dan niet op deze manier? anders is het toch getransponeerd?
              if (mode .eq. 'bai')  bai(j,i) = dimag(dens(j+l1,i))
            end do
         end do
      end if
c
      if (mode.ne.'aai' .and. mode.ne.'ba' .and. mode.ne.'bai' .and.
     +    mode.ne.'bb' .and. mode.ne.'bbi' .and. mode.ne.'aa')
     +          call caserr('invalid mode in extractd_z')
c
      return
      end
c
c**********************************************************************
c
      subroutine asymd_z(d,l1)
c
c antisymmetrize matrix d (square l1 x l1) (spinorbit zora densities)
c
      REAL q,d,pt5,th,a
      dimension d(l1,l1)
c
      data pt5,th /0.5d0,1.0d-06/ 
c
      do i=1,l1
         do j=1,l1
            a = (d(j,i) - d(i,j))*pt5
c           if (dabs(a) .lt. th) a = 0.0d0
            d(j,i) = a
            d(i,j) = -a
         end do
      end do
c 
      return
      end
c
c
c***********************************************************************
c
      subroutine so_zora(q,dens,oscale,iter)
c
c spinorbit zora corrections
c
      implicit REAL    (a-h,p-z), integer (i-n), logical (o)
      character*1 xn, xt
c 20/6/07      parameter (m100 =400000)
      parameter (m100 =900000)
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/mapper)
INCLUDE(common/zorac)
c
      dimension q(*),dens(*)
c
      data xn,xt/'n','t'/
c
      COMPLEX mat1,mat2
      common/somatrix1/mat1(m100),mat2(m100)
      common/somatrix2/rmat1(m100),rmat2(m100)
c
      c1 = cspeed
c
      if (.not.opre_zora) then
         call caserr('zora not prepared for spin-orbit correction')
      end if
c
      ofirst = iter.eq.0.or.iter.eq.999999999
      numold = num
c
      call swap_zora('internal')
      if (m100.lt.num*num) call caserr('dimensioning error in so')
c
      l1 = num
      l2 = num*(num+1)/2
      l3 = num*num
c      if (oscale) write(iwr,*) 'computing scaling factors'
c
c
c...  set pointers for remainder
      idens = igmem_alloc(l3)
      isd1  = igmem_alloc(l3)
      iout2 = igmem_alloc(l3)
      isd2  = igmem_alloc(l3)
      isd3  = igmem_alloc(l3)
      isd4  = igmem_alloc(l3)
      it1   = igmem_alloc(l3)
      ivc   = igmem_alloc(l3)
      is2   = igmem_alloc(l3)
      is4   = is2
      isx   = igmem_alloc(l3)
      isy   = igmem_alloc(l3)
      isz   = igmem_alloc(l3)
      is    = igmem_alloc(l3)

c...  the coulomb matrix was just calculated 
c...  (in scalar zora) so read it 
      call rdedx(q(isd1),nwv_z,ibcoul_z,num8)
      
c...  even for get atom option this is okay; however,
c...  extra work can be avoided since coulomb matrix is 
c...  constant sf
c
c compute vc
c     using a coulomb matrix at q(isd1)
c
      cc2=2.0d0*c1**2            
c
c compute v
c      
      call rdedx(q(ivc),nwv_z,iblv_z,num8)
      call rdedx(q(is),nwv_z,ibls_z,num8)
      if (.not. oscale) then
          do 198 i=1,l2
             q(ivc+i-1)=-(q(ivc+i-1)+q(isd1 +i-1))/cc2 + q(is+i-1)
198       continue
          if (op2zora) write(iwr,*) 'vc matrix, no scaling'
          if (op2zora) call prtri(q(ivc),num)
      else
          do 298 i=1,l2
             q(ivc+i-1)=(-q(ivc+i-1)-q(isd1 + i-1) + q(is+i-1)*cc2)/c1
298       continue
          if (op2zora) write(iwr,*) 'vc matrix, with scaling'
          if (op2zora) call prtri(q(ivc),num)
      endif
c
c transform to orthogonal basis:
c transform matrix
c 
      call rdedx(q(is2),nwmat_z,ibsmin_z,num8)
      call mult2(q(is2),q(isd2),q(ivc),num,num,num)
c
c compute zora corrections
c
c
      call square(q(ivc),q(isd2),num,num)
c
      d=0.d0
      call minvrt(q(ivc),num,d,q(isd2),q(isd1))
c
c potential independent, need to do this only once
c
c real part
c
         call rdedx(q(isx),nwmat_z,ibsx_z,num8)
         call rdedx(q(isy),nwmat_z,ibsy_z,num8)
         call rdedx(q(isz),nwmat_z,ibsz_z,num8)
c
             call mxmaa(q(isz),num,1,q(isx),num,1,q(isd1),num,1,
     +         num,num,num)
             call mxmaa(q(isx),num,1,q(isz),num,1,q(isd2),num,1,
     +         num,num,num)
             n=0
             do i=1,num
                do j=1,num
                   n=n+1
                   q(isd1+n-1)=q(isd1+n-1)-q(isd2+n-1)
                enddo
             enddo
c
c imaginary part
c
             call mxmaa(q(isz),num,1,q(isy),num,1,q(isd2),num,1,
     +         num,num,num)
             call mxmaa(q(isy),num,1,q(isz),num,1,q(isd3),num,1,
     +         num,num,num)
             n=0
             do i=1,num
                do j=1,num
                   n=n+1
                   q(isd2+n-1)=q(isd2+n-1)-q(isd3+n-1)
                enddo
             enddo
             n=0
c
             n=0
             do i=1,num
                do j=1,num
                   n=n+1
                   mat1(n)=q(isd1+n-1)+(0.0d0,-1.0d0)*q(isd2+n-1)
                enddo
             enddo
c
c compute (imaginary) SO-z (pot. indep.) part
c
             call mxmaa(q(isx),num,1,q(isy),num,1,q(isd1),num,1,
     +          num,num,num)
             call mxmaa(q(isy),num,1,q(isx),num,1,q(isd2),num,1,
     +          num,num,num)
             n=0
             do i=1,num
                do j=1,num
                   n=n+1
                   rmat1(n)=q(isd1+n-1)-q(isd2+n-1)
                enddo
             enddo
c
c sf end of 'ofirst' if block
c
c
c compute off diagonal spinorbit corrections to fock matrix
c (potential dependent part)
c
         if (oscale) then
             call mxmaa(q(ivc),num,1,q(ivc),num,1,q(isd1),num,1,
     +                  num,num,num)
             call dcopy(num*num,q(isd1),1,q(ivc),1)
         endif

         call mxmaa(q(isz),num,1,q(ivc),num,1,q(isd1),num,1,
     +         num,num,num)
         call mxmaa(q(isd1),num,1,q(isx),num,1,q(isd2),num,1,
     +         num,num,num)
         call mxmaa(q(isx),num,1,q(ivc),num,1,q(isd1),num,1,
     +         num,num,num)
         call mxmaa(q(isd1),num,1,q(isz),num,1,q(isd3),num,1,
     +         num,num,num)
         n=0
         do i=1,num
            do j=1,num
               n=n+1
               q(isd1+n-1)=q(isd2+n-1)-q(isd3+n-1)
c volgorde
            enddo
         enddo
c
c imaginary part
c
         call mxmaa(q(isy),num,1,q(ivc),num,1,q(isd2),num,1,
     +         num,num,num)
         call mxmaa(q(isd2),num,1,q(isz),num,1,q(isd3),num,1,
     +         num,num,num)
         call mxmaa(q(isz),num,1,q(ivc),num,1,q(isd2),num,1,
     +         num,num,num)
         call mxmaa(q(isd2),num,1,q(isy),num,1,q(isd4),num,1,
     +         num,num,num)
         n=0
         do i=1,num
            do j=1,num
               n=n+1
               q(isd2+n-1)=q(isd3+n-1)-q(isd4+n-1)
            enddo
         enddo
c
         pt5=0.5d0
         if (.not. oscale) then
            n=0
            do i=1,num
               do j=1,num
                  n=n+1
                  mat2(n)=(q(isd1+n-1)+(0.0d0,1.0d0)*q(isd2+n-1)-
     +                 mat1(n))*pt5
               enddo
            enddo
         else
c    c: was commentaar; rest voorlopig
            do i=1,num*num
c    c              mat2(i)=(sdummy1(i)+(0.0d0,1.0d0)*sdummy2(i))*pt5
               mat2(i)=(q(isd1+i-1)+(0.0d0,1.0d0)*q(isd2+i-1))
            enddo
         endif
c
c
c compute potential dependent so-z part
c
         call mxmaa(q(isx),num,1,q(ivc),num,1,q(isd1),num,1,
     +       num,num,num)
         call mxmaa(q(isd1),num,1,q(isy),num,1,q(isd2),num,1,
     +       num,num,num)
c
         call mxmaa(q(isy),num,1,q(ivc),num,1,q(isd1),num,1,
     +       num,num,num)
         call mxmaa(q(isd1),num,1,q(isx),num,1,q(isd3),num,1,
     +       num,num,num)
         n=0
         if (.not. oscale) then
            do i=1,num
               do j=1,num
                  n=n+1
c                 rmat2(n)=sdummy2(n)-sdummy3(n)
                  rmat2(n)=pt5*(q(isd2+n-1)-q(isd3+n-1)-rmat1(n))
               enddo
            enddo
         else
c    c: was commentaar; rest voorlopig
            do i=1,num*num
c    c               rmat2(i)=pt5*(sdummy2(i)-sdummy3(i))
                rmat2(i) = q(isd2+i-1)-q(isd3+i-1)
             enddo
         endif
c
c  ... transform back to origininal basis
c
c         if (oscale .and. op2zora) then
c            print *,'rmat before transformations:'
c            call prsq(rmat2,num,num,num)
c         endif
c
         call rdedx(q(is4),num*num,ibsplus_z,num8)
c
         call mxmaa(q(is4),num,1,rmat2,num,1,q(isd1),num,1,
     +       num,num,num)
         call mxmaa(q(isd1),num,1,q(is4),num,1,rmat2,num,1,
     +       num,num,num)
   
      call rdedx(q(iout2),nwtrat_z,ibtrout_z,num8)
c
c transform back to original external basis
c
         call dgemm(xn,xn,num,numold,num
     +    ,1.0d0,rmat2
     +    ,num,q(iout2),num,0.0d0,q(isd1),num)
c
         call dgemm(xt,xn,numold,numold,num
     +    ,1.0d0,q(iout2)
     +    ,num,q(isd1),num,0.0d0,rmat2,numold)
c
c    
c follow same procedure for mat2 (2 times)
c
c real part first
c
         do i=1,num*num
            q(isd1+i-1) = dreal(mat2(i))
         enddo
c          
         call mxmaa(q(is4),num,1,q(isd1),num,1,q(isd2),num,1,
     +       num,num,num)
         call mxmaa(q(isd2),num,1,q(is4),num,1,q(isd1),num,1,
     +       num,num,num)
c
c transform back to original external basis
c
         call dgemm(xn,xn,num,numold,num
     +    ,1.0d0,q(isd1)
     +   ,num,q(iout2),num,0.0d0,q(isd2),num)
c
         call dgemm(xt,xn,numold,numold,num
     +    ,1.0d0,q(iout2)
     +    ,num,q(isd2),num,0.0d0,q(isd1),numold)
c
c imaginary part
c
         do i=1,num*num
            q(isd2+i-1) = dimag(mat2(i))
         enddo
c         
         do i=1,numold*numold
            mat2(i) = (1.0d0,0.0d0)*q(isd1+i-1)
         enddo 
c
         call mxmaa(q(is4),num,1,q(isd2),num,1,q(isd1),num,1,
     +       num,num,num)
         call mxmaa(q(isd1),num,1,q(is4),num,1,q(isd2),num,1,
     +       num,num,num)
c
c transform back to original external basis
c
         call dgemm(xn,xn,num,numold,num
     +    ,1.0d0,q(isd2)
     +    ,num,q(iout2),num,0.0d0,q(isd1),num)
c
         call dgemm(xt,xn,numold,numold,num
     +    ,1.0d0,q(iout2)
     +    ,num,q(isd1),num,0.0d0,q(isd2),numold)
c
         do i=1,numold*numold
            mat2(i) = mat2(i) + (0.0d0,1.0d0)*q(isd2+i-1)
         enddo
c
c     num = numold
c
      call swap_zora('external')

      call gmem_free(is)
      call gmem_free(isz)
      call gmem_free(isy)
      call gmem_free(isx)
      call gmem_free(is2)
      call gmem_free(ivc)
      call gmem_free(it1)
      call gmem_free(isd4)
      call gmem_free(isd3)
      call gmem_free(isd2)
      call gmem_free(iout2)
      call gmem_free(isd1)
      call gmem_free(idens)
c
      return
9088  format (/20x,23('*')/20x,'symmetrized matrix'/20x,23(
     +     '-'))
5030  format (2i4,3f16.9)
      end
c
c***********************************************************************
c
      subroutine etot_cs(q,ehf1,tester,vshift,iter_z,oconvd,ocz,ocrit)
c
c... calculates total energy for unscaled zora spin orbit (closed shell)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/zorac)
INCLUDE(common/scfopt)
INCLUDE(common/funct)
INCLUDE(common/restri)
INCLUDE(common/scra7)
      common/ezora/eunscaled,escaled,sum,sum2,sum2e
      COMPLEX sum2e
      logical ocrit
      dimension q(*)
      data m1,m10,m16/1,10,16/
c
c
      save e_unsc_old
c
c first build total density matrix
c
      en = enucf(nat,czan,c)
c
      e_unsc = ehf1 
      e_scaled = ehf1 - eunscaled + escaled 
      de = e_unsc - e_unsc_old
c      oconvd = (dabs(de) .lt. 1.0d-5) 
c
      oconvd = (tester .lt. acurcy .and. (iter_z .gt. 1)
     +         .and.ocrit)
c
      if (iter.eq.1.and..not.op2zora) write(iwr,100)
      if(op2zora) write(iwr,100)
      if (iter.lt.10.or.10*(iter/10).eq.iter) then !only show 1 in 10
        write(iwr,101) iter,e_scaled + en, e_scaled,de,tester,vshift
      endif   !avoids overloading the screen -XJ
c...  this is written in the same way as a normal SCF
c...  so the standard output can be used in optim.m:4119
      enrgy = e_scaled + en
      i10 = igmem_alloc(m10)     
      q(1-1+i10) = en
      q(2-1+i10) = e_scaled
      q(3-1+i10) = e_scaled + en
      do i=4,10
       q(i-1+i10) = 0.
      enddo
      call secput(isect(494),16,1,iblk16)
      call wrt3(q(i10),m10,iblk16,idaf)
      call gmem_free(i10)      

c      if (oconvd.and.ocz) then
c         write(iwr,90)en
c         write(iwr,91)e_scaled+en
c         write(iwr,92)e_unsc+en
c         write(iwr,93)e_scaled
c         write(iwr,94)e_unsc
c         write(iwr,95)
c      end if
      e_unsc_old = e_unsc
c
      return
5     format (a17,f14.7)
6     format ('sum of imag part',f14.5)
7     format ('total',f25.7)
8     format (5x,i4,2f22.7,2(1pe10.2),0pf10.3)
c
100   format(1x,'======================================================'
     +         ,'===================================='
     +      /1x,'  cycle            total         electronic           '
     +         ,' e conv.         tester    virtual'
     +      /1x,'                  energy             energy           '
     +         ,'                             shift'
     +      /1x,'======================================================'
     +         ,'====================================')
101   format(1x,i5,f19.8,f19.8,f19.8,2x,f13.8,5x,f6.3)
c
90    format(/15x,'****************************',22('*')
     +       /15x,'*nuclear energy'            ,16x,f18.6,'*')
91    format (15x,'*scaled total energy'       ,11x,f18.6,'*')
92    format (15x,'*unscaled total energy'     ,9x, f18.6,'*')
93    format (15x,'*scaled electronic energy'  ,6x, f18.6,'*')
94    format (15x,'*unscaled electronic energy',4x, f18.6,'*')
95    format (15x,'****************************',22('*'))
c
      end
c 
c**********************************************************************
c
      REAL function etot_cs2(dens,fock,f_s,soaai,sobac,hF,q,l1)
c
c    calculate total HF energy for spin orbit zora 0.5*P(h+F)
c
c    hF : space for total 1e (t+v) fock
c    fock  : total fock matrix (input)
c    f_s   : scalar corrected 1e fock (l2)
c    dens  : total density matrix (complex)
c    soaai : spinorbit zora correction (imaginary aa part)
c    sobac : complex so zora correction (alpha beta part)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(common/zorac)
INCLUDE(common/iofile)
      common/ezora/eunscaled,escaled,sum,sum2
      COMPLEX fock,hF,sobac,ehf,zdotc,dens
      dimension fock(2*l1,2*l1),f_s(*),hF(2*l1,2*l1),q(*)
     +          ,soaai(l1,l1),sobac(l1,l1),dens(2*l1,2*l1)
c
      itri(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)  
c
      do i=1,l1
         do j=1,l1
            hF(j,i) = fock(j,i) + 
     +                dcmplx(f_s(itri(j,i)),soaai(j,i))  
            hF(j+l1,i+l1) = fock(j+l1,i+l1) + 
     +                      dcmplx(f_s(itri(j,i)),-soaai(j,i))
            hF(j+l1,i) = fock(j+l1,i) + dconjg(sobac(j,i))
c           hF(j,i+l1) = dconjg(fock(j+l1,i)) + sobac(j,i)
            hF(j,i+l1) = fock(j,i+l1) + sobac(i,j)
         end do
      end do
c
      if (op2zora) then
c          write(iwr,*) 'etot_cs2 prints h+f',l1
         write(iwr,*) 'etot_cs2 checks hermicity q(ifock)'
         call hermi_check(q(ifock),2*l0)         
c         write(iwr,*)'hermicity:',ohm
c         print*,itol
c          call prevz(q,hf,l1)
c          write(iwr,*) 'etot_cs2 prints density matrix'
c         call herm(q(ifock),2*l0,ohm,itol)         
c         write(iwr,*)' hermiticity:',ohm
c         print*,itol
c          call prevz(q,dens,2*l1)
      end if 
c
      ehf = zdotc(4*l1*l1,dens,1,hF,1)
c
      if (dabs(dimag(ehf)/4.0d0) .gt. 10.0d-15) 
     +    write(iwr,*) '*** WARNING *** imaginary part of energy',
     +                 dimag(ehf)/2.0d0
      etot_cs2 = dreal(ehf)/2.0d0
      return
      end
c
c**********************************************************************
c 
      subroutine totals2(s2,s22,l0,l1)
c
c  build total (2*l1 x 2*l1) orthonormalising transformation from s2
c
      implicit REAL  (a-h,p-z),integer(i-n),logical(o)
c 
      COMPLEX s22
      dimension s2(l1,l1),s22(2*l1,2*l1)
c clear matrix with vclr, all elements set to zero 
      call vclr(s22,1,2*2*l1*2*l1)
c        ... dcmplx because cmplx goes through real
      do i=1,l0
         do j=1,l1
           s22(j,i) = dcmplx(s2(j,i))
           s22(j+l1,i+l0) = dcmplx(s2(j,i))
         end do
      end do
c
      return
      end
c
c**********************************************************************
c
      subroutine prevz(q,ev,num)
c
c prints complex matrices from spin orbit zora calculation
c
      REAL q
      COMPLEX  ev
      integer i,j,num,kk
      dimension q(*),ev(num,num)
INCLUDE(common/iofile)
c
      id1 =igmem_alloc (num*num)
c
      kk = id1
      do i=1,num
       do j=1,num
         q(kk)=dreal(ev(j,i))
         kk = kk+1
       enddo
      enddo    
c 
      write(iwr,*)'real part'
      call prsq(q(id1),num,num,num)
c
      kk = id1
      do i=1,num
       do j=1,num
         q(kk)=dimag(ev(j,i))
         kk = kk+1
       enddo
      enddo   
c
      write(iwr,*) 'imaginary part'
      call prsq(q(id1),num,num,num) 
c
      call gmem_free(id1)
      return
      end
c

c**************************************************************************
c
      subroutine hstar_z2(d,exch,nopk)
c
c computes anti symmetric exchange matrix for spin orbit zora
c calculation (based on hstar); uses exch_zora.
c
c     ----- ***********************************  ------
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
      common/filel/nfile,notape(20),iblk(20),lblk(20)
      common/disc/isel(5),ipos(maxlfn)
      common/iofile/ ir,iw,ip,main,iblkm,idaf,iblkd,num8,iblk21,
     * iblkqa,iblkpa,iblkea,iblkqb,iblkpb,iblkeb,iblkg,iblkhs,
     * iblksp(3),nav
      dimension d(*),exch(*)
      common/infoa/nat,ich,mul,num,nx,ne,na,nb
      common/blkin/gin(510),nint
      common/mapper/ia(maxorb),ikyp(1)
INCLUDE(common/sortp)
c
INCLUDE(common/maxlen) 
      common/vinteg/q(1)
      dimension mww(2)
      equivalence (mww(1),dumm)
c
      data pt5,two/0.5d0,2.0d0/
      nzora=num*num
c
c     l2 = (num*(num+1))/2
c
c     m=0
c     do i=1,num
c        do j=1,num
c           m=m+1
c           if (i.eq.j) d(m)=d(m)/2.0d0
c        enddo
c     enddo
c
c
      if (nopk.ne.1) then
         call caserr('integrals in supermatrix form')
         return
      endif 
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c 
      n2=num*num
      call vclr(exch,1,n2)
      do 1000 ifile=1,nfile
         main=notape(ifile)
         if (omem(main)) then
c
c original code: call to sgmtmm
c
                call caserr('Error in hstar_z')
c
         else
             call search(iblk(ifile),main)
             call find(main)
             jblock=lblk(ifile)
1004         jblock=jblock+1
             call get(gin,nw)
             if (nw .eq. 0) go to 1005
             if (jblock .ne. 0)  call find(main)
c
             call exch_zora(d,exch,num)
c
             if (jblock) 1004,1000,1004
1005         lblk(ifile)=iblk(ifile)-ipos(main)+1
         endif
1000  continue
c
c     m=0
c     do i=1,num
c        do j=1,num
c           m=m+1
c           if (i.eq.j) d(m)=d(m)+d(m)
c        enddo
c     enddo
c
c vgl. originele dichtheid (aa)
c
c     print *,'hstar:'
c     call prev_z(exch,num,num,num)
      pt25 = 0.25d0
      call dscal(nzora,pt5,exch,1)
c     call prev_z(exch,num,num,num)
c
c
      return
      end
c
c*********************************************************************** 
c
      subroutine exch_zora(p,exch,num)
c
c build exchange matrix from square density matrix
c
c p: density (input); exch: exchange (output) 
c 
      implicit REAL (a-h,o-z)
c
      dimension p(*),exch(*)
c
INCLUDE(common/atmblk)
      common/blkin/gg(510),mword
      common/craypk/integ(1)
c
      call unpack(gg(num2e+1),lab816,integ,numlab)
c
c     write(*,2)
2     format('i',4x,'j',4x,'k',4x,'l',4x,'integral')
c
      iword = 1
      do 100 iw=1,mword
         i=integ(iword)
         j=integ(iword+1)
         k=integ(iword+2)
         l=integ(iword+3)
         gik=gg(iw)
         g2=gik+gik
         gil=gik
c        if (i .eq. k .or. j .eq. l) gik=g2
c        if (j .eq. k) gil=g2
c        write(*,1)i,j,k,l,gik
c
         il = (l-1)*num + i
         kj = (j-1)*num + k
         jl = (l-1)*num + j
         ki = (i-1)*num + k
         ik = (k-1)*num + i
         lj = (j-1)*num + l
         li = (i-1)*num + l
         jk = (k-1)*num + j   

c        write(*,3)il,p(il),kj,p(kj),jl,p(jl),ki,p(ki)
c        write(*,3)ik,p(ik),lj,p(lj),li,p(li),jk,p(jk)
3        format(4(i2,f8.4,2x))

c
         if ((i.eq.j).and.(j.eq.k).and.(k.eq.l)) then
c
             exch(il) = gil*p(jk) + exch(il)
c
         else if ((j.eq.k).and.(k.eq.l)) then
c
             eil = gil*p(jk) + exch(il)              !zero
             ejl = gik*p(ik) + exch(jl)
             eki = gik*p(lj) + exch(ki)              !zero
             ekj = gil*p(li) + exch(kj)
c
             exch(il)=eil
             exch(jl)=ejl
             exch(kj)=ekj
             exch(ki)=eki
c
         else if ((i.eq.j).and.(j.eq.k)) then
c
             eil = gil*p(jk) + exch(il)
             eik = gik*p(jl) + exch(ik)
             ekj = gil*p(li) + exch(kj)
             elj = gik*p(ki) + exch(lj)
c
             exch(il)=eil
             exch(ik)=eik
             exch(kj)=ekj
             exch(lj)=elj
c
         else if ((i.eq.j) .and. (k.eq.l)) then
c
             eil = gil*p(jk) + exch(il)
             ekj = gil*p(li) + exch(kj)
c
             exch(il)=eil
             exch(kj)=ekj
c
         else if ((i.eq.k).and.(j.eq.l)) then
c
             eil = gil*p(jk) + exch(il) 
             ejl = gik*p(ik) + exch(jl)              !zero
             eik = gik*p(jl) + exch(ik)
             ejk = gil*p(il) + exch(jk)
c
             exch(il)=eil
             exch(jl)=ejl
             exch(ik)=eik
             exch(jk)=ejk
c
         else if (i.eq.j) then
c
             eil = gil*p(jk) + exch(il)
             eik = gik*p(jl) + exch(ik)
             ekj = gil*p(li) + exch(kj)
             elj = gik*p(ki) + exch(lj)
c
             exch(il)=eil
             exch(ik)=eik
             exch(kj)=ekj
             exch(lj)=elj
c
         else if (k.eq.l) then
c
             eil = gil*p(jk) + exch(il)              !zero
             ejl = gik*p(ik) + exch(jl)
             eki = gik*p(lj) + exch(ki)              !zero
             ekj = gil*p(li) + exch(kj)
c
             exch(il)=eil
             exch(jl)=ejl
             exch(kj)=ekj
             exch(ki)=eki
c
         else
c            print *,'different?',i,j,k,l
             eil = gil*p(jk) + exch(il)
             ejl = gik*p(ik) + exch(jl)
             eik = gik*p(jl) + exch(ik)
             ejk = gil*p(il) + exch(jk)
c
             ekj = gil*p(li) + exch(kj)
             eki = gik*p(lj) + exch(ki)
             elj = gik*p(ki) + exch(lj)
             eli=  gil*p(kj) + exch(li)
c
             exch(il)=eil
             exch(jl)=ejl
             exch(ik)=eik
             exch(jk)=ejk 
             exch(kj)=ekj
             exch(ki)=eki
             exch(lj)=elj
             exch(li)=eli
c
        endif
c        
1        format(4(i2,2x),f8.4)
         iword=iword+4
100   continue
c
      return
c 
      end
c
c**********************************************************************
c
      subroutine occ_nr
c
c... sets up vector of occupation numbers for spin orbit zora calculation
c... as it us unrestricted nd runs over ball orbitals
c
      implicit REAL (a-h,p-z), integer(i-n), logical(o)
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/zorac)
INCLUDE(common/occnr_z)
INCLUDE(common/iofile)
c
c     if (op2zora) print *,xoccav_z,nd_z,no_z

      do i=1,2*num
         xocc_z(i)=0.0d0
      end do
c 
      if (2*num.gt.nocc_z) call  caserr('xocc_z overflow')
c
      nd_z = ne
      no_z =0
c
      do i=1,nd_z
         xocc_z(i)=1.0d0
      enddo
      if (no_z .gt. 0) then
         call caserr('occ_nr not reaqdy')
         do i=nd_z+1,nd_z+no_z
            xocc_z(i) = xoccav_z 
         enddo
      end if
c
      sum = 0.0d0
      do i=1,2*num
          sum = sum +xocc_z(i)
      enddo
      if (opzora) then
          write(iwr,*)' occupation numbers for kramers degenerate pairs'
          do i=1,2*num
            write(iwr,667) i,xocc_z(i)
          enddo
      end if
      if (dabs(dfloat(ne)-sum) .gt. 1.d-05)
     +        call caserr('confusion in spin orbit occupation')
667   format(i3,f14.5)
      return
      end

c
c***********************************************************************
c
      subroutine dens_z(q,dens)
c
c...  generate alien density matrix for zora (is_z > 0)
c...  or read it from section -is_z    (is_z < 0)
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension q(*),dens(*)
c
INCLUDE(common/sizes)
INCLUDE(common/zorac)
INCLUDE(common/tran)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/machin)
c
      common/atmol3/mina(2),mouta,moutb
      common/junk/ilifc2(maxorb),nterm2(maxorb),
     + it2(mxorb3),ct2(mxorb3),otran2
c
      if (is_z.lt.0) then
c
c...     read density matrix directly
c...     only working for cas yet (see anal routine dipmat)
c
         call secget(-is_z,990,ibl)
         call rdedx(dens,num*(num+1)/2,ibl,idaf)
         write(iwr,100) -is_z,ibl,idaf
100      format(' read density matrix from section',i4,' block',i6,
     1          ' of unit',i3)
         return
      end if
c....
c
c...  note is_z may not be the current output sectors, due to possible
c...  overwriting of the occupations
c
       if (is_z.eq.mouta.or.is_z.eq.moutb) 
     1 call caserr('zora vectors section may not be the current output')
      l1 = num
      l2 = num*(num+1)/2
      l3 = num*num
c
      i10 = igmem_alloc(l3+2*l1)
      i20 = i10 + l3
      i30 = i20 + l1
c
      call getq(q(i10),q(i20),q(i30),nbas,newb,3,ieig,ipop,is_z,'zora')
      if (nbas.ne.num.or.newb.ne.num) call caserr('bad vectors dens_z')
      if (ipop.ne.1)call caserr('vectors without occupations in dens_z')
c
c...  print the occupations
c
      n = 0
      nn = 0
      do i=1,l1
         if (q(i30+i-1).ne.2.0d0) go to 1
         n = n + 1
      end do
1     do i=l1,1,-1
         if (q(i30+i-1).ne.0.0d0) go to 2
         nn = nn + 1
      end do
2     write(iwr,6) n,l1-n-nn,nn
      if (l1-n-nn.gt.0) write(iwr,7) (q(i30+i-1),i=n+1,l1-nn)
6     format(' # Doubly occupied',i5,' # Variably occupied',i5,
     1       ' # Empty orbitals',i5)
7     format(' Occupations ',5f12.8,/,(13x,5f12.8))
c
c...   go to non-adapted
c
      call secget(is_z,3,ioctr2)
      nav = lenwrd()
      ioctr2 = ioctr2 + lensec(2*maxorb+1+6/nav) + 1
      call readi(ilifc2,mach(9)*nav,ioctr2,idaf)
      call tdown2(q(i10),l1,l1,ilifc2,nterm2,it2,ct2,otran2)
c
      call vclr(dens,1,l2)
      do i=1,l1
         pop = q(i30+i-1)
         ii = i10 + (i-1)*l1 - 1
         if (dabs(pop).gt.1.0d-13) then
            kl = 0
            do k=1,l1
               do l=1,k
                  kl = kl + 1
                  dens(kl) = dens(kl) + pop*q(ii+k)*q(ii+l)
               end do
            end do
         end if
      end do
c
      call gmem_free(i10)
c
      return
      end
c
c*********************************************************************
c
      subroutine zshift_cs(fock,vshift,ndim)
c
c shift virtual levels of closed shell spin orbit ZORA  fock 
c matrix up by vshift 
c
INCLUDE(common/occnr_z)
c
      COMPLEX fock
      REAL vshift,test
      integer ndim
c
      dimension fock(*)
c
      do i=1,ndim
         ii = (i-1)*ndim +i
         if (dabs(xocc_z(i)) .lt. 1.0d-5) then
            fock(ii) = fock(ii) + vshift*(1.0d0,0.0d0)
         end if
       end do
c
       return
       end
c
c**************************************************************
c
      subroutine scale_so_cs(q,qz,f,ev,ew,oprint,oconvd,ocz)
c
c build total scaling factor matrix and perform scaling
c 
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      complex *16 ev,sum,mat1,mat2,qz
c
c 20/6/07      parameter (m100 =400000)
      parameter (m100 =900000)
      common/somatrix1/mat1(m100),mat2(m100)
      common/somatrix2/rmat1(m100),rmat2(m100)
INCLUDE(common/occnr_z)
      common/ezora/eunscaled,escaled
INCLUDE(common/sizes)
INCLUDE(common/infoa) 
INCLUDE(common/zorac)
INCLUDE(common/iofile)
c     complex *16 ev,sum,mat1,mat2,qz
c
      dimension q(*),f(*),ev(2*num,*),qz(*),ew(*)
c
      l2 = num*(num+1)/2
      l3 = num*num
      isd1 = igmem_alloc(12*l3 +2*num)
      isd2 = isd1 + 4*l3
      isd3 = isd2 + 4*l3
      iscal = isd3 + 4*l3
c
      if (opzora) call prtri(f,num) 
c
c     call mult2(s2,q(isd1),f,num,num,num)
c     call dcopy(l2,q(isd1),1,f,1)
c
      if (opzora) then
         write(iwr,*) 'scaling factors'
         call prtri(f,num)
      endif
      call vclr(q(isd1),1,4*l3)
      n=-num
      mn=0
      do i=1,num
         n=n+num
         do j=1,num
            n=n+1
            mn=mn+1
            imax = max(i,j)
            ij = min(i,j) + imax*(imax-1)/2
c           qz(n) = (1.0d0,0.0d0)*f(ij) !+ (0.0d0,1.0d0)*rmat2(mn)
            qz(n) = (1.0d0,0.0d0)*f(ij) - (0.0d0,1.0d0)*rmat2(mn)
c           qz(n) = (0.0d0,1.0d0)*rmat2(mn)
c           q(isd1+n-1) = f(ij)
         enddo
      enddo
c
      mn=0
      n=2*num*num
      do i=1,num
         n=n+num
         do j=1,num
            mn=mn+1
            n=n+1
            imax = max(i,j)
            ij = min(i,j) + imax*(imax-1)/2
c           qz(n)=(1.0d0,0.0d0)*f(ij) !- (0.0d0,1.0d0)*rmat2(mn)
            qz(n)=(1.0d0,0.0d0)*f(ij) + (0.0d0,1.0d0)*rmat2(mn)
c           qz(n)=- (0.0d0,1.0d0)*rmat2(mn)
c           q(isd1+n-1) = f(ij)
         enddo
      enddo
c
      if (opzora)  then
c
      write(iwr,*)'scaling factor matrix (scalar part)'
      n = 2*num
c     do i=1,4*num*num
c        q(isd1+i-1) = dreal(qz(i))
c     enddo
      call dcopy (4*num*num,qz,2,q(isd1),1)
      call prsq(q(isd1),n,n,n)
      write(iwr,*) 'so corrections to scaling factor matrix'
      write(iwr,*) 'rmat2'
      call prsq(rmat2,num,num,num)
      write(iwr,*) 'mat2 real part'
      do i=1,num*num
         q(isd2+i-1) = dreal(mat2(i))
      enddo
      call prsq(q(isd2),num,num,num)
      do i=1,num*num
         q(isd2+i-1) = dimag(mat2(i))
      enddo
      write(iwr,*) 'imaginary part'
      call prsq(q(isd2),num,num,num)
c
      end if
c
c
c ok ? without mat2
c
c
      oflop = .true.
      if (oflop) then
      n=0
      mn = 0
      do i=1,num
         n=n+num
         do j=1,num
            n=n+1
            mn=mn+1
            qz(n) = qz(n) - dconjg(mat2(mn))
         enddo
      enddo
c
c
      n=2*num*num-num
      mn = 0
      do i=1,num
         n=n+num
         do j=1,num
            n=n+1
            mn=mn+1
             qz(n) = qz(n) + (mat2(mn))
         enddo
      enddo
      endif
c
      if (opzora) then
          write(iwr,*) 'total matrix of scaling factors'
          call prevz(q,qz,num)
      endif
c
c
c ... perform matrix multiplication in order to get scaling factors
c ... in MO basis
c
c ... real part first
c
c
      n=2*num
c     print *,'orbital and scalingfactor'
      do k=1,n
         sum= (0.0d0,0.0d0)
         do i=1,n
            do j=1,n
               mn = (i-1)*n + j
c              ij = min(i,j) + iky(max(i,j))
c
c ok in combination with only rmat2
c
c 
c              if (.not. oflop) then
c              sum = sum - dconjg(ev(j,k))*ev(i,k)*dreal(qz(mn)) 
c   
c              else
               sum = sum - dconjg(ev(j,k))*ev(i,k)*(qz(mn)) 
c              endif
            
c              sum = sum - dconjg(ev(j,k))*ev(i,k)*q(isd1+mn-1)
            enddo
         enddo
c        print 66,k,sum 
         q(iscal+k-1) = dreal(sum)
      enddo 
      if (oconvd.and.ocz)
     + write(iwr,*) 'unscaled and scaled eigenvalues'
      eunscaled = 0.0d0
      escaled = 0.d0
      if (oconvd.and.ocz) then
          do i=1,ne+2
             write(iwr,67) ew(i),ew(i)/(1.0d0+q(iscal+i-1))
          enddo
      end if 
      do i=1,2*num
c        print 67,ew(i),ew(i)/(1.0d0+scal(i))
         eunscaled = eunscaled + xocc_z(i)*ew(i)
         escaled = escaled + xocc_z(i)*ew(i)/(1.0d0+q(iscal+i-1))
      enddo
      if (opzora) print 666,eunscaled,escaled
c
      call gmem_free(isd1)
      return
66    format (i3,'factor',2f15.8)
67    format (f14.8,5x,f14.8)
666   format ('sum of unscaled and scaled orb. eng.',2f14.7)
      end 
c
c***************************************************************
c
      subroutine block_z(s,ndim,mode)
c
      implicit REAL  (a-h,p-z),integer   (i-n),logical    (o)
c
c...   subroutine to block 1-integral matrix over atoms 
c...   for atomic startup
c
c...   if mode = rect: s rectagular num2*num 
c                      num2 internal basis (big) = ndim
c...                   num exterbal basis (small)
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
      common/junk/ilifc3(maxorb),ntran3(maxorb),itran3(mxorb3),
     * ctran3(mxorb3),itrn3(2),
     *ex2(mxprim),cs2(mxprim),cp2(mxprim),cd2(mxprim),cf2(mxprim),
     *cg2(mxprim),cspace(maxat),
     +kstrt2(mxshel),ktom2(mxshel),ktype2(mxshel),kng2(mxshel),
     +kloc2(mxshel),kmin2(mxshel),kmax2(mxshel),nshel2,nspace(3),
     *nprin1,itol1,icut1,normf1,normp1
     *,nprini(695),iovmat,iosvec,ioproj,iorthg,c2(3,maxat)
c
      character*(*) mode
      dimension s(ndim,*)
c
      itrian(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j) 
c
c...  loop over shells 
c
      if (mode .eq. 'tri') then
c
        do ii = 1 , nshell
            do jj = 1 , ii-1 
               if (katom(ii).ne.katom(jj)) then
                  nni = kmax(ii)-kmin(ii)
                  nnj = kmax(jj)-kmin(jj)
                  do iorb = kloc(ii) , kloc(ii)+nni
                     do jorb =  kloc(jj) , kloc(jj)+nnj
                        s(itrian(iorb,jorb),1) = 0.0d0
                     end do
                  end do
               end if
            end do
        end do
c
      else if (mode .eq. 'rect') then 
c
        call caserr(' obsolete option')
c
        do ii = 1 , nshell
            do jj = 1 , nshel2
               if (katom(ii).ne.ktom2(jj)) then
                  nni = kmax(ii)-kmin(ii)
                  nnj = kmax2(jj)-kmin2(jj)
                  do iorb = kloc(ii) , kloc(ii)+nni
                     do jorb =  kloc2(jj) , kloc2(jj)+nnj
                        s(iorb,jorb) = 0.0d0
                     end do
                  end do
               end if
            end do
        end do
      else
         call caserr ('invalid mode in block_z')
      end if 
c
      return
      end
c
c**********************************************************************
      subroutine dens_eweight(A,B,C,l,ne) !renamed from dichtheid_gewogen -XJ
c...  this subroutine computes the energy weighed density matrix
c...  A is the matrix of eigenvectors, B is the eigenvalues, C is the output
      implicit   none
      integer    l,ne,i,j,k,m,n
      COMPLEX    A(l,l),C(l,l),zero,temp
      parameter   (zero=(0.0d+0, 0.0d+0))
      REAL       B(l)
      do k=1,l
       do i=1,l
        C(k,i) = zero
        do j=1,ne
         temp = -dcmplx(B(j),0.0d0)
         C(k,i) = C(k,i)+temp*conjg(A(k,j))*A(i,j)
        enddo
       enddo
      enddo
      return
      end
c**********************************************************************
      subroutine densdecomplex(A,B,l) !renamed from dichtheid_test1 -XJ
c...  this subroutine computes the density as suggested by Simon Faas
      implicit   none
      COMPLEX    A(2*l,2*l)
      REAL       B(l,l)
      integer     l,i,k
      do k=1,l
       do i=1,l
        B(k,i) = real(A(k,i))+real(A(k+l,i+l))
       enddo
      enddo     
      return
      end
c*********************************************************************
      subroutine ortho_check_kolom(A,B,l)
c...  this subroutine computes B=A*conjg(A'). A and B are of size l*l
c...  the resulting matrix should be the unit matrix if the columns are
c...  orthogonal
      implicit    none
      COMPLEX     A(l,l),B(l,l),zero
      parameter   (zero=(0.0d+0, 0.0d+0))
      integer     l,i,j,k
      real        restreal,matrimag,diagreal
INCLUDE(common/iofile)
      restreal = 0.
      matrimag = 0.
      diagreal = 0.
      do k=1,l
       do i=1,l
        B(k,i) = zero
        do j=1,l
         B(k,i) = B(k,i)+A(k,j)*conjg(A(i,j))
        enddo
       enddo
      enddo
c...  matrix B now should be a complex unit matrix size l*l      
      do k=1,l
       do i=k,l
        if(i.eq.k)diagreal = diagreal + abs(real(B(k,i)))-1.
        if(i.ne.k)restreal = restreal + abs(real(B(k,i)))
        matrimag = matrimag + abs(imag(B(k,i)))
       enddo
      enddo
      write(iwr,*)'test for orthogonality, residues:'
      write(iwr,*)'restreal,diagreal,matrimag',
     +restreal,diagreal,matrimag
      return 
      end
c*********************************************************************
      subroutine nonrel_atom(addzora)
c
c...  allow for non-relativistic atoms (only in zora atomic)
c...  by clearing the corresponding elements of the zora corrections
c...  Thus the zora corrections are calculated but ignored
c
      implicit REAL  (a-h,p-z),integer   (i-n),logical    (o)
c
      dimension addzora(*)
c
INCLUDE(common/sizes)
INCLUDE(common/iofile)
INCLUDE(common/infoa)
INCLUDE(common/nshel)
INCLUDE(common/runlab)
INCLUDE(common/datgue)
c
      save ofirst
      data ofirst/.true./
c
      itrian(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j) 
c
      if (natconf.eq.0) return
c
c...  loop over nonrelativistic atoms
c
      do iat=1,natconf
         if (nonrel(iat)) then
            if (ofirst) write(iwr,600) zatconf(iat)
600         format(4x,'**NONRELATIVISTIC** ATOM ',a8)
            ofirst = .false.
c
c...        loop over all atoms
c
            do iatom=1,nat
c
               if (zaname(iatom).eq.zatconf(iat)) then
c
c...               loop over shells 
c
                  do ii = 1 , nshell
                     if (katom(ii).eq.iatom) then
                        nni = kmax(ii)-kmin(ii)
                        do jj= 1 , nshell
                           if (katom(jj).eq.iatom) then
                              nnj = kmax(jj)-kmin(jj)
                              do iorb = kloc(ii) , kloc(ii)+nni
                                 do jorb = kloc(jj) , kloc(jj)+nnj
                                    addzora(itrian(iorb,jorb)) = 0.0d0
                                 end do
                              end do
                           end if
                        end do
                     end if
                  end do
               end if
            end do
         end if
      end do
c
      ofirst = .false.
c
      return
      end
c *********************************************************************

