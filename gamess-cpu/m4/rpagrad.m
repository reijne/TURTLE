_IF(rpagrad)
      block data initrpagrad
      implicit none
c     Huub van Dam, 1999
INCLUDE(common/sizes)
INCLUDE(common/rpadcom)
      data orpagrad/.false./
c     data orpaenrgy/.false./
c     data orpaegrad/.false./
      end
c
c-----------------------------------------------------------------------
c
      subroutine rpdinit()
      implicit none
c
c...  Reads the input specific to RPA gradients
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/rpacom)
INCLUDE(common/rpadcom)
INCLUDE(common/work)
INCLUDE(common/iofile)
INCLUDE(common/gjs)
c
c...  Local variables
c
      integer i,j, ntotstat
      character*4 ytest
c
      orpagrad  = .true.
c     orpaenrgy = .false.
c     orpaegrad = .false.
c
      ntotstat = 0
      do i = 1, 8
         rpad_nstates(i) = 0
         do j = 1, max_rpad_states
            rpad_istates(j,i) = 0
         enddo
      enddo
c
      jrec = 0
      call inpa4(ytest)
      if (ytest.ne.'rpag') then
         jrec = jrec - 1
      else
 10      call input
         call inpa4(ytest)
         if (ytest.eq.'end') then
         else if (ytest.eq.'symm') then
            call inpi(i)
            if (i.lt.1.or.i.gt.8) then
               write(iwr,630)1,nirr,i
               call caserr('Illegal symmetry specified')
            endif
 20         if (jrec.lt.jump) then
               ntotstat = ntotstat + 1
               if (ntotstat.gt.max_rpad_states) then
                  call caserr('# gradients requested exceeds parameter m
     +ax_rpad_states')
               endif
               rpad_nstates(i) = rpad_nstates(i)+1
               call inpi(rpad_istates(rpad_nstates(i),i))
               goto 20
            endif
         else
            call caserr('invalid RPA gradients directive')
         endif
         if (ytest.ne.'end') go to 10
      endif
      nrpastate = ntotstat
c
c...  Check input data
c
      do i = 1, 8
         call consort(rpad_istates(1,i),rpad_nstates(i),0)
         do j = 1, rpad_nstates(i)-1
            if (rpad_istates(j,i).eq.rpad_istates(j+1,i)) then
               write(iwr,600)rpad_istates(j,i),i
               call caserr('State specified more than once in RPA gradie
     +nts directive')
            endif
         enddo
         if (rpad_nstates(i).ne.0) then
            if (rpad_istates(1,i).lt.nevlo(i).or.
     +          rpad_istates(rpad_nstates(i),i).gt.nevhi(i)) then
               write(iwr,610)i,nevlo(i),nevhi(i)
               write(iwr,620)
               call caserr('States for RPA gradients must be a subset of
     + the RPA states')
            endif
         endif
      enddo
 600  format('*** ERROR: State ',i2,' of symmetry ',i1,
     +       ' specified more than once')
 610  format('*** ERROR: Some RPA gradient states for symmetry ',i1,
     +       ' are not among the RPA states ',i2,'-',i2)
 620  format('*** ERROR: Gradients are implemented only for a subset ',
     +       'of the specified RPA states.')
 630  format('*** ERROR: Only irreps ',i1,' to ',i1,' are possible. ',
     +       'Irrep ',i3,' is out of bounds.')
      end
c
c-----------------------------------------------------------------------
c
      subroutine drpadr(core)
c
c     Driver for the RPA+Gradients calculations
c
c     1. compute 1&2 electron integrals
c     2. carryout scf (rpa,dirrpa) 
c     3. perform integral transformation
c     4. perform response calculation
c     5. evaluate (50-55,58) in Ortiz (i.e. Gamma, I, J and T)
c     6. compute stuff for CPHF equations including equation (63)
c     7. solve CPHF equations
c     8. evaluate (57) and add (57) and (63).
c
c     Huub van Dam, 1999
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
INCLUDE(common/maxlen)
INCLUDE(common/cslosc)
INCLUDE(common/atmol3)
INCLUDE(common/restar)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/restrj)
INCLUDE(common/symtry)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/timez)
INCLUDE(common/runlab)
INCLUDE(common/funct)
INCLUDE(common/rpadcom)
INCLUDE(common/vectrn)
INCLUDE(common/nshel)
      common/scfblk/enuc,etotal,ehfock,sh1(17),osh1(7),ish1(7)
      common/blkin/potnuc(10)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      dimension core(*),zcas(2)
c
      character*9 fnm
      character*6 snm
      data fnm,snm/"rpagrad.m","drpadr"/
c
      data oyes,ono / .true.,.false./
      data zcas  / 'casscf','mcscf'/
      data m10,m16/10,16/
c
      dum = 0.0d0
      idum = 0
      ntsave = nt
      omccas = zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2)
      iofsym = 0
      nopk = 1
      if (.not.odirpa) then
         nt = 1
         iofsym = 1
      endif
      enrgy = 0.0d0
      isecvv = isect(8)
      itypvv = 8
      if (opass2) then
         if(nopk.eq.0.or.nopk.ne.nopkr.or.iofsym.eq.1.or.
     +      iofrst.ne.iofsym) then
            write (iwr,130)
            opass2 = .false.
         endif
      endif
 130  format(/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/
     + 1x,'= integrals must NOT be in supermatrix form ='/
     + 1x,'=        requested bypass is ignored        ='/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/)
      if (.not.opass1) then
         if (irest.lt.1) then
            call standv(0,core)
            call putstv(core)
            if (tim.ge.timlim) go to 80
         end if
      end if
      ionsv = 1
c
      if (opass2) then
         if(.not.ognore.and..not.odirpa) 
     &        call checkinteg(nopk,nopkr,iofsym,iofrst)
      else
         if(.not.odscf) then
            if (irest.le.1) then
               iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                                  "prefac",IGMEM_DEBUG)
               call rdmake(core(iprefa))
               call jandk(zscftp,core,core(inull),core(inull),
     +                    core(inull),core(inull),core(inull),
     +                    core(iprefa),core(inull))
               call gmem_free_inf(iprefa,fnm,snm,"prefac")
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6010)
                  go to 80
               end if
            end if
         end if
      end if
c
      if (.not.opass4) then
         if (irest.le.2) then
            call adapti(core,ono)
            if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6020)
               go to 80
            end if
         end if
      end if
      nt = ntsave
c
      if (.not.opass5) then
         if (irest.le.4) then
            call scf(core)
            if (odscf. and. 
     +      (guess.eq.'atoms'))go to 1000
c           if(.not.omccas) then
               call putdev(core,mouta,7,1)
c
               call secget(isect(494),m16,iblk34)
               call rdedx(potnuc,m10,iblk34,idaf)
               enuc = potnuc(1)
               ehfock = potnuc(2)
               etotal = potnuc(3)
c
               call secget(isect(13),13,iblok)
               call wrt3(enuc,lds(isect(13)),iblok,idaf)
c           endif
1000        if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6030)
               go to 80
            end if
         end if
      end if
c
      if (.not.opass5) then
         if (.not.opass6) then
            if (irest.le.8) then
               call adapti(core,oyes)
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6040)
                  go to 80
               end if
            end if
         end if
      end if
c
c     linear response
c
      call rpdriv(core,core)
c
      call pre_rpa_grad(core)
c
      fkder  = 'fockder'
      call hfgrdn(core)
      ibase  = igmem_alloc_all(mword)
      maxq   = mword
      lword4 = mword
      call indxsy(core(ibase))
      call indx2t(core(ibase))
      call indx4t(core(ibase))
      call trnfkd(core(ibase))
c     call chfndr(core(ibase),core(ibase))
      call gmem_free(ibase)
      call chfndr(core,core)
c
      call rpa_grad(core)
      call pr_rpa_grad()
c
c     update restart data
c
      opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
      opass6 = ono
 80   itask(mtask) = irest
      nt = ntsave
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
c
      return
 6010 format (/1x,'termination due to incomplete integral',
     +        ' evaluation')
 6020 format (/1x,'termination due to incomplete integral',
     +        ' adaptation')
 6030 format (/1x,'termination due to incomplete scf')
 6040 format (/1x,'termination due to incomplete integral',
     +        ' transformation')
      end
c
c-----------------------------------------------------------------------
c
      subroutine drpadr_e(core)
c
c     Driver for the RPA geometry optimisations: energy calculation.
c
c     This routine is based on drpadr. It is assumed to be used in 
c     geometry optimisations. In these calculations it computes only the 
c     energy of the state of interest. The gradient part is taken care
c     of by drpadr_g.
c
c     1. compute 1&2 electron integrals
c     2. carryout scf (rpa,dirrpa) 
c     3. perform integral transformation
c     4. perform response calculation
c     5. evaluate the total energy
c
c     Huub van Dam, 2000
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
INCLUDE(common/gmempara)
c
INCLUDE(common/maxlen)
INCLUDE(common/cslosc)
INCLUDE(common/atmol3)
INCLUDE(common/restar)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/restrj)
INCLUDE(common/symtry)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/timez)
INCLUDE(common/runlab)
INCLUDE(common/funct)
INCLUDE(common/rpadcom)
INCLUDE(common/vectrn)
INCLUDE(common/nshel)
      common/scfblk/enuc,etotal,ehfock,sh1(17),osh1(7),ish1(7)
      common/blkin/potnuc(10)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      dimension core(*),zcas(2)
c
      character*9 fnm
      character*8 snm
      data fnm,snm/"rpagrad.m","drpadr_e"/
c
      data oyes,ono / .true.,.false./
      data zcas  / 'casscf','mcscf'/
      data m10,m16/10,16/
c
      if (nrpastate.ne.1) then
         call caserr('RPA-optimise: only 1 state at a time')
      endif
      dum = 0.0d0
      idum = 0
      ntsave = nt
      omccas = zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2)
      iofsym = 0
      nopk = 1
      if (.not.odirpa) then
         nt = 1
         iofsym = 1
      endif
      enrgy = 0.0d0
      isecvv = isect(8)
      itypvv = 8
      if (opass2) then
         if(nopk.eq.0.or.nopk.ne.nopkr.or.iofsym.eq.1.or.
     +      iofrst.ne.iofsym) then
            write (iwr,130)
            opass2 = .false.
         endif
      endif
 130  format(/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/
     + 1x,'= integrals must NOT be in supermatrix form ='/
     + 1x,'=        requested bypass is ignored        ='/
     + 1x,'= = = = = = = = = = = = = = = = = = = = = = ='/)
c
      if (.not.opass1) then
         if (irest.lt.1) then
            call standv(0,core)
            call putstv(core)
            if (tim.ge.timlim) go to 80
         end if
      end if
      ionsv = 1
c
      if (opass2) then
         if(.not.ognore.and..not.odirpa) 
     &        call checkinteg(nopk,nopkr,iofsym,iofrst)
      else
         if(.not.odscf) then
            if (irest.le.1) then
               iprefa = igmem_alloc_inf(nshell*(nshell+1)/2,fnm,snm,
     +                                  "prefac",IGMEM_DEBUG)
               call rdmake(core(iprefa))
               call jandk(zscftp,core,core(inull),core(inull),
     +                    core(inull),core(inull),core(inull),
     +                    core(iprefa),core(inull))
               call gmem_free_inf(iprefa,fnm,snm,"prefac")
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6010)
                  go to 80
               end if
            end if
         end if
      end if
c
      if (.not.opass4) then
         if (irest.le.2) then
            call adapti(core,ono)
            if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6020)
               go to 80
            end if
         end if
      end if
      nt = ntsave
c
      if (.not.opass5) then
         if (irest.le.4) then
            call scf(core)
            if (odscf. and. 
     +      (guess.eq.'atoms'))go to 1000
            call putdev(core,mouta,7,1)
c
            call secget(isect(494),m16,iblk34)
            call rdedx(potnuc,m10,iblk34,idaf)
            enuc = potnuc(1)
            ehfock = potnuc(2)
            etotal = potnuc(3)
c
            call secget(isect(13),13,iblok)
            call wrt3(enuc,lds(isect(13)),iblok,idaf)
1000        if (tim.ge.timlim .or. irest.ne.0) then
               write (iwr,6030)
               go to 80
            end if
         end if
      end if
c
      if (.not.opass5) then
         if (.not.opass6) then
            if (irest.le.8) then
               call adapti(core,oyes)
               if (tim.ge.timlim .or. irest.ne.0) then
                  write (iwr,6040)
                  go to 80
               end if
            end if
         end if
      end if
c
c     linear response
c
      call rpdriv(core,core)
c
c     store the total energy for geometry optimisation
c
      call drpa_enrgy()
c
c     update restart information
c
      opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
      opass6 = ono
 80   itask(mtask) = irest
      nt = ntsave
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
c
      return
 6010 format (/1x,'termination due to incomplete integral',
     +        ' evaluation')
 6020 format (/1x,'termination due to incomplete integral',
     +        ' adaptation')
 6030 format (/1x,'termination due to incomplete scf')
 6040 format (/1x,'termination due to incomplete integral',
     +        ' transformation')
      end
c
c-----------------------------------------------------------------------
c
      subroutine drpadr_g(core)
c
c     Driver for the RPA geometry optimisations: Gradient part
c
c     This subroutine computes the gradient of an excited state within
c     the RPA formalism. It is assumed that the energy has already been
c     calculated by drpadr_e, and that all required data has been
c     stored appropriately. This routine executes only the additional
c     steps needed to obtain the gradient.
c
c     5. evaluate (50-55,58) in Ortiz (i.e. Gamma, I, J and T)
c     6. compute stuff for CPHF equations including equation (63)
c     7. solve CPHF equations
c     8. evaluate (57) and add (57) and (63).
c     9. compute the total gradient and store in egrad.
c
c     Huub van Dam, 2000
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/sizes)
c
INCLUDE(common/maxlen)
INCLUDE(common/cslosc)
INCLUDE(common/atmol3)
INCLUDE(common/restar)
      common/restri/jjfile(63),lds(508),isect(508),ldsect(508),
     +              iacsct(508)
INCLUDE(common/cndx40)
INCLUDE(common/cndx41)
INCLUDE(common/restrj)
INCLUDE(common/symtry)
INCLUDE(common/infoa)
INCLUDE(common/iofile)
INCLUDE(common/timez)
INCLUDE(common/runlab)
INCLUDE(common/funct)
INCLUDE(common/rpadcom)
INCLUDE(common/vectrn)
      common/scfblk/enuc,etotal,ehfock,sh1(17),osh1(7),ish1(7)
      common/blkin/potnuc(10)
c
      character *8 title,scftyp,runtyp,guess,conf,fkder
      character *8 scder,dpder,plder,guesc,rstop,charsp
      common/restrz/title(10),scftyp,runtyp,guess,conf,fkder,
     + scder,dpder,plder,guesc,rstop,charsp(30)
c
      dimension core(*),zcas(2)
c
      data oyes,ono / .true.,.false./
      data zcas  / 'casscf','mcscf'/
      data m10,m16/10,16/
c
      if (nrpastate.ne.1) then
         call caserr('RPA-optimise: only 1 state at a time')
      endif
      dum = 0.0d0
      idum = 0
      ntsave = nt
      omccas = zscftp.eq.zcas(1) .or. zscftp.eq.zcas(2)
c
      call pre_rpa_grad(core)
c
      fkder  = 'fockder'
      call hfgrdn(core)
      ibase  = igmem_alloc_all(mword)
      maxq   = mword
      lword4 = mword
      call indxsy(core(ibase))
      call indx2t(core(ibase))
      call indx4t(core(ibase))
      call trnfkd(core(ibase))
c     call chfndr(core(ibase),core(ibase))
      call gmem_free(ibase)
      call chfndr(core,core)
c
      call rpa_grad(core)
      call pr_rpa_grad()
c
c     we need to store the gradient in common/funct/egrad
c     for the geometry optimisation.
c
      call drpa_egrad()
c
c     update restart information
c
      opass1 = ono
      opass2 = ono
      opass3 = ono
      opass4 = ono
      opass5 = ono
      opass6 = ono
 80   itask(mtask) = irest
      nt = ntsave
      call wrrec1(idum,idum,idum,idum,c,c,dum,dum,dum,c,c)
c
      return
 6010 format (/1x,'termination due to incomplete integral',
     +        ' evaluation')
 6020 format (/1x,'termination due to incomplete integral',
     +        ' adaptation')
 6030 format (/1x,'termination due to incomplete scf')
 6040 format (/1x,'termination due to incomplete integral',
     +        ' transformation')
      end
c
c-----------------------------------------------------------------------
c
      subroutine pre_rpa_grad(q)
      implicit none
      REAL q(*)
c
c...  This subroutine calculates the quantities going into the RPA 
c...  gradient that don't require the derivatives of the integrals
c...  nor the CPHF solutions.
c
c...  The RPA gradient is given by 
c...  (see J.V. Ortiz, J. Chem. Phys., Vol 101 (1994), No 8, 6743-6749)
c...  (MO indeces ijkl are occupied MOs and abcd are unoccupied MOs.)
c...  Beware of the error in Ortiz's paper, the diagonal elements of
c...  Sa should be included in the summations!
c
c...  Ea = EaMO + EaAO
c
c...  EaMO = Sum(ij,Gamma_ij Qa_ij)
c...       + Sum(ab,Gamma_ab Qa_ab)
c...       + Sum(ij,Sa_ij I_ij)
c...       + Sum(ab,Sa_ab I_ab)
c...       + Sum(ai,Sa_ai I_ai)
c...       + Sum(ai,U_ai T_ai)                    (see Ortiz page 6747)
c
c...  Gamma_ij = - Sum(a, Z_ai Z_aj + Y_ai Y_aj)
c
c...  Gamma_ab = Sum(i, Z_ai Z_bi + Y_ai Y_bi)
c
c...  I_ij = Sum(a,- Z_ai Z_aj ( Epole + e_j - e_a) 
c...               - Y_ai Y_aj (-Epole + e_j - e_a) )
c
c...  I_ab = Sum(a,- Z_ai Z_bi ( Epole + e_i - e_b) 
c...               - Y_ai Y_bi (-Epole + e_i - e_b) )
c
c...  I_ai = -2 Sum(jkb, (Z_bj Y_ak + Y_bj Z_ak)<kj||bi>
c...                    +(Z_bj Z_ak + Y_bj Y_ak)<ji||bk> )
c
c...  T_ai =   2 Sum(jbc, (Z_bj Y_ci + Y_bj Z_ci)<aj||bc>
c...                     +(Z_bj Z_ci + Y_bj Y_ci)<jc||ba> )
c...         - 2 Sum(jkb, (Z_bj Y_ak + Y_bj Z_ak)<kj||bi>
c...                     +(Z_bj Z_ak + Y_bj Y_ak)<ji||bk> )
c...         + 2 Sum(k<l, Gamma_kl <ki||la>)
c...         + 2 Sum(b<c, Gamma_bc <bi||ca>)
c
c...  EaAO = Sum(pqrs, <pq|rs>a g_pqrs)           (see Ortiz page 6748)
c
c...  g_pqrs = (Z_rp Z_qs + Y_rp Y_qs)
c...         - (Z_sp Z_qr + Y_sp Y_qr)
c...         + (Z_rq Y_sp + Y_rq Z_sp)
c...         - (Z_sq Y_rp + Y_sq Z_rp)
c
c...  Z_pq = sum(ia,C_pa C_qi Z_ai)
c
c...  Y_pq = sum(ia,C_pa C_qi Y_ai)
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/rpadcom)
INCLUDE(common/rpacom)
INCLUDE(common/machin)
INCLUDE(common/infoa)
INCLUDE(common/sector)
INCLUDE(common/atmol3)
INCLUDE(common/restri)
INCLUDE(common/mapper)
INCLUDE(common/iofile)
      integer icount,kcore,kvirt,ksymm
      integer lcore,lvirt,iofset,npair
      integer nbox,nbox0,jstart,iofsym
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
c
c...  File stuff
c
      integer irpatyp
      parameter(irpatyp =120)
      integer irpatyp1, irpatyp2
      parameter(irpatyp1=121)
      parameter(irpatyp2=122)
      integer irpasec, irpasec1, irpasec2
      parameter(irpasec =351)
      parameter(irpasec1=352)
      parameter(irpasec2=353)
      integer irpablk
      integer nwords(8), lensect
      REAL rwords
      equivalence(rwords,nwords(1))
c
c...  Symmetry stuff
c
      integer isym
c
c...  RPA results
c
      REAL epole(2), dnrm
      integer ipo, ipz, ipy
c
c...  Some counters
c
      integer i,istate,nrpaln2,irpast
c
c...  Some addresses
c
      integer ipgab, ipgij, igij, igab, iij, iab, iblk1, iai, itai
      integer izsu, iysu
      integer ic, icblk
      integer ieps, iepsblk
c
c...  Some sizes
c
      integer nrpalen, isizeiij, isizeiab
      integer ncore
c
c...  Some functions
c
      integer lensec, lenint, igmem_alloc
      REAL dnrsq
c
c...  Perform a symmetry check
c
      do isym = nirr+1, 8
         if (rpad_nstates(isym).ne.0) then
            write(iwr,600)isym,1,nirr
            call caserr('Gradients requested for non-present symmetry')
         endif
      enddo
 600  format('Gradients for irrep ',i1,' requested but only irreps ',
     +       i1,' to ',i1,' present.')
c
c...  Zero the excitation energy gradients
c
      call vclr(rpade,1,3*maxat*max_rpad_states)
c
c...  Load the MO coefficients
c
      ic   = igmem_alloc(num*num)
      call secget(mouta,3,icblk)
      icblk = icblk+1+lensec(mach(8))+lensec(mach(9))
      call rdedx(q(ic),num*num,icblk,numdu)
      call tdown(q(ic),ilifq,q(ic),ilifq,num)
c
c...  Load the orbital energies
c
      ieps = igmem_alloc(num)
      call secget(isect(9),9,iepsblk)
      call rdedx(q(ieps),num,iepsblk,numdu)
c
c...  Load the Z and Y component mapping
c
      ipo = igmem_alloc(lenint(nbasq))
      call secget(irpasec,irpatyp,irpablk)
      call rdedx(q(ipo),lenint(nbasq),irpablk,numdu)
      nrpalen = lensec(lenint(nbasq))
c
c...  Allocate memory for the Gamma_ab and Gamma_ij component mappings
c
      ncore = na
      ipgab = igmem_alloc(lenint(nvirt*nvirt))
      ipgij = igmem_alloc(lenint(ncore*ncore))
      call setsto(nvirt*nvirt,0,q(ipgab))
      call setsto(ncore*ncore,0,q(ipgij))
c
c...  Allocate memory for Gamma_ab and Gamma_ij itself
c
      isizegij = 0
      isizegab = 0
      do i = 1, nirr
         isizegij = isizegij+(kcore(i)+1)*kcore(i)/2
         isizegab = isizegab+(kvirt(i)+1)*kvirt(i)/2
      enddo
      igij = igmem_alloc(isizegij)
      igab = igmem_alloc(isizegab)
c
c...  Allocate memory for the I_ab and I_ij terms
c
      isizeiij = ncore*ncore
      isizeiab = nvirt*nvirt
      iij = igmem_alloc(isizeiij)
      iab = igmem_alloc(isizeiab)
c
c...  Allocate memory for the I_ab and I_ij terms
c
      iai  = igmem_alloc(nvirt*ncore)
      itai = igmem_alloc(nvirt*ncore)
c
c...  Allocate memory for Z_su and Y_su (i.e. Z and Y in AO basis)
c
      izsu = igmem_alloc(num*num)
      iysu = igmem_alloc(num*num)
c
c...  Create a new section for storing all intermediate results.
c...  The length of this section is not know yet.
c
      lensect = 0
      call secput(irpasec1,irpatyp1,lensect,iblk1)
c
c...  Per symmetry we go through all required excited states and 
c...  and calculate the things we need for the gradients
c
      nrpastate = 0
      do isym = 1, nirr
c
c...     Allocate memory to store the RPA eigenvectors
c
         ipz = igmem_alloc(icount(isym))
         ipy = igmem_alloc(icount(isym))
c
c...     Loop over all require RPA states of this symmetry
c
         nrpastate = nrpastate+rpad_nstates(isym)
chvd     do istate = 1, rpad_nstates(isym)
         irpast = 1
         do 10 istate = 0, nev(isym)-1
            if (rpad_istates(irpast,isym).ne.istate+nevlo(isym)) 
     +          go to 10
            irpast = irpast + 1
c
c...        Load the excitation energy, Z-vector and Y-vector
c
chvd        nrpaln2 = (rpad_istates(istate,isym)-1)*
chvd +                (2*lensec(icount(isym))+lensec(2))
            nrpaln2 = (istate)*
     +                (2*lensec(icount(isym))+lensec(2))
            call rdedx(epole,2,irpablk+nrpalen+nrpaln2,numdu)
            nrpaln2 = nrpaln2+lensec(2)
            call rdedx(q(ipz),icount(isym),irpablk+nrpalen+nrpaln2,
     +                 numdu)
            nrpaln2 = nrpaln2+lensec(icount(isym))
            call rdedx(q(ipy),icount(isym),irpablk+nrpalen+nrpaln2,
     +                 numdu)
            dnrm = dnrsq(icount(isym),q(ipy),1)
     +           - dnrsq(icount(isym),q(ipz),1)
            dnrm = dabs(dnrm)
            do i = 0, icount(isym)-1
               q(ipy+i) = q(ipy+i)/dsqrt(dnrm)
               q(ipz+i) = q(ipz+i)/dsqrt(dnrm)
            enddo
c
c...        Compute the RPA energy (just to check if we are doing the 
c...        right things)
c
            call rpa_energy(q(ipo),isym,q(ipz),q(ipy),q(ieps))
c
c...        Compute Gamma_ab and Gamma_ij.
c...        Note that Gamma_ab = Gamma_ba and that Gamma_ij = Gamma_ji
c...        so we compute only triangular matrices.
c
            call gamma_ij_ab(q(ipo),num,isym,q(ipz),q(ipy),
     +                       q(ipgij),q(igij),q(ipgab),q(igab))
c
c...        Compute I_ab and I_ij.
c
            call i_ij_ab(q(ipo),num,isym,q(ieps),epole,q(ipz),q(ipy),
     +                   q(iij),q(iab))
c
c...        Compute I_ai and T_ai
c
            call i_t_ai(q(ipo),num,isym,q(ipz),q(ipy),
     +                  q(ipgij),q(igij),q(ipgab),q(igab),
     +                  q(iai),q(itai))
c
c...        Compute Z_su and Y_su
c
            call z_y_su(q(ipo),num,isym,q(ipz),q(ipy),q(ic),
     +                  q(izsu),q(iysu))
c
c...        Store the results:
c...        Gamma_ij, Gamma_ab, I_ij, I_ab, I_ai, T_ai, Z_su, Y_su
c
            nwords(1) = isizegij
            nwords(2) = isizegab
            nwords(3) = isizeiij
            nwords(4) = isizeiab
            nwords(5) = nvirt*ncore
            nwords(6) = nvirt*ncore
            nwords(7) = num*num
            nwords(8) = num*num
            call wrt3(nwords,lenint(8),iblk1+lensect,numdu)
            lensect = lensect+1
            call wrt3(q(igij),nwords(1),iblk1+lensect,numdu)
            lensect = lensect+lensec(nwords(1))
            call wrt3(q(igab),nwords(2),iblk1+lensect,numdu)
            lensect = lensect+lensec(nwords(2))
            call wrt3(q(iij),nwords(3),iblk1+lensect,numdu)
            lensect = lensect+lensec(nwords(3))
            call wrt3(q(iab),nwords(4),iblk1+lensect,numdu)
            lensect = lensect+lensec(nwords(4))
            call wrt3(q(iai),nwords(5),iblk1+lensect,numdu)
            lensect = lensect+lensec(nwords(5))
            call wrt3(q(itai),nwords(6),iblk1+lensect,numdu)
            lensect = lensect+lensec(nwords(6))
            call wrt3(q(izsu),nwords(7),iblk1+lensect,numdu)
            lensect = lensect+lensec(nwords(7))
            call wrt3(q(iysu),nwords(8),iblk1+lensect,numdu)
            lensect = lensect+lensec(nwords(8))
 10      continue
         nrpalen = nrpalen+nev(isym)*(2*lensec(icount(isym))+lensec(2))
c
c...     Free memory for RPA vectors storage
c
         call gmem_free(ipy)
         call gmem_free(ipz)
      enddo
c
c...  Update length of irpasec1
c
      call secput(irpasec1,irpatyp1,lensect,iblk1)
c
c...  Now store the Gamma_ab and Gamma_ij component mappings
c
      lensect = lensec(lenint(ncore*ncore)) 
      lensect = lensec(lenint(nvirt*nvirt)) + lensect
      call secput(irpasec2,irpatyp2,lensect,iblk1)
      call wrt3(q(ipgij),lenint(ncore*ncore),iblk1,numdu)
      lensect = lensec(lenint(ncore*ncore))
      call wrt3(q(ipgab),lenint(nvirt*nvirt),iblk1+lensect,numdu)
c
c...  Free all memory 
c
      call gmem_free(iysu)
      call gmem_free(izsu)
      call gmem_free(itai)
      call gmem_free(iai)
      call gmem_free(iab)
      call gmem_free(iij)
      call gmem_free(igab)
      call gmem_free(igij)
      call gmem_free(ipgij)
      call gmem_free(ipgab)
      call gmem_free(ipo)
      call gmem_free(ieps)
      call gmem_free(ic)
      end
c
c-----------------------------------------------------------------------
c
      subroutine gamma_ij_ab(ipo,nbasis,irrep,z,y,ipgij,gij,ipgab,gab)
      implicit none
c
c...  Compute the core-core and virtual-virtual parts of Gamma
c...  (See subroutine pre_rpa_grad for the equations)
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      integer isym,isize,neig,iperm,isspac,maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
INCLUDE(common/infoa)
      integer nbasis,irrep,ipo(nbasis,nbasis)
      REAL z(*), y(*), gij(*), gab(*)
      integer ipgij(na,na), ipgab(nvirt,nvirt)
c
c...  Local variables
c
      integer ii, jj, k, ia, ib
c
      k = 0
      do jj = 1, na
         do ii = 1, jj
            if (isymmo(ii).eq.isymmo(jj)) then
               k = k + 1
               ipgij(jj,ii) = k
               ipgij(ii,jj) = k
               gij(k) = 0.0d0
               do ia = na+1, num
                  if (iperm(isymmo(ii),isymmo(ia)).eq.irrep.and.
     +                iperm(isymmo(jj),isymmo(ia)).eq.irrep) then
                     gij(k) = gij(k)
     +                      - z(ipo(ii,ia))*z(ipo(jj,ia))
     +                      - y(ipo(ii,ia))*y(ipo(jj,ia))
                  endif
               enddo
            endif
         enddo
      enddo
c
      k = 0
      do ib = na+1, num
         do ia = na+1, ib
            if (isymmo(ia).eq.isymmo(ib)) then
               k = k + 1
               ipgab(ia-na,ib-na) = k
               ipgab(ib-na,ia-na) = k
               gab(k) = 0.0d0
               do ii = 1, na
                  if (iperm(isymmo(ii),isymmo(ia)).eq.irrep.and.
     +                iperm(isymmo(ii),isymmo(ib)).eq.irrep) then
                     gab(k) = gab(k)
     +                      + z(ipo(ii,ia))*z(ipo(ii,ib))
     +                      + y(ipo(ii,ia))*y(ipo(ii,ib))
                  endif
               enddo
            endif
         enddo
      enddo
      end
c
c-----------------------------------------------------------------------
c
      subroutine i_ij_ab(ipo,nbasis,irrep,e,epole,z,y,iij,iab)
      implicit none
c
c...  Compute the core-core and virtual-virtual parts of I
c...  (See subroutine pre_rpa_grad for the equations)
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      integer isym,isize,neig,iperm,isspac,maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
INCLUDE(common/infoa)
      integer nbasis,irrep,ipo(nbasis,nbasis)
      REAL e(num), epole, z(*), y(*), iij(na,na), iab(nvirt,nvirt)
c
c...  Local variables
c
      integer jj, ii, k, ia, ib
c
      do jj = 1, na
         do ii = 1, na
            iij(ii,jj) = 0.0d0
            if (isymmo(ii).eq.isymmo(jj)) then
               do ia = na+1, num
                  if (iperm(isymmo(ii),isymmo(ia)).eq.irrep) then
                     iij(ii,jj) = iij(ii,jj)
     +                      -  z(ipo(ii,ia))*z(ipo(jj,ia))
     +                        *(epole+e(ii)-e(ia))
     +                      -  y(ipo(ii,ia))*y(ipo(jj,ia))
     +                        *(-epole+e(ii)-e(ia))
                  endif
               enddo
            endif
         enddo
      enddo
c
      do ib = na+1, num
         do ia = na+1, num
            iab(ia-na,ib-na) = 0.0d0
            if (isymmo(ia).eq.isymmo(ib)) then
               do ii = 1, na
                  if (iperm(isymmo(ii),isymmo(ia)).eq.irrep) then
                     iab(ia-na,ib-na) = iab(ia-na,ib-na) 
     +                      -  z(ipo(ii,ia))*z(ipo(ii,ib))
     +                        *(epole+e(ii)-e(ib))
     +                      -  y(ipo(ii,ia))*y(ipo(ii,ib))
     +                        *(-epole+e(ii)-e(ib))
                  endif
               enddo
            endif
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine i_t_ai(ipo,nbasis,irrep,z,y,ipgij,gij,ipgab,gab,
     +                  gai,tai)
      implicit none
c
c...  This routine computes Gamma_ai and T_ai. 
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/atmblk)
      REAL gin
      common/blkin /gin(511)
      integer i205
      common/junk2 /i205(4,340)
      integer nprint,n6file,n6tape,n6blk,n6last
      common/restar/nprint(142),n6file,n6tape(20),n6blk(20),n6last(20)
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      integer isym,isize,neig,iperm,isspac,maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),
     +              maxvec,mvec,mymaxv
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      integer mword
      equivalence(mword,gin(511))
c
c...  Parameters
c
      integer nbasis,ipo(nbasis,nbasis),irrep
      integer ipgij(na,na),ipgab(nvirt,nvirt)
      REAL gai(nvirt,na), tai(nvirt,na), z(*), y(*), gij(*), gab(*)
c
c...  Local variables
c
      integer ned,ifile,nw,iint,ii,jj,kk,ll,ie,mm,i
      integer l6blki,l6blk(20)
      integer isij,isik,isil,isjl,isim,isjm,iskj,iskl,iskm,islm
      integer isej,isel,isek,isjk,islj
      REAL gg
_IF(linux)
      external fget
_ENDIF
c
      call vclr(gai,1,nvirt*na)
      call vclr(tai,1,nvirt*na)
c
      do i = 1, n6file
         l6blk(i) = n6blk(i) - n6last(i)
      enddo
      call setsto(1360,0,i205)
      do ifile = 1, n6file
         ned = n6tape(ifile)
         call search(n6blk(ifile),ned)
         call find(ned)
         l6blki = l6blk(ifile)
10       l6blki = l6blki + 1
_IF(c)
         call get(gin,nw)
_ELSE
         call fget(gin,nw,ned)
_ENDIF
         if (nw.gt.0) then
            call unpack(gin(num2e+1),lab816,i205,numlab)
_IF(c)
            if (l6blki.ne.0) call find(ned)
_ENDIF
            do iint = 1, mword
_IF(littleendian)
               jj = i205(1,iint)
               ii = i205(2,iint)
               ll = i205(3,iint)
               kk = i205(4,iint)
_ELSE
               ii = i205(1,iint)
               jj = i205(2,iint)
               kk = i205(3,iint)
               ll = i205(4,iint)
_ENDIF
c
chvd           Note that the integrals are stored such that:
chvd           indxj <= indxi
chvd           indxl <= indxk
chvd           (indxk-1)indxk/2+indxl <= (indxi-1)indxi/2+indxj
chvd           which implies:
chvd           indxk <= indxi
c
               if (ii.le.na) then
c...              [ij|kl] (not used in I_ai or T_ai so skip)
               else
c...              [a*|**]
                  if (jj.le.na) then
c...                 [aj|**]
                     if (kk.le.na) then
c...                    [aj|kl]
                        gg = gin(iint)
                        if (kk.eq.ll) gg = 0.5d0*gg
                        do ie = 1, nvirt
                           mm   = na+ie
                           isel = iperm(isymmo(mm),isymmo(ll))
                           isij = iperm(isymmo(ii),isymmo(jj))
                           isej = iperm(isymmo(mm),isymmo(jj))
                           isil = iperm(isymmo(ii),isymmo(ll))
                           isek = iperm(isymmo(mm),isymmo(kk))
                           isik = iperm(isymmo(ii),isymmo(kk))
                           if (isel.eq.irrep.and.isij.eq.irrep) then
                              gai(ie,kk) = gai(ie,kk)+4.0d0*gg*(
     +                                   - z(ipo(mm,ll))*z(ipo(ii,jj))
     +                                   - y(ipo(mm,ll))*y(ipo(ii,jj))
     +                                   + z(ipo(mm,ll))*y(ipo(ii,jj))
     +                                   + y(ipo(mm,ll))*z(ipo(ii,jj)))
                           endif
                           if (isek.eq.irrep.and.isij.eq.irrep) then
                              gai(ie,ll) = gai(ie,ll)+4.0d0*gg*(
     +                                   - z(ipo(mm,kk))*z(ipo(ii,jj))
     +                                   - y(ipo(mm,kk))*y(ipo(ii,jj))
     +                                   + z(ipo(mm,kk))*y(ipo(ii,jj))
     +                                   + y(ipo(mm,kk))*z(ipo(ii,jj)))
                           endif
                           if (isej.eq.irrep.and.isil.eq.irrep) then
                              gai(ie,kk) = gai(ie,kk)+2.0d0*gg*(
     +                                   - z(ipo(mm,jj))*y(ipo(ii,ll))
     +                                   - y(ipo(mm,jj))*z(ipo(ii,ll)))
                           endif
                           if (isej.eq.irrep.and.isik.eq.irrep) then
                              gai(ie,ll) = gai(ie,ll)+2.0d0*gg*(
     +                                   - z(ipo(mm,jj))*y(ipo(ii,kk))
     +                                   - y(ipo(mm,jj))*z(ipo(ii,kk)))
                           endif
                           if (isel.eq.irrep.and.isik.eq.irrep) then
                              gai(ie,jj) = gai(ie,jj)+2.0d0*gg*(
     +                                   + z(ipo(mm,ll))*z(ipo(ii,kk))
     +                                   + y(ipo(mm,ll))*y(ipo(ii,kk)))
                           endif
                           if (isek.eq.irrep.and.isil.eq.irrep) then
                              gai(ie,jj) = gai(ie,jj)+2.0d0*gg*(
     +                                   + z(ipo(mm,kk))*z(ipo(ii,ll))
     +                                   + y(ipo(mm,kk))*y(ipo(ii,ll)))
                           endif
                        enddo
                        if (isymmo(kk).eq.isymmo(ll)) then
                           tai(ii-na,jj)=tai(ii-na,jj)
     +                              +8.0d0*gg*gij(ipgij(kk,ll))
                        endif
                        if (isymmo(jj).eq.isymmo(ll)) then
                           tai(ii-na,kk)=tai(ii-na,kk)
     +                              -2.0d0*gg*gij(ipgij(jj,ll))
                        endif
                        if (isymmo(jj).eq.isymmo(kk)) then
                           tai(ii-na,ll)=tai(ii-na,ll)
     +                              -2.0d0*gg*gij(ipgij(jj,kk))
                        endif
                     else
c...                    [aj|c*]
                        if (ll.le.na) then
c...                       [aj|cl] (not used in I_ai or T_ai so skip)
                        else
c...                       [aj|cd]
                           gg = gin(iint)
                           if (kk.eq.ll) gg = 0.5d0*gg
                           do mm = 1, na
                              iskm = iperm(isymmo(kk),isymmo(mm))
                              isij = iperm(isymmo(ii),isymmo(jj))
                              isim = iperm(isymmo(ii),isymmo(mm))
                              iskj = iperm(isymmo(kk),isymmo(jj))
                              islj = iperm(isymmo(ll),isymmo(jj))
                              islm = iperm(isymmo(ll),isymmo(mm))
                              if (iskm.eq.irrep.and.isij.eq.irrep) then
                                 tai(ll-na,mm)=tai(ll-na,mm)+4.0d0*gg*(
     +                                 + z(ipo(kk,mm))*z(ipo(ii,jj))
     +                                 + y(ipo(kk,mm))*y(ipo(ii,jj))
     +                                 - z(ipo(kk,mm))*y(ipo(ii,jj))
     +                                 - y(ipo(kk,mm))*z(ipo(ii,jj)))
                              endif
                              if (islm.eq.irrep.and.isij.eq.irrep) then
                                 tai(kk-na,mm)=tai(kk-na,mm)+4.0d0*gg*(
     +                                 + z(ipo(ll,mm))*z(ipo(ii,jj))
     +                                 + y(ipo(ll,mm))*y(ipo(ii,jj))
     +                                 - z(ipo(ll,mm))*y(ipo(ii,jj))
     +                                 - y(ipo(ll,mm))*z(ipo(ii,jj)))
                              endif
                              if (isim.eq.irrep.and.iskj.eq.irrep) then
                                 tai(ll-na,mm)=tai(ll-na,mm)+2.0d0*gg*(
     +                                 + z(ipo(ii,mm))*y(ipo(kk,jj))
     +                                 + y(ipo(ii,mm))*z(ipo(kk,jj)))
                              endif
                              if (isim.eq.irrep.and.islj.eq.irrep) then
                                 tai(kk-na,mm)=tai(kk-na,mm)+2.0d0*gg*(
     +                                 + z(ipo(ii,mm))*y(ipo(ll,jj))
     +                                 + y(ipo(ii,mm))*z(ipo(ll,jj)))
                              endif
                              if (iskm.eq.irrep.and.islj.eq.irrep) then
                                 tai(ii-na,mm)=tai(ii-na,mm)+2.0d0*gg*(
     +                                 - z(ipo(kk,mm))*z(ipo(ll,jj))
     +                                 - y(ipo(kk,mm))*y(ipo(ll,jj)))
                              endif
                              if (islm.eq.irrep.and.iskj.eq.irrep) then
                                 tai(ii-na,mm)=tai(ii-na,mm)+2.0d0*gg*(
     +                                 - z(ipo(ll,mm))*z(ipo(kk,jj))
     +                                 - y(ipo(ll,mm))*y(ipo(kk,jj)))
                              endif
                           enddo
                           if (isymmo(kk).eq.isymmo(ll)) then
                              tai(ii-na,jj)=tai(ii-na,jj)
     +                              +8.0d0*gg*gab(ipgab(kk-na,ll-na))
                           endif
                           if (isymmo(ii).eq.isymmo(kk)) then
                               tai(ll-na,jj)=tai(ll-na,jj)
     +                              -2.0d0*gg*gab(ipgab(ii-na,kk-na))
                           endif
                           if (isymmo(ii).eq.isymmo(ll)) then
                               tai(kk-na,jj)=tai(kk-na,jj)
     +                              -2.0d0*gg*gab(ipgab(ii-na,ll-na))
                           endif
                        endif
                     endif
                  else
c...                 [ab|**]
                     if (kk.le.na) then
c...                    [ab|kl] (not used in I_ai or T_ai so skip)
                     else
c...                    [ab|c*]
                        if (ll.le.na) then
c...                       [ab|cl]
                           gg = gin(iint)
                           if (ii.eq.jj) gg = 0.5d0*gg
                           do mm = 1, na
                              isim = iperm(isymmo(ii),isymmo(mm))
                              iskl = iperm(isymmo(kk),isymmo(ll))
                              iskm = iperm(isymmo(kk),isymmo(mm))
                              isil = iperm(isymmo(ii),isymmo(ll))
                              isjl = iperm(isymmo(jj),isymmo(ll))
                              isjm = iperm(isymmo(jj),isymmo(mm))
                              if (isim.eq.irrep.and.iskl.eq.irrep) then
                                 tai(jj-na,mm)=tai(jj-na,mm)+4.0d0*gg*(
     +                                 + z(ipo(ii,mm))*z(ipo(kk,ll))
     +                                 + y(ipo(ii,mm))*y(ipo(kk,ll))
     +                                 - z(ipo(ii,mm))*y(ipo(kk,ll))
     +                                 - y(ipo(ii,mm))*z(ipo(kk,ll)))
                              endif
                              if (isjm.eq.irrep.and.iskl.eq.irrep) then
                                 tai(ii-na,mm)=tai(ii-na,mm)+4.0d0*gg*(
     +                                 + z(ipo(jj,mm))*z(ipo(kk,ll))
     +                                 + y(ipo(jj,mm))*y(ipo(kk,ll))
     +                                 - z(ipo(jj,mm))*y(ipo(kk,ll))
     +                                 - y(ipo(jj,mm))*z(ipo(kk,ll)))
                              endif
                              if (iskm.eq.irrep.and.isil.eq.irrep) then
                                 tai(jj-na,mm)=tai(jj-na,mm)+2.0d0*gg*(
     +                                 + z(ipo(kk,mm))*y(ipo(ii,ll))
     +                                 + y(ipo(kk,mm))*z(ipo(ii,ll)))
                              endif
                              if (iskm.eq.irrep.and.isjl.eq.irrep) then
                                 tai(ii-na,mm)=tai(ii-na,mm)+2.0d0*gg*(
     +                                 + z(ipo(kk,mm))*y(ipo(jj,ll))
     +                                 + y(ipo(kk,mm))*z(ipo(jj,ll)))
                              endif
                              if (isim.eq.irrep.and.isjl.eq.irrep) then
                                 tai(kk-na,mm)=tai(kk-na,mm)+2.0d0*gg*(
     +                                 - z(ipo(ii,mm))*z(ipo(jj,ll))
     +                                 - y(ipo(ii,mm))*y(ipo(jj,ll)))
                              endif
                              if (isjm.eq.irrep.and.isil.eq.irrep) then
                                 tai(kk-na,mm)=tai(kk-na,mm)+2.0d0*gg*(
     +                                 - z(ipo(jj,mm))*z(ipo(ii,ll))
     +                                 - y(ipo(jj,mm))*y(ipo(ii,ll)))
                              endif
                           enddo
                           if (isymmo(ii).eq.isymmo(jj)) then
                              tai(kk-na,ll)=tai(kk-na,ll)
     +                                 +8.0d0*gg*gab(ipgab(ii-na,jj-na))
                           endif
                           if (isymmo(ii).eq.isymmo(kk)) then
                              tai(jj-na,ll)=tai(jj-na,ll)
     +                                 -2.0d0*gg*gab(ipgab(ii-na,kk-na))
                           endif
                           if (isymmo(jj).eq.isymmo(kk)) then
                              tai(ii-na,ll)=tai(ii-na,ll)
     +                                 -2.0d0*gg*gab(ipgab(jj-na,kk-na))
                           endif
                        else
c...                       [ab|cd] (not used in I_ai or T_ai so skip)
                        endif
                     endif
                  endif
               endif
            enddo
            if (l6blki.lt.0) go to 10
         endif
      enddo
c
      do i = 1, na
         do ie = 1, nvirt
            tai(ie,i) = tai(ie,i) + gai(ie,i)
         enddo
      enddo  
      end
c
c-----------------------------------------------------------------------
c
      subroutine z_y_su(ipo,nbasis,irrep,z,y,cscf,zsu,ysu)
      implicit none
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      integer isym,isize,neig,iperm,isspac,maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
INCLUDE(common/infoa)
c
c...  Parameters
c
      integer nbasis,ipo(nbasis,nbasis),irrep
      REAL z(*), y(*), cscf(nbasis,*)
      REAL zsu(nbasis,nbasis), ysu(nbasis,nbasis)
c
c...  Local variables
c
      integer iu, is, ii, ia
c
      do iu = 1, nbasis
         do is = 1, nbasis
            zsu(is,iu) = 0.0d0
            ysu(is,iu) = 0.0d0
            do ii = 1, na
               do ia = na+1, num
                  if (iperm(isymmo(ii),isymmo(ia)).eq.irrep) then
                     zsu(is,iu) = zsu(is,iu) + 
     +                  cscf(is,ia)*cscf(iu,ii)*z(ipo(ii,ia))
                     ysu(is,iu) = ysu(is,iu) + 
     +                  cscf(is,ia)*cscf(iu,ii)*y(ipo(ii,ia))
                  endif
               enddo
            enddo
         enddo
      enddo
      end
c
c-----------------------------------------------------------------------
c
      subroutine rpa_grad(q)
      implicit none
      REAL q(*)
c
c...  Evaluates expression (57) in Ortiz.
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/rpacom)
INCLUDE(common/rpadcom)
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
INCLUDE(common/infoa)
c
c...  Locate variables:
c
c...  Counters
c
c     integer i
c
c...  Sizes of datastructures
c
      integer mnij,mnab
c
c...  Addresses
c
      integer igij, ipgij, iqxij, igab, ipgab, iqxab, isx, iij, iab, iai
      integer itai, iuai
c
c...  Functions
c
      integer igmem_alloc, lenint
c
      mnij   = na*na
      mnab   = nvirt*nvirt
c
c...  Sum(ij, Gamma_ij Q(x)_ij)
c
c     isizegij = 0
c     do i = 1, nirr
c        isizegij = isizegij+(kcore(i)+1)*kcore(i)/2
c     enddo
      igij  = igmem_alloc(isizegij)
      ipgij = igmem_alloc(lenint(na*na))
      iqxij = igmem_alloc(mnij)
      call add_gq_ij(q(ipgij),q(igij),q(iqxij))
      call gmem_free(iqxij)
      call gmem_free(ipgij)
      call gmem_free(igij)
c
c...  Sum(ab, Gamma_ab Q(x)_ab)
c
      igab  = igmem_alloc(isizegab)
      ipgab = igmem_alloc(lenint(nvirt*nvirt))
      iqxab = igmem_alloc(mnab)
      call add_gq_ab(q(ipgab),q(igab),q(iqxab))
      call gmem_free(iqxab)
      call gmem_free(ipgab)
      call gmem_free(igab)
c
c...  Sum(ij,S(x)_ij I_ij)+Sum(ab,S(x)_ab I_ab)+Sum(ai,S(x)_ai I_ai)
c
      isx = igmem_alloc(num*(num+1)/2 )
      iij = igmem_alloc(na*na)
      iab = igmem_alloc(nvirt*nvirt)
      iai = igmem_alloc(nvirt*na)
      call add_si(q(isx),q(iij),q(iab),q(iai))
      call gmem_free(iai)
      call gmem_free(iab)
      call gmem_free(iij)
      call gmem_free(isx)
c
c...  Sum(ai,U(x)_ai T_ai)
c
      isx  = igmem_alloc(num*(num+1)/2 )
      itai = igmem_alloc(nvirt*na)
      iuai = igmem_alloc(nvirt*na)
      call add_ut(q(iuai),q(itai),q(isx))
      call gmem_free(iuai)
      call gmem_free(itai)
      call gmem_free(isx)
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine add_gq_ij(ipgij,gij,qxij)
      implicit none
c
c...  Sum(ij, Gamma_ij Q(x)_ij) Was in Ortiz but wrong!!!
c
c...  Has to be
c
c     Sum(i, Gamma_ii Q(x)_ii)
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/rpadcom)
INCLUDE(common/sector)
INCLUDE(common/infoa)
INCLUDE(common/restar)
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
c
c...  Parameters
c
      integer ipgij(na,na)
      REAL gij(*), qxij(na,na)
c
c...  Local variables
c
      integer irpasec1, irpatyp1, irpablk1
      parameter(irpasec1=352,irpatyp1=121)
      integer irpasec2, irpatyp2, irpablk2
      parameter(irpasec2=353,irpatyp2=122)
      integer nwords(8), lensect, nstates, icoord, isym, istate, iat
      integer ii, jj, iblkb
      REAL rwords
      equivalence(rwords,nwords(1))
c
c...  Functions
c
      integer lenint, lensec
c
      call secget(irpasec2,irpatyp2,irpablk2)
      call rdedx(ipgij,lenint(na*na),irpablk2,numdu)
      irpablk2 = irpablk2 + lensec(lenint(na*na))
c
      call secget(irpasec1,irpatyp1,irpablk1)
      lensect = 0
      nstates = 0
      do isym = 1, nirr
         do istate = 1, rpad_nstates(isym)
            nstates = nstates + 1
            call rdedx(nwords,lenint(8),irpablk1+lensect,numdu)
            lensect = lensect+1
            call rdedx(gij,nwords(1),irpablk1+lensect,numdu)
            lensect = lensect+lensec(nwords(1))+lensec(nwords(2))+
     +         lensec(nwords(3))+lensec(nwords(4))+lensec(nwords(5))+
     +         lensec(nwords(6))+lensec(nwords(7))+lensec(nwords(8))
            iblkb = iblks
            do iat = 1, nat
               do icoord = 1, 3
                  call rdedx(qxij,na*na,iblkb,ifils)
                  iblkb = iblkb+lensec(na*na)
                  do jj = 1, na
                     do ii = 1, jj-1
                        if (isymmo(ii).eq.isymmo(jj)) then
                           rpade(icoord,iat,nstates) = 
     +                        rpade(icoord,iat,nstates) + 
     +                        gij(ipgij(ii,jj))*qxij(ii,jj) +
     +                        gij(ipgij(ii,jj))*qxij(jj,ii)
                        endif
                     enddo
                     rpade(icoord,iat,nstates) = 
     +                  rpade(icoord,iat,nstates) + 
     +                  gij(ipgij(jj,jj))*qxij(jj,jj) 
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine add_gq_ab(ipgab,gab,qxab)
      implicit none
c
c...  Sum(ab, Gamma_ab Q(x)_ab) Was in Ortiz but wrong !!!
c
c     Has to be:
c
c...  Sum(a, Gamma_aa Q(x)_aa)
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/rpadcom)
INCLUDE(common/sector)
INCLUDE(common/restar)
INCLUDE(common/infoa)
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
c
c...  Parameters
c
      integer ipgab(nvirt,nvirt)
      REAL gab(*), qxab(nvirt,nvirt)
c
c...  Local variables
c
      integer irpasec1, irpatyp1, irpablk1
      parameter(irpasec1=352,irpatyp1=121)
      integer irpasec2, irpatyp2, irpablk2
      parameter(irpasec2=353,irpatyp2=122)
      integer nwords(8), lensect, nstates, icoord, isym, istate, iat
      integer ii, jj, iblkb
      REAL rwords
      equivalence(rwords,nwords(1))
c
c...  Functions
c
      integer lenint, lensec
c
      call secget(irpasec2,irpatyp2,irpablk2)
      irpablk2 = irpablk2 + lensec(lenint(na*na))
      call rdedx(ipgab,lenint(nvirt*nvirt),irpablk2,numdu)
c
      call secget(irpasec1,irpatyp1,irpablk1)
      lensect = 0
      nstates = 0
      do isym = 1, nirr
         do istate = 1, rpad_nstates(isym)
            nstates = nstates + 1
            call rdedx(nwords,lenint(8),irpablk1+lensect,numdu)
            lensect = lensect+1+lensec(nwords(1))
            call rdedx(gab,nwords(2),irpablk1+lensect,numdu)
            lensect = lensect+lensec(nwords(2))+
     +         lensec(nwords(3))+lensec(nwords(4))+lensec(nwords(5))+
     +         lensec(nwords(6))+lensec(nwords(7))+lensec(nwords(8))
            iblkb = iblks + 3*nat*lensec(na*na)
            do iat = 1, nat
               do icoord = 1, 3
                  call rdedx(qxab,nvirt*nvirt,iblkb,ifils)
                  iblkb = iblkb+lensec(nvirt*nvirt)
                  do jj = 1, nvirt
                     do ii = 1, jj-1
                        if (isymmo(na+ii).eq.isymmo(na+jj)) then
                           rpade(icoord,iat,nstates) = 
     +                        rpade(icoord,iat,nstates) + 
     +                        gab(ipgab(ii,jj))*qxab(ii,jj) +
     +                        gab(ipgab(ii,jj))*qxab(jj,ii)
                        endif
                     enddo
                     rpade(icoord,iat,nstates) = 
     +                  rpade(icoord,iat,nstates) + 
     +                  gab(ipgab(jj,jj))*qxab(jj,jj) 
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine add_si(s,iij,iab,iai)
      implicit none
c
c...  Sum(ij,S(x)_ij I_ij)+Sum(ab,S(x)_ab I_ab)+Sum(ai,S(x)_ai I_ai)
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/sector)
INCLUDE(common/rpadcom)
INCLUDE(common/cndx41)
INCLUDE(common/restar)
      integer nvirtx,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirtx,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
c
c...  Parameters
c
      REAL s(*), iij(na,na), iab(nvirt,nvirt), iai(nvirt,na)
c
c...  Local variables
c
      integer irpasec1, irpatyp1, irpablk1
      parameter(irpasec1=352,irpatyp1=121)
      integer lensect, nstates, ns, nwords(8), ibs, iat, icoord, isym
      integer istate, ii, jj
      REAL rwords
      equivalence(rwords,nwords(1))
c
c...  Functions
c
      integer lenint, lensec
c
      call secget(irpasec1,irpatyp1,irpablk1)
      lensect = 0
      nstates = 0
      ns      = num*(num+1)/2
      do isym = 1, nirr
         do istate = 1, rpad_nstates(isym)
            nstates = nstates + 1
            call rdedx(nwords,lenint(8),irpablk1+lensect,numdu)
            lensect = lensect+1+lensec(nwords(1))+lensec(nwords(2))
            call rdedx(iij,nwords(3),irpablk1+lensect,numdu)
            lensect = lensect+lensec(nwords(3))
            call rdedx(iab,nwords(4),irpablk1+lensect,numdu)
            lensect = lensect+lensec(nwords(4))
            call rdedx(iai,nwords(5),irpablk1+lensect,numdu)
            lensect = lensect+lensec(nwords(5))+
     +         lensec(nwords(6))+lensec(nwords(7))+lensec(nwords(8))
            ibs = iochf(14)
            do iat = 1, nat
               do icoord = 1, 3
                  call rdedx(s,ns,ibs,ifockf)
                  ibs = ibs+lensec(ns)
                  do jj = 1, na
                     do ii = 1, jj-1
                        rpade(icoord,iat,nstates) = 
     +                     rpade(icoord,iat,nstates) + 
     +                     s(jj*(jj-1)/2+ii)*iij(ii,jj) +
     +                     s(jj*(jj-1)/2+ii)*iij(jj,ii) 
                     enddo
                     rpade(icoord,iat,nstates) =
     +                  rpade(icoord,iat,nstates) +
     +                  s(jj*(jj-1)/2+jj)*iij(jj,jj) 
                  enddo
                  do jj = 1, nvirtx
                     do ii = 1, jj-1
                        rpade(icoord,iat,nstates) = 
     +                     rpade(icoord,iat,nstates) + 
     +                     s((na+jj)*(na+jj-1)/2+na+ii)*
     +                     iab(ii,jj) +
     +                     s((na+jj)*(na+jj-1)/2+na+ii)*
     +                     iab(jj,ii) 
                     enddo
                     rpade(icoord,iat,nstates) =
     +                  rpade(icoord,iat,nstates) +
     +                  s((na+jj)*(na+jj-1)/2+na+jj)*
     +                  iab(jj,jj) 
                  enddo
                  do jj = 1, na
                     do ii = 1, nvirtx
                        rpade(icoord,iat,nstates) = 
     +                     rpade(icoord,iat,nstates) + 
     +                     s((na+ii)*(na+ii-1)/2+jj)*iai(ii,jj)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine add_ut(uai,tai,s)
      implicit none
c
c...  Sum(ij,U(x)_ai T_ai)
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/rpadcom)
INCLUDE(common/infoa)
INCLUDE(common/sector)
INCLUDE(common/cndx41)
INCLUDE(common/restar)
      integer nvirtx,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirtx,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
c
c...  Parameters
c
      REAL uai(na,nvirt), tai(nvirt,na), s(num*(num+1)/2)
c
c...  Local variables
c
      integer irpasec1, irpatyp1, irpablk1
      parameter(irpasec1=352,irpatyp1=121)
      integer lensect, nstates, ns, nwords(8), iat, icoord
      integer isym, istate, iblku, ii, jj, ibs
      REAL rwords
      equivalence(rwords,nwords(1))
c
c...  Functions
c
      integer lenint, lensec
c
      call secget(irpasec1,irpatyp1,irpablk1)
      lensect = 0
      nstates = 0
      ns      = num*(num+1)/2
      do isym = 1, nirr
         do istate = 1, rpad_nstates(isym)
            nstates = nstates + 1
            call rdedx(nwords,lenint(8),irpablk1+lensect,numdu)
            lensect = lensect+1+lensec(nwords(1))+lensec(nwords(2))+
     +         lensec(nwords(3))+lensec(nwords(4))+lensec(nwords(5))
            call rdedx(tai,nwords(6),irpablk1+lensect,numdu)
            lensect = lensect+lensec(nwords(6))+lensec(nwords(7))+
     +         lensec(nwords(8))
            iblku = iblks + 3*nat*lensec(na*na)
     +                    + 3*nat*lensec(nvirt*nvirt)
     +                    + 3*nat*lensec(nvirt*na)
            ibs = iochf(14)
            do iat = 1, nat
               do icoord = 1, 3
                  call rdedx(uai,nvirt*na,iblku,ifils)
                  iblku = iblku+lensec(nvirt*na)
                  call rdedx(s,ns,ibs,ifockf)
                  ibs = ibs+lensec(ns)
                  do ii = 1, nvirt
                     do jj = 1, na
                        rpade(icoord,iat,nstates) = 
     +                     rpade(icoord,iat,nstates) + 
     +                     uai(jj,ii)*tai(ii,jj)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine pr_rpa_grad()
      implicit none
c
c...  Evaluates expression (57) in Ortiz.
c
c     Huub van Dam, 1999
c
INCLUDE(common/sizes)
INCLUDE(common/rpacom)
INCLUDE(common/rpadcom)
INCLUDE(common/grad2)
INCLUDE(common/infoa)
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
c
c...  Local variables
c
      integer ist, isym, istate, iat, i
c
      ist = 0
      do isym = 1, nirr
         do istate = 1, rpad_nstates(isym)
            ist = ist + 1
            write(*,*)'*** Gradient RPA ',isym,rpad_istates(istate,isym)
            do iat = 1, nat
               write(*,'(i3,3f16.8)')iat,(rpade(i,iat,ist),i=1,3)
            enddo
         enddo
      enddo
c
      ist = 0
      do isym = 1, nirr
         do istate = 1, rpad_nstates(isym)
            ist = ist + 1
            write(*,*)'*** Total Gradient RPA ',isym,
     +                rpad_istates(istate,isym)
            do iat = 1, nat
               write(*,'(i3,3f16.8)')iat,
     +                 (rpade(i,iat,ist)+de(i,iat),i=1,3)
            enddo
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine drpa_egrad()
      implicit none
c
c...  Evaluates expression (57) in Ortiz and stores the resulting
c...  total gradient in common/funct/egrad for use in the internal
c...  coordinate geometry optimisation (subroutine minit).
c
c     Huub van Dam, 2000
c
INCLUDE(common/sizes)
INCLUDE(common/rpacom)
INCLUDE(common/rpadcom)
INCLUDE(common/grad2)
INCLUDE(common/infoa)
INCLUDE(common/funct)
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
c
c...  Local variables
c
      integer iat, i
c
      if (nrpastate.ne.1) then
c
c        We can only sensibly perform a geometry optimisation on 
c        1 state at a time
c
         call caserr(
     +        'Geometry optimisation of 1 state at a time please')
      endif
c
      do iat = 1, nat
         do i = 1, 3
            egrad((iat-1)*3+i)=rpade(i,iat,1)+de(i,iat)
         enddo
      enddo
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine drpa_enrgy()
      implicit none
c
c...  Evaluates the RPA total energy of the excited state. This is 
c...  needed in some geometry optimisation routines. The energy gets
c...  stored in common/funct/enrgy
c
c     Huub van Dam, 2000
c
INCLUDE(common/sizes)
INCLUDE(common/rpacom)
INCLUDE(common/rpadcom)
INCLUDE(common/grad2)
INCLUDE(common/infoa)
INCLUDE(common/funct)
INCLUDE(common/sector)
INCLUDE(common/restri)
INCLUDE(common/iofile)
      integer icount,kcore,kvirt,ksymm
      integer lcore,lvirt,iofset,npair
      integer nbox,nbox0,jstart,iofsym
      common/craypk/icount(8),kcore(8),kvirt(8),ksymm(8),
     +              lcore(8),lvirt(8),iofset(8,8),npair(8),
     +              nbox(8),nbox0(8),jstart(8),iofsym(8)
      integer nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri
      integer nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s
      integer inozer,nozer,idone
      common/dims  /nvirt,nvc,maxr2,nvc2,nc1,nbasq,nvcnev,nbatri,
     +              nstart,nvcnst,nevlo1,nirr,nirs,maxrs,maxr2s,
     +              inozer,nozer(3),idone
c
c...  Functions
c
      integer lensec, lenint
c
c...  Local variables
c
      integer iat, i
      integer ir, ist
      REAL epole(2), ehf(10)
c
      integer irpatyp, irpasec,irpablk, nrpalen
      parameter(irpatyp=120)
      parameter(irpasec=351)
c
      integer ihftyp,ihfblk
      parameter(ihftyp=16)
c
      if (nrpastate.ne.1) then
c
c        We can only sensibly perform a geometry optimisation on 
c        1 state at a time
c
         call caserr(
     +        'Geometry optimisation of 1 state at a time please')
      endif

      call secget(irpasec,irpatyp,irpablk)
      nrpalen = lensec(lenint(nbasq))
      do ir = 1, nirr
         do ist = 0, nev(ir)-1
            if (rpad_nstates(ir).eq.1.and.
     +          rpad_istates(1,ir).eq.ist+nevlo(ir)) then
               call rdedx(epole,2,irpablk+nrpalen,numdu)
            else
               nrpalen = nrpalen + lensec(2) + 2*lensec(icount(ir))
            endif
         enddo
      enddo
c
      call secget(isect(494),ihftyp,ihfblk)
      call rdedx(ehf,10,ihfblk,idaf)
c
      enrgy = ehf(3)+epole(1)
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine prZai(ipo,z,zsq,nbas,irrep)
      implicit none
      integer nbas, ipo(nbas,nbas), nocc, irrep
      REAL z(*), zsq(nbas,nbas)
c
c     Print routine for printing RPA vectors 
c     Purpose is to facilitate debugging.
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      integer isym,isize,neig,iperm,isspac,maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
c
c...  Local variables
c
      integer ii, jj, k, ia, ib
c
      call vclr(zsq,1,nbas*nbas)
      do ii = 1, na
         do ia = na+1, nbas
            if (iperm(isymmo(ii),isymmo(ia)).eq.irrep) then
               zsq(ii,ia) = z(ipo(ii,ia))
            endif
         enddo
      enddo
      call prsqm(zsq,nbas,nbas,nbas,6)
      end
c
c-----------------------------------------------------------------------
c
      subroutine prGij(ipgij,gij,gsq,nocc)
      implicit none
      integer nocc, ipgij(nocc,nocc)
      REAL gij(*), gsq(nocc,nocc)
c
c     Print routine for Gij to facilitate debugging
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      integer isym,isize,neig,iperm,isspac,maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
c
c...  Local variables
c
      integer ii, jj, k, ia, ib
c
      call vclr(gsq,1,nocc*nocc)
      do ii = 1, nocc
         do jj = 1, nocc
            if (isymmo(ii).eq.isymmo(jj)) then
               gsq(ii,jj) = gij(ipgij(ii,jj))
               gsq(jj,ii) = gij(ipgij(jj,ii))
            endif
         enddo
      enddo
      call prsqm(gsq,nocc,nocc,nocc,6)
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine prGab(ipgab,gab,gsq,nocc)
      implicit none
c
c     Print routine for Gab to facilitate debugging.
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
c
      integer nocc, ipgab(num-nocc,num-nocc)
      REAL gab(*), gsq(num-nocc,num-nocc)
c
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      integer isym,isize,neig,iperm,isspac,maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
c
c...  Local variables
c
      integer ii, jj, k, ia, ib, nvirt
c
      nvirt = num - nocc
      call vclr(gsq,1,nvirt*nvirt)
      do ia = 1, nvirt
         do ib = 1, nvirt
            if (isymmo(nocc+ia).eq.isymmo(nocc+ib)) then
               gsq(ia,ib) = gab(ipgab(ia,ib))
               gsq(ib,ia) = gab(ipgab(ib,ia))
            endif
         enddo
      enddo
      call prsqm(gsq,nvirt,nvirt,nvirt,6)
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine rpa_energy(ipo,irrep,z,y,eps)
      implicit none
c
c     This routine evaluates the RPA energy. It can be used to check 
c     that we have assigned the correct meaning to the data structures.
c
INCLUDE(common/sizes)
      integer mmmm,isymao,isymmo
      common/craypm/mmmm(65),isymao(maxorb),isymmo(maxorb)
      integer isym,isize,neig,iperm,isspac,maxvec,mvec,mymaxv
      common/symmet/isym,isize,neig,iperm(8,8),isspac(10),maxvec,
     +              mvec,mymaxv
INCLUDE(common/infoa)
INCLUDE(common/atmblk)
      REAL gin
      common/blkin /gin(511)
      integer i205
      common/junk2 /i205(4,340)
      integer nprint,n6file,n6tape,n6blk,n6last
      common/restar/nprint(142),n6file,n6tape(20),n6blk(20),n6last(20)
      integer mword
      equivalence(mword,gin(511))
c
      integer irrep,ipo(num,num)
      REAL eps(num), z(*), y(*), exc
c
      integer ned,ifile,nw,iint,ii,jj,kk,ll,ie,im,i
      integer l6blki,l6blk(20)
      logical isij, iskl, isil, isjk, isik, isjl
      REAL gg
_IF(linux)
      external fget
_ENDIF
c
      exc = 0.0d0
c
      do i = 1, na
         do ie = na+1, num
            if (iperm(isymmo(i),isymmo(ie)).eq.irrep) then
               exc = exc + (eps(ie)-eps(i))*(
     +               z(ipo(ie,i))*z(ipo(ie,i))
     +             + y(ipo(ie,i))*y(ipo(ie,i))  )
            endif
         enddo
      enddo
c
      do i = 1, n6file
         l6blk(i) = n6blk(i) - n6last(i)
      enddo
      call setsto(1360,0,i205)
      do ifile = 1, n6file
         ned = n6tape(ifile)
         call search(n6blk(ifile),ned)
         call find(ned)
         l6blki = l6blk(ifile)
10       l6blki = l6blki + 1
_IF(c)
         call get(gin,nw)
_ELSE
         call fget(gin,nw,ned)
_ENDIF
         if (nw.gt.0) then
            call unpack(gin(num2e+1),lab816,i205,numlab)
_IF(c)
            if (l6blki.ne.0) call find(ned)
_ENDIF
            do iint = 1, mword
_IF(littleendian)
               jj = i205(1,iint)
               ii = i205(2,iint)
               ll = i205(3,iint)
               kk = i205(4,iint)
_ELSE
               ii = i205(1,iint)
               jj = i205(2,iint)
               kk = i205(3,iint)
               ll = i205(4,iint)
_ENDIF
            if (ii.le.na) then
c              (ij|kl) not used here
            else
               if (jj.le.na) then
c                 (aj|**)
                  if (kk.gt.na.and.ll.le.na) then
c                    (aj|cl)
                     gg = gin(iint)
                     if (ii.eq.kk.and.jj.eq.ll) gg = 0.5d0*gg
                     isij = iperm(isymmo(ii),isymmo(jj)).eq.irrep
                     iskl = iperm(isymmo(kk),isymmo(ll)).eq.irrep
                     isil = iperm(isymmo(ii),isymmo(ll)).eq.irrep
                     isjk = iperm(isymmo(jj),isymmo(kk)).eq.irrep
                     if (isij.and.iskl) then
                        exc = exc + 4.0d0*gg*(
     +                        z(ipo(ii,jj))*z(ipo(kk,ll))
     +                      + y(ipo(ii,jj))*y(ipo(kk,ll))
     +                      - z(ipo(ii,jj))*y(ipo(kk,ll))
     +                      - y(ipo(ii,jj))*z(ipo(kk,ll))  )
                     endif
                     if (isil.and.isjk) then
                        exc = exc + 2.0d0*gg*(
     +                        z(ipo(ii,ll))*y(ipo(kk,jj))
     +                      + y(ipo(ii,ll))*z(ipo(kk,jj))  )
                     endif
                  else
c                    (aj|kl) or (aj|cd) not used here
                  endif
               else
c                 (ab|**)
                  if (kk.le.na) then
c                    (ab|kl)
                     gg = gin(iint)
                     if (ii.eq.jj) gg = 0.5d0*gg
                     if (kk.eq.ll) gg = 0.5d0*gg
                     isik = iperm(isymmo(ii),isymmo(kk)).eq.irrep
                     isjl = iperm(isymmo(jj),isymmo(ll)).eq.irrep
                     isil = iperm(isymmo(ii),isymmo(ll)).eq.irrep
                     isjk = iperm(isymmo(jj),isymmo(kk)).eq.irrep
                     if (isik.and.isjl) then
                        exc = exc - 2.0d0*gg*(
     +                        z(ipo(ii,kk))*z(ipo(jj,ll))
     +                      + y(ipo(ii,kk))*y(ipo(jj,ll))  )
                     endif
                     if (isil.and.isjk) then
                        exc = exc - 2.0d0*gg*(
     +                        z(ipo(ii,ll))*z(ipo(jj,kk))
     +                      + y(ipo(ii,ll))*y(ipo(jj,kk))  )
                     endif
                  else
c                    (ab|cl) or (ab|cd) not used here
                  endif
               endif
            endif
            enddo
            if (l6blki.lt.0) go to 10
         endif
      enddo
c
      end
_ENDIF
      subroutine ver_rpagrad(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/rpagrad.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
