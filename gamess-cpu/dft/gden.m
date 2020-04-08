c
c Code to implement the setup of a frozen density
c expressed as a linear combination of gaussian functions
c
c It is assumed that GAMESS-UK is running and that the
c energy and force contributions will be computed there
c via the normal bq terms (for point charges) and the 
c additional calls to the 3-centre integral and gradient
c integral codes (implemented here).
c
c  - issues
c      - what happens when the coulomb fit code is also
c        in use ???
c            a) hold both sets of gaussians in the same list..
c            b) set up two basis sets and compute 2 independent
c               contributins to the energy/potential
c
c The following routines are provided
c
c Used gamini(CH) -> gamess(GAM) -> chmdat2(CH) 
c
c  gden_init
c
c  gden_add_gauss(ncentr, iout, debug)
c
c    Called once, to define the exponents. 
c    It is assumed that the coordinates have already been added
c    to the atom lists
c
c  Called from within GAMESS-UK to get terms in the 
c  energy/forces
c
c  gden_energy()
c    Compute the energy contributions
c
c  gden_forces()
c    Compute the contributions, to both AO and expansion
c    centres
c
c**********************************************************************
c
c  gden_init
c
c  by this point any control settings should have been made
c
c**********************************************************************

      function gden_init(memory_fp,iout,debug)

      implicit none

INCLUDE(common/dft_parameters)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_iofile)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_intctl)
INCLUDE(common/dft_basis_api)
INCLUDE(common/ccpdft.hf77)

INCLUDE(common/mden)

      REAL memory_fp(*)
      integer iout
      logical debug
      integer allocate_memory

      integer ierror, i

      logical opg_root
      logical rhf_sw, opshell_sw

      rhf_sw = .true.
      opshell_sw = .false.
      rks_sw  = .true.

c      if (.not. active_sw) then
c        gden_init = 0
c        return
c      endif

      icut_dft  = 9
      itol_dft  = 20
      icutd_dft = 10
      itold_dft = 20

      icut_2c = 0  
      itol_2c = 0
      icutd_2c = 0
      itold_2c = 0
      icut_3c = 0 
      itol_3c = 0
      icutd_3c = 0
      itold_3c = 0
      icut_4c = 0
      itol_4c = 0
      icutd_4c = 0
      itold_4c = 0
c
c  set default tolerances for integral codes
c  if the user hasnt explicit set them
c
      if(icut_2c .eq. 0)icut_2c = icut_dft
      if(itol_2c .eq. 0)itol_2c = itol_dft
      if(icut_3c .eq. 0)icut_3c = icut_dft
      if(itol_3c .eq. 0)itol_3c = itol_dft
      if(icut_4c .eq. 0)icut_4c = icut_dft
      if(itol_4c .eq. 0)itol_4c = itol_dft
      if(icutd_2c .eq. 0)icutd_2c = icutd_dft
      if(itold_2c .eq. 0)itold_2c = itold_dft
      if(icutd_3c .eq. 0)icutd_3c = icutd_dft
      if(itold_3c .eq. 0)itold_3c = itold_dft
      if(icutd_4c .eq. 0)icutd_4c = icutd_dft
      if(itold_4c .eq. 0)itold_4c = itold_dft

      if(debug)write(iout,*)'3-centre tolerances',icut_3c, itol_3c
c
c  print control
c
      do i=1,MAX_DEBUG
        print_sw(i) = .false.
      enddo
      debug_sw = .false.
c
c  print stack
c
      print_stack_depth = 1
      current_print_level(1) = PRINT_DEFAULT
c
c output stream
c
      iwr = iout
C 
C Write out basis sets, if requested.
c
c !!! Need to set these tags more carefully 
c
      ao_bas = 1
c model density will be basis 2
      cd_bas = 2
c nuclear density will be basis 3
      nu_bas = 3

      call order_fill
      call BL_init
      call interface_gamess(rhf_sw,opshell_sw,ao_bas,iout)
C 
C Write out AO basis sets, if requested.
c
ccc   if(opg_root() .and. print_sw(DEBUG_AOBAS))then

      if(debug)then
         write(iout,*)'AO  Basis'
         ierror=BL_write_basis(ao_bas,iout)
      endif

      if(debug)write(iout,*)'3-centre tolerances end',icut_3c, itol_3c

      gden_init = 0
      return
      end

c**********************************************************************
c
c  gden_add_gauss(ncent, iout, debug)
c
c   ncent   (in)   number of sites to consider
c   iout    (in)   unit number for listing 
c   debug   (in)   debug output flag 
c
c  in common/blur
c   alpha   (in)   exponents (negative values denote no model fn
c                  on the site)
c   wght    (in)   weight - i.e. charge being blurred
c  
c  by this point any control settings should have been made
c
c**********************************************************************

      subroutine gden_add_gauss(ncent,iout,debug)

      implicit none

      integer ncent,iout
      logical debug

INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_mol_info)

INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/blur)

INCLUDE(common/dft_intctl)
INCLUDE(common/mden)

      integer ibas, i, type, atomno, angm, hybr
      integer nprm_count, itype, ierror
      REAL ex(10),cs(10),cp(10),cd(10),cf(10),cg(10)

      integer itest, idum, icount, icount2

      integer inuc, inullnuc, inullcd

c this is kept internally just for houskeeping
      REAL expon(max_atype)
c
c External functions
c
      logical opg_root

      if(debug) write(iout,*) 'Entered gden_add_gauss'
c
      ierror = BL_clear_basis_set(cd_bas)
      ierror = BL_clear_basis_set(nu_bas)

c
c Create three atom tags 
c
c inullnuc - nuclear basis tag for atoms which to 
c            not have nuclear charges (blurred centres)
c
c inuc     - a single S function for nuclear charges
c            represented as a nuclear basis
c 
c inullcd  - a CD basis for atoms with no CD basis
c            functions - all un-blurred
c            point charges and real atoms
c

      atomno = -1

      inullnuc=BL_create_atomtag(nu_bas,atomno)
      inuc=BL_create_atomtag(nu_bas,atomno)

      inullcd=BL_create_atomtag(cd_bas,atomno)

c expon is used to remove the need for duplicate
c CD tags if atoms share a blur exponent
c the null tag won't be used for this

      expon(inullcd) = -1.0

      if(debug)then
         write(iout,*)'atom tags'
         write(iout,*)'null cd tag  = ',inullcd
         write(iout,*)'nuclear tag  = ',inullcd
         write(iout,*)'null nuc tag = ',inullnuc
      endif
c
c Check in a single S basis function, we will use
c this for all nuclei
c
      nprm_count = 1
      angm = 1
      hybr = 1
      cs(nprm_count)=1.0d0
      cp(nprm_count)=1.0d0
      cd(nprm_count)=1.0d0
      cf(nprm_count)=1.0d0
      cg(nprm_count)=1.0d0
      ex(nprm_count)=100.0d0
      ierror=BL_import_shell(nu_bas,inuc,nprm_count,angm,hybr,
     &     ex,cs,cp,cd,cf,cg)

      icount = 0
      icount2 = 0

      do i = 1,ncent
c
c Is there a blurred charge on the centre
c
         if(debug)write(iout,*)'atom',i,' blur ',blexpo(i)

         if(blexpo(i).gt.-1.0d-10)then
c
c Check if this exponent has been used before
c
            type = -1
            do itype = 1,BL_num_types(cd_bas)
               if (abs(expon(itype)-blexpo(i)) .lt. 1.0d-10)then
                  type=itype
               endif
            enddo

            if(type .eq. -1)then
c
c Create new atom type (use -1 as this is not a real atom)
c
               atomno = -1
               type=BL_create_atomtag(cd_bas,atomno)
               if(debug)
     &         write(iout,*) 'creating CD function exp=',blexpo(i),type
c
c Check in a single S basis function
c Use weight as the coefficient
c
               nprm_count = 1
               angm = 1
               hybr = 1
               cs(nprm_count)=1.0d0
               cp(nprm_count)=1.0d0
               cd(nprm_count)=1.0d0
               cf(nprm_count)=1.0d0
               cg(nprm_count)=1.0d0
               ex(nprm_count)=blexpo(i)
               ierror=BL_import_shell(cd_bas,type,nprm_count,angm,hybr,
     &              ex,cs,cp,cd,cf,cg)
               expon(type) = blexpo(i)
            endif
            if(debug)write(iout,*) 'assigning cd type', type
            idum = BL_assign_type(cd_bas,i,type)
            icount = icount + 1
            bltab(icount) = i
         else
c
c assign null type 
c     
            if(debug)write(iout,*) 'assigning null cd'
            idum = BL_assign_type(cd_bas,i,inullcd)
         endif

c
c assign nuclear density basis
c
         if(ian(i) .gt. 0)then
            if(debug)write(iout,*) 'assigning nuc type',i
            idum = BL_assign_type(nu_bas,i,inuc)
            icount2 = icount2 + 1
            nutab(icount2) = i
         else
            if(debug)write(iout,*) 'assigning null nuc type',i
            idum = BL_assign_type(nu_bas,i,inullnuc)
         endif
      enddo

      write(iout,*)'blurred functions assigned=',icount
      write(iout,*)'nuclear functions assigned=',icount2

      call checkin_basis(cd_bas, debug)
      call checkin_basis(nu_bas, debug)

      if(debug)then
         write(iout,*)'Blurred charge  Basis'
         if(opg_root())ierror=BL_write_basis(cd_bas,iout)

         write(iout,*)'Nuclear charge  Basis'
         if(opg_root())ierror=BL_write_basis(nu_bas,iout)

         ierror = BL_summarise()
      endif

      return
      end

***********************************************************************
*
*  gden_energy
*
* Description:                                                        *
*                                                                     *
* Form the Kohn-Sham matrix contributions from the blurred charges    *
* Also construct contributions to the energy                          *
*
*  e_nuc ... blurred-charge... nuclei
*  e_elec .. blurred-charge... electrons
*
*
* NB - e_elec is diagnostic only as ehf1 will include this            *
* 
***********************************************************************

      integer function gden_energy(
     &     coords,
     &     kma,kmb,adenm,bdenm,e_elec,e_nuc,
     &     memory_int,memory_fp,print,debug,iout)

      implicit none
c
C In variables
c
      REAL coords(3,*)
      REAL kma(*),kmb(*),adenm(*),bdenm(*)
      logical print, debug
      integer iout
c
c Out variables
c
      REAL e_elec, e_nuc

C *Scratch space and pointers
      integer memory_int(*)
      REAL memory_fp(*)

INCLUDE(common/dft_parameters)
INCLUDE(common/dft_api)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_memory_info)
INCLUDE(common/dft_intctl)
INCLUDE(common/dft_mol_info)
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/blur)
INCLUDE(common/mden)

c to suppres output from integral routines
INCLUDE(../m4/common/restar)

C *Local variables
      integer i,j,iscr,iscrb,ltri
      REAL grad_dum
      integer xbfn_num
      integer idum, nsh, n_ao, n_cd, n_nu
      integer gout_pt,eri_pt,iso_pt,coef_pt, coef2_pt
      integer imode
      REAL dum
      REAL gamma, fac, charge
      integer i2c2e_pt
      integer nprint_keep
C *Functions
      integer allocate_memory
      logical opg_root
      REAL tracep      
      integer cd_update_geom

      nprint_keep = nprint

      if(debug)then
         write(iout,*) 'Entering gden_energy'
         write(iout,*)'3-centre tolerances',icut_3c, itol_3c
         write(iout,*)'2-centre tolerances',icut_2c, itol_2c
      else
         nprint = -5
      endif
c
c  Refresh coordinates
c
      idum = cd_update_geom(coords)

      nsh = BL_max_shell_count()
      n_ao = BL_basis_size(ao_bas)
      n_cd = BL_basis_size(cd_bas)
      n_nu = BL_basis_size(nu_bas)

      ltri=((n_ao+1)*n_ao)/2

      if(opg_root() .and. debug)then
         write(iout,*)'** Input Density **'
         call prtri(adenm,n_ao)
      endif

      iscr = allocate_memory(ltri,'d')
      call aclear_dp(memory_fp(iscr),ltri,0.0d0)

      iso_pt  = allocate_memory(nsh,'d')
      gout_pt = allocate_memory(50625,'d')
      eri_pt  = allocate_memory(225*n_cd,'d')
      coef_pt = allocate_memory(n_cd,'d')
      coef2_pt = allocate_memory(n_nu,'d')
c 
c Note -1 so that model density is stored in units of
c the electron charge
c
      do i = 1, n_cd
         gamma = blexpo(bltab(i))
         fac = (1.0d0*gamma/3.1415927d0)**1.5d0
         memory_fp(coef_pt + i -1) = - fac*blwght(bltab(i))
         if(debug)write(iout,*)'Coef ',i, memory_fp(coef_pt + i -1)
      enddo

c
      do i = 1, n_nu
c @@ pseudo alert
	 if(debug)write(iout,*)'Nutab',i,nutab(i)
         charge = ian(nutab(i))
         if(debug)write(iout,*)'ian',i,ian(nutab(i)),charge
         gamma = 100.0d0
         fac = (1.0d0*gamma/3.1415927d0)**1.5d0
         memory_fp(coef2_pt + i -1) = -fac*charge
         if(debug)write(iout,*)'Coef2 ',i, memory_fp(coef2_pt + i -1)
      enddo

      i2c2e_pt = allocate_memory(n_cd*n_nu,'d')
      imode = 6
      call jkint_dft(memory_fp(iso_pt),
     &     memory_fp(gout_pt),nsh,
     &     cd_bas, -1, nu_bas, -1, 
     &     imode,memory_fp(i2c2e_pt),dum,
     &     dum,dum,
     &     dum,dum,dum)
c

      if(debug)then
      write(iout,*)'2-centre integrals'
      do i = 1,n_cd
         write(iout,100) (memory_fp(i2c2e_pt + 
     &        n_cd*(j-1) + (i-1) ) *
     &        memory_fp(coef_pt + i -1) * 
     &        memory_fp(coef2_pt + j -1),
     &        j=1,n_nu)
      enddo
      endif

      e_nuc = 0.0d0
      do i = 1,n_cd
         do j = 1,n_nu
            e_nuc= e_nuc + 
     &           memory_fp(i2c2e_pt + n_cd*(j-1) + (i-1) ) *
     &           memory_fp(coef_pt + i-1) * 
     &           memory_fp(coef2_pt + j-1)
         enddo
      enddo

      if(debug)write(iout,*)'e_nuc',e_nuc

 100  format(1x,10f15.8)

      call free_memory(i2c2e_pt,'d')

      imode = 3
      call jkint_dft(memory_fp(iso_pt),
     &     memory_fp(gout_pt),nsh,
     &     ao_bas, ao_bas, cd_bas, -1, 
     &     imode,memory_fp(coef_pt),dum,
     &     memory_fp(iscr),dum,
     &     memory_fp(eri_pt),dum,dum)

      call free_memory(coef2_pt,'d')
      call free_memory(coef_pt,'d')
      call free_memory(eri_pt,'d')
      call free_memory(gout_pt,'d')
      call free_memory(iso_pt,'d')
c
c Compute 3C energy contribution
c
      e_elec=tracep(adenm,memory_fp(iscr),n_ao)
c
c print Vxc
c
      if(opg_root() .and. debug)then
         if(rks_sw)then
            write(iout,*)'** V coul **'
            call prtriprec(memory_fp(iscr),n_ao)
         else
            write(iout,*)'** Vxc alpha **'
            call prtri(memory_fp(iscr),n_ao)
            write(iout,*)'** Vxc beta **'
            call prtri(memory_fp(iscrb),n_ao)
         endif
      endif

c
c sum in Vxc
c
      call daxpy(ltri,1.0d0,memory_fp(iscr),1,kma,1)
      call free_memory(iscr,'d')

      if(.not. rks_sw)then
         call daxpy(ltri,1.0d0,memory_fp(iscrb),1,kmb,1)
         call free_memory(iscrb,'d')
      endif
c
c exit code if requested
c

      if(abort_sw)then
         if(opg_root())then
            write(iout,*)'EXIT REQUESTED'
            call list_time_periods(.false.,.false.)
         endif
         call pg_end
         call exitc(0)
      endif

      GDEN_energy = 0

      if(debug)call gmem_check_guards('end gden_energy')

      if(debug)then
         nprint = nprint_keep
      endif

      end


      subroutine prtriprec(d,n)
c
c     ----- print out a triangular matrix -----
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/runlab)
      dimension d(*),dd(12)
      mmax=7
      imax = 0
  100 imin = imax+1
      imax = imax+mmax
      if (imax .gt. n) imax = n
      write (iwr,9008)
      write (iwr,9028) (i,i = imin,imax)
      do 160 j = 1,n
      k = 0
      do 140 i = imin,imax
      k = k+1
      m = max(i,j)*(max(i,j)-1)/2 + min(i,j)
  140 dd(k) = d(m)
      write (iwr,9048) j,(dd(i),i = 1,k)
  160 continue
      if (imax .lt. n) go to 100
      return
 9008 format(/)
 9028 format(6x,7(6x,i3,6x))
 9048 format(i5,1x,7f15.10)
 8028 format(6x,12(3x,i3,3x))
 8048 format(i5,1x,12f9.4)
      end


C *****************************************************************************
c
c   integer function gden_forces - Driver for Force evaluation
c
c    Arguments:
c
c     coords,        I   current coordinates
c     adenm,         I   alpha or total density 
c     bdenm,         I   beta density, not used.
c     memory_int,    S   base for dynamic memory allocation (integer)
c     memory_fp,     S   base for dynamic memory allocation (real*8)
c     grad,         I/O  gradient 
c
C *****************************************************************************

      integer function gden_forces(coords,adens,bdens,
     &     memory_int,memory_fp,grad,print,debug,iout)
      implicit none

      REAL coords(3,*)
      REAL adens(*),bdens(*)
      integer memory_int(*)
      REAL memory_fp(*)
C *Out variables
      REAL grad(3,*)

      logical print, debug
      integer iout

INCLUDE(common/dft_parameters) 
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_basis_api)
INCLUDE(common/dft_api)

cccINCLUDE(common/dft_mbasis)

INCLUDE(../m4/common/symtry)      ! for iso
INCLUDE(../m4/common/iofile)      ! to access idaf

INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/blur)
INCLUDE(../m4/common/restar)
INCLUDE(common/mden)
     
c
c Local variables
c
      integer lbasf,i,j,col,ico
      REAL gamma, fac, charge, dum
      integer coef_pt, coef2_pt
      integer idum
      integer iiso, nsh, n_ao, n_cd, n_nu
      character*8 zrhf
      REAL grad2(3,maxat)
      REAL grad3(3,maxat)
      integer nprint_keep

c
c Functions
c
      integer allocate_memory
      integer cd_update_geom
      logical opg_root
      integer ipg_nodeid

      nprint_keep = nprint
      if(debug)then
         write(iout,*) 'Entering model density grad'
      else
         nprint = -5
      endif
c
C *End declarations
c
C    ===== Store current coordinates ====
c
      idum = cd_update_geom(coords)
c
C    ===== Form the gradients due to the coulomb part   =====

      n_cd=BL_basis_size(cd_bas)
      n_ao=BL_basis_size(ao_bas)
      n_nu=BL_basis_size(nu_bas)
C
C Assign atomic charges as fitting coefficients
C (probably will need a normalisation term also)
c
      coef_pt = allocate_memory(n_cd,'d')
      coef2_pt = allocate_memory(n_nu,'d')

      do i = 1, n_cd
         gamma = blexpo(bltab(i))
         fac = (1.0d0*gamma/3.1415927d0)**1.5d0
         memory_fp(coef_pt + i -1) = fac*blwght(bltab(i))
         if(debug)write(iout,*)'Coef ',i, memory_fp(coef_pt + i -1)
      enddo

      do i = 1, n_nu
c @@ pseudo alert
         charge = ian(nutab(i))
         gamma = 100.0d0
         fac = (1.0d0*gamma/3.1415927d0)**1.5d0
         memory_fp(coef2_pt + i -1) = -fac*charge
         if(debug)write(iout,*)'Coef2 ',i, memory_fp(coef2_pt + i -1)
      enddo

c
c NB iso should be redundant as symmetry is not used, but array is
c still referenced
c
      nsh = BL_max_shell_count()
      iiso = allocate_memory(nsh,'d')

      if(opg_root() .and. print_sw(DEBUG_FORCES) )then
         write(iout,*) 'Forces before blurred charge terms:'
         write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
         do i=1,natoms
            write(iout,1000)i,(grad(j,i),j=1,3)
         enddo
 1000    format(i3,3f16.8)
      endif

      do i=1,natoms
         do j = 1,3
            grad2(j,i) = 0.0d0
            grad3(j,i) = 0.0d0
         enddo
      enddo
c
c 3 centre term
c
      zrhf = 'rhf'
      call jkder_dft(zrhf,
     &     memory_fp,
     &     memory_fp(iiso),
     &     nsh,
     &     ao_bas, ao_bas, cd_bas, -1,
     &     adens,bdens,
     &     memory_fp(coef_pt),dum,
     &     grad3)

cc      if(opg_root() .and. print_sw(DEBUG_FORCES) )then
        if(debug)then
         write(iout,*) '3-centre blurred charge forces:'
         write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
         do i=1,natoms
            write(iout,1000)i,(grad3(j,i),j=1,3)
         enddo
        endif
c
c 2 centre term (blur, nuclear)
c
      zrhf = 'rhf'
      call jkder_dft(zrhf,
     &     memory_fp,
     &     memory_fp(iiso),
     &     nsh,
     &     cd_bas, -1, nu_bas, -1,
     &     adens,bdens,
     &     memory_fp(coef_pt),
     &     memory_fp(coef2_pt),
     &     grad2)

      if(debug)then
        write(iout,*) '2-centre blurred charge terms:'
        write(iout,*) 
     &        'atom         de/dx             de/y          de/dz'
        do i=1,natoms
           write(iout,1000)i,(grad2(j,i),j=1,3)
        enddo
      endif

      do i=1,natoms
         do j = 1,3
            grad(j,i) = grad(j,i) + grad2(j,i) - grad3(j,i)
c            grad(j,i) = grad(j,i) + grad2(j,i)
         enddo
      enddo

      if(debug)then
        write(iout,*) 'cumulated forces:'
        write(iout,*) 
     &     'atom         de/dx             de/y          de/dz'
        do i=1,natoms
           write(iout,1000)i,(grad(j,i),j=1,3)
       enddo
      endif

      call free_memory(iiso,'d')
      call free_memory(coef2_pt,'d')
      call free_memory(coef_pt,'d')

      gden_forces = 0

      if(debug)call gmem_check_guards('end gmem_forces')

      nprint = nprint_keep

      return
      end

      subroutine ver_dft_gden(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/gden.m,v $
     +     "/
      data revision /
     +     "$Revision: 6286 $"
     +      /
      data date /
     +     "$Date: 2013-04-05 23:36:00 +0200 (Fri, 05 Apr 2013) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
