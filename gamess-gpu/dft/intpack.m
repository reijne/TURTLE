c---- memory counting routines -----------------------------------------
      subroutine memreq_intPack_drv(memory_int,
     &                              memory_fp,
     &                              matrix_scr,
     &                              expand_sw,Vform_sw,
     &                              oe1c_int_sw,
     &                              oe2c_ove_sw,oe3c_ove_sw,
     &                              te2c_rep_sw,te3c_rep_sw,
     &                              mat_mult_sw,
     &                              matrix_out)
C **********************************************************************
C *Description                                                         *
C *Top level driver for intPack.                                       *
C **********************************************************************
      implicit none
C **********************************************************************
C *Declarations                                                        *
C *                                                                    *
C *Parameters                                                          *
INCLUDE(common/dft_parameters)
C *In variables                                                        *
INCLUDE(common/dft_basis)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_order_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_memory_info)
      REAL matrix_scr(*)
      logical expand_sw,Vform_sw
      logical oe1c_int_sw
      logical oe2c_ove_sw,oe3c_ove_sw
      logical te2c_rep_sw,te3c_rep_sw
      logical mat_mult_sw
C *Scratch space and pointers                                          *
      integer memory_int(*)
      REAL memory_fp(*)
      integer nprm_pt,angm_pt,hyb_pt,ctre_pt,pstr_pt,alp_pt,cc_pt
      integer expons_pt,d_pt,e_pt,f_pt
      integer sigma_pt
      integer bigP_pt,ss_pt,smallp_pt
      integer s_pt
C *Out variables                                                       *
      REAL matrix_out(*)
C *Local variables                                                     *
      integer size_required
      integer max_angm,max_pr,prim_sum,lsum,lprod
      integer ncentres,a_tag,b_tag
      integer memory_alloc,memory_rel
C *Functions                                                           *
      integer incr_memory
C *End declarations                                                    *
C **********************************************************************
      if(debug_sw) then
        write(6,*)
        write(6,*) 'Entered into intPack'
        write(6,*) 'Integrals to be calculated'
        write(6,*) 'One centre Gaussian integration  :',oe1c_int_sw
        write(6,*) 'Two centre overlap integrals     :',oe2c_ove_sw
        write(6,*) 'Three centre overlap integrals   :',oe3c_ove_sw
        write(6,*) 'Two centre repulsion integrals   :',te2c_rep_sw
        write(6,*) 'Three centre repulsion integrals :',te3c_rep_sw
      endif
      memory_alloc=0
      memory_rel=0

      if(expand_sw) then
C **********************************************************************
C *Expanded basis set
C *
C * Pointers for the expanded basis set
C *
C * Pointer             Length				Memory type
C * nprm_pt		num_bset*natoms*nshells		integer
C * angm_pt		num_bset*natoms*nshells		integer
C * hyb_pt		num_bset*natoms*nshells		integer
C * ctre_pt		num_bset*natoms*nshells		integer
C * pstr_pt		num_bset*natoms*nshells		integer
C * alp_pt		num_bset*natoms*tot_prm		double
C * cc_pt		num_bset*natoms*tot_prm		double
C *
C *Expand the basis set into the integer scratch arrays
C *
        nprm_pt = incr_memory(num_bset*maxi_shlA,'i')
        angm_pt = incr_memory(num_bset*maxi_shlA,'i')
        hyb_pt  = incr_memory(num_bset*maxi_shlA,'i')
        ctre_pt = incr_memory(num_bset*maxi_shlA,'i')
        pstr_pt = incr_memory(num_bset*maxi_shlA,'i')
        alp_pt  = incr_memory(num_bset*maxi_primA,'d')
        cc_pt   = incr_memory(num_bset*maxi_primA,'d')

        if(debug_sw) then
          write(6,*) 'Memory requirements for expanding basis sets'
c         write(6,*) ' Start of integer memory ',start_ifreemem
          write(6,*) ' Integer memory required (per block)',
     &                 maxi_shlA
c         write(6,*) ' Start of floating point memory ',start_dfreemem
          write(6,*) ' Floating point memory required (per block)',
     &                 maxi_primA
          write(6,*) 'Pointer		Value'
          write(6,*) 'nprm_pt		',nprm_pt
          write(6,*) 'angm_pt		',angm_pt
          write(6,*) 'hyb_pt		',hyb_pt
          write(6,*) 'ctre_pt		',ctre_pt
          write(6,*) 'pstr_pt		',pstr_pt
          write(6,*) 'alp_pt		',alp_pt
          write(6,*) 'cc_pt		',cc_pt
        endif
c       call expand_toshells(memory_int(nprm_pt),
c    &                       memory_int(angm_pt),
c    &                       memory_int(hyb_pt),
c    &                       memory_int(ctre_pt),
c    &                       memory_int(pstr_pt),
c    &                       memory_fp(alp_pt),
c    &                       memory_fp(cc_pt))
C *
C **********************************************************************
      endif
      if(oe1c_int_sw) then
C **********************************************************************
C *One electron one centre gaussian integral
C *
C * pointers to the arrays used in the primitive integral routines
C *
C * Pointer     Length                                  Memory type
C * expons	a_nprim					double
C * s		a_nprim					double
C * d		a_nprim*la*la				double
C * e		a_nprim*la*la				double
C * f		a_nprim*la*la				double
C * sigma	a_nprim					double
C *
        a_tag=2
        ncentres=1
        call scratch_size(a_tag,max_pr,max_angm)
        lsum=max_angm
        lprod=max_angm
        prim_sum=max_pr
        size_required=prim_sum*lsum*lprod

        s_pt      = incr_memory(prim_sum,'d')
        expons_pt = incr_memory(prim_sum,'d')
        d_pt      = incr_memory(size_required,'d')
        e_pt      = incr_memory(size_required,'d')
        f_pt      = incr_memory(size_required,'d')
        sigma_pt  = incr_memory(prim_sum,'d')
        memory_alloc=memory_alloc+size_required*3+
     &               prim_sum*3

c       call oe_integral_1c(memory_int(nprm_pt),
c    &                      memory_int(angm_pt),
c    &                      memory_int(hyb_pt),
c    &                      memory_int(ctre_pt),
c    &                      memory_int(pstr_pt),
c    &                      memory_fp(alp_pt),
c    &                      memory_fp(cc_pt),
c    &                      memory_fp(s_pt),
c    &                      memory_fp(expons_pt),
c    &                      memory_fp(d_pt),
c    &                      memory_fp(e_pt),
c    &                      memory_fp(f_pt),
c    &                      memory_fp(sigma_pt),
c    &                      a_tag,matrix_out)
C *
C *Free memory
C *
        call decr_memory(sigma_pt,'d')
        call decr_memory(f_pt,'d')
        call decr_memory(e_pt,'d')
        call decr_memory(d_pt,'d')
        call decr_memory(expons_pt,'d')
        call decr_memory(s_pt,'d')

C *
C **********************************************************************
      endif
      if(oe2c_ove_sw) then
C **********************************************************************
C *One electron two centre overlap integrals
C *
C * pointers to the arrays used in the primitive integral routines
C *
C * Pointer     Length                                  Memory type
C * expons      ab_nprim*2				double
C * ss          ab_nprim				double
C * bigQ?P        ab_nprim*3				double
C * d           ab_nprim*lsum*lprod			double
C * e           ab_nprim*lsum*lprod			double
C * f           ab_nprim*lsum*lprod			double
C * sigma	ab_nprim				double
C * smallP	ab_nprim*3				double
C *
        a_tag=2
        b_tag=2
        ncentres=2
        call scratch_size(a_tag,max_angm,max_pr)
        lsum=max_angm
        lprod=max_angm*max_angm
        prim_sum=max_pr+max_pr
        size_required=prim_sum*lsum*lprod
        write(6,*) 'size_required ',size_required

        ss_pt     = incr_memory(prim_sum,'d')
        expons_pt = incr_memory(2*prim_sum,'d')
        bigP_pt   = incr_memory(3*prim_sum,'d')
        d_pt      = incr_memory(size_required,'d')
        e_pt      = incr_memory(size_required,'d')
        f_pt      = incr_memory(size_required,'d')
        sigma_pt  = incr_memory(prim_sum,'d')
        smallP_pt = incr_memory(3*prim_sum,'d')

c       call oe_overlap_2c(memory_int(nprm_pt),
c    &                     memory_int(angm_pt),
c    &                     memory_int(hyb_pt),
c    &                     memory_int(ctre_pt),
c    &                     memory_int(pstr_pt),
c    &                     memory_fp(alp_pt),
c    &                     memory_fp(cc_pt),
c    &                     memory_fp(ss_pt),
c    &                     memory_fp(expons_pt),
c    &                     memory_fp(bigP_pt),
c    &                     memory_fp(d_pt),
c    &                     memory_fp(e_pt),
c    &                     memory_fp(f_pt),
c    &                     memory_fp(sigma_pt),
c    &                     memory_fp(smallP_pt),
c    &                     a_tag,b_tag,
c    &                     matrix_out)
C *
C *Free memory
C *
        call decr_memory(smallP_pt,'d')
        call decr_memory(sigma_pt,'d')
        call decr_memory(f_pt,'d')
        call decr_memory(e_pt,'d')
        call decr_memory(d_pt,'d')
        call decr_memory(bigP_pt,'d')
        call decr_memory(expons_pt,'d')
        call decr_memory(ss_pt,'d')
C *
C *
C **********************************************************************
      endif
      if(te2c_rep_sw) then
C **********************************************************************
C *Two electron two centre repulsion integrals
C *

         call caserr('obsolete intpack  fn called')
ccc        call te_repulsion_2c(memory_fp,
ccc     &                       matrix_out)

C *
C **********************************************************************
      endif
      if(te3c_rep_sw) then
C **********************************************************************
C *Two electron three centre repulsion integrals

        call caserr('obsolete intpack  fn called')
c        call te_repulsion_3c(memory_fp,
c     &                       Vform_sw,
c     &                       matrix_scr,
c     &                       matrix_out)
C *
C **********************************************************************
      endif
C *
C *Free remaining memory
C *
      if(expand_sw) then

        call decr_memory(cc_pt,'d')
        call decr_memory(alp_pt,'d')
        call decr_memory(pstr_pt,'i')
        call decr_memory(ctre_pt,'i')
        call decr_memory(hyb_pt,'i')
        call decr_memory(angm_pt,'i')
        call decr_memory(nprm_pt,'i')

      endif

      return
      end
      subroutine memreq_o3driver(ao_tag,kf_tag,
     &                           memory_int,memory_fp,Ckfit,
     &                           kma,kmb)   
C **********************************************************************
C *Description:	                                                       *
C *One electron three centre overlap integral driver.                  *
C **********************************************************************
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis)
      integer max_pshell
      parameter(max_pshell=6)
c
      integer ao_tag,kf_tag
      integer memory_int(*)
      REAL memory_fp(*)
      REAL CKfit(*)
      REAL kma,kmb
c
      integer nprm_pt,angm_pt,hyb_pt,ctre_pt,pstr_pt
      integer alp_pt,cc_pt
      integer size_prm,size_coef
      integer sss_pt,expons_pt,bigQ_pt,d_pt,e_pt,f_pt
      integer sigma_pt,smallq_pt,sabc_pt
c
      integer ao_basfn
      integer nu,mu,cf
c
      integer incr_memory
C
C Expand the basis set into scratch arrays
      nprm_pt = incr_memory(num_bset*maxi_shlA,'i')
      angm_pt = incr_memory(num_bset*maxi_shlA,'i')
      hyb_pt  = incr_memory(num_bset*maxi_shlA,'i')
      ctre_pt = incr_memory(num_bset*maxi_shlA,'i')
      pstr_pt = incr_memory(num_bset*maxi_shlA,'i')
      alp_pt  = incr_memory(num_bset*maxi_primA,'d')
      cc_pt   = incr_memory(num_bset*maxi_primA*max_ang,'d')
c     call expand_toshells(memory_int(nprm_pt),
c    &                     memory_int(angm_pt),
c    &                     memory_int(hyb_pt),
c    &                     memory_int(ctre_pt),
c    &                     memory_int(pstr_pt),
c    &                     memory_fp(alp_pt),
c    &                     memory_fp(cc_pt))
      write(6,*) 'num_bset:',num_bset
      write(6,*) 'maxi_shlA:',maxi_shlA
      write(6,*) 'maxi_primA:',maxi_primA
      write(6,*) 'Checking guards after expansion'
      size_prm  = max_pshell**3
      sss_pt    = incr_memory(size_prm,'d')
      expons_pt = incr_memory(3*size_prm,'d')
      bigQ_pt   = incr_memory(3*size_prm,'d')
      size_coef = size_prm*(max_pshell*3)*max_ang**3
      d_pt      = incr_memory(size_coef,'d')
      e_pt      = incr_memory(size_coef,'d')
      f_pt      = incr_memory(size_coef,'d')
      sigma_pt  = incr_memory(size_prm,'d')
      smallQ_pt = incr_memory(3*size_prm,'d')
      Sabc_pt   = incr_memory(size_prm,'d')
      ao_basfn=totbfn(ao_tag)
c     do nu=1,totshl(ao_tag)
c       do mu=1,nu
c         do cf=1,totshl(kf_tag)
c           call oe_overlap_3c(ao_tag,kf_tag,
c    &                         nu,mu,cf,
c    &                         memory_int(nprm_pt),
c    &                         memory_int(angm_pt),
c    &                         memory_int(hyb_pt),
c    &                         memory_int(ctre_pt),
c    &                         memory_int(pstr_pt),
c    &                         memory_fp(alp_pt),
c    &                         memory_fp(cc_pt),
c    &                         memory_fp(sss_pt),
c    &                         memory_fp(expons_pt),
c    &                         memory_fp(bigQ_pt),
c    &                         memory_fp(d_pt),
c    &                         memory_fp(e_pt),
c    &                         memory_fp(f_pt),
c    &                         memory_fp(sigma_pt),
c    &                         memory_fp(smallQ_pt),
c    &                         memory_fp(Sabc_pt))
c         write(6,*) 'Sabc:',memory_fp(Sabc_pt)
c         enddo 
c         call dgemv('N',ao_basfn,ao_basfn,1.0d0,CKfit,ao_basfn,
c    &               memory_fp(Sabc_pt),1,1.0d0,kma,1)
c       enddo
c     enddo
C *
C *Free memory
C *
      call decr_memory(Sabc_pt,'d')
      call decr_memory(smallQ_pt,'d')
      call decr_memory(sigma_pt,'d')
      call decr_memory(f_pt,'d')
      call decr_memory(e_pt,'d')
      call decr_memory(d_pt,'d')
      call decr_memory(bigQ_pt,'d')
      call decr_memory(expons_pt,'d')
      call decr_memory(sss_pt,'d')

      call decr_memory(cc_pt,'d')
      call decr_memory(alp_pt,'d')
      call decr_memory(pstr_pt,'i')
      call decr_memory(ctre_pt,'i')
      call decr_memory(hyb_pt,'i')
      call decr_memory(angm_pt,'i')
      call decr_memory(nprm_pt,'i')
 
      return
      end
c---- the routines that do the real work -------------------------------
      subroutine intPack_drv(memory_int,
     &                       memory_fp,
     &                       matrix_scr,
     &                       expand_sw,Vform_sw,
     &                       oe1c_int_sw,
     &                       oe2c_ove_sw,oe3c_ove_sw,
     &                       te2c_rep_sw,te3c_rep_sw,
     &                       mat_mult_sw,
     &                       matrix_out)
C **********************************************************************
C *Description                                                         *
C *Top level driver for intPack.
C **********************************************************************
      implicit none
C **********************************************************************
C *Declarations                                                        *
C *                                                                    *
C *Parameters                                                          *
INCLUDE(common/dft_parameters)
C *In variables                                                        *
INCLUDE(common/dft_basis)
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_order_info)
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_memory_info)
      REAL matrix_scr(*)
      logical expand_sw,Vform_sw
      logical oe1c_int_sw
      logical oe2c_ove_sw,oe3c_ove_sw
      logical te2c_rep_sw,te3c_rep_sw
      logical mat_mult_sw
C *Scratch space and pointers                                                 *
      integer memory_int(*)
      REAL memory_fp(*)
      integer nprm_pt,angm_pt,hyb_pt,ctre_pt,pstr_pt,alp_pt,cc_pt
      integer expons_pt,d_pt,e_pt,f_pt
      integer sigma_pt
      integer bigP_pt,ss_pt,smallp_pt
      integer s_pt
C *Out variables                                                              *
      REAL matrix_out(*)
C *Local variables                                                            *
      integer size_required
      integer max_angm,max_pr,prim_sum,lsum,lprod
      integer ncentres,a_tag,b_tag
      integer memory_alloc,memory_rel
C *Functions                                                                  *
      integer allocate_memory
C *End declarations                                                           *
C *****************************************************************************
      if(debug_sw) then
        write(6,*)
        write(6,*) 'Entered into intPack'
        write(6,*) 'Integrals to be calculated'
        write(6,*) 'One centre Gaussian integration  :',oe1c_int_sw
        write(6,*) 'Two centre overlap integrals     :',oe2c_ove_sw
        write(6,*) 'Three centre overlap integrals   :',oe3c_ove_sw
        write(6,*) 'Two centre repulsion integrals   :',te2c_rep_sw
        write(6,*) 'Three centre repulsion integrals :',te3c_rep_sw
      endif
      memory_alloc=0
      memory_rel=0

      if(expand_sw) then
C *****************************************************************************
C *Expanded basis set
C *
C * Pointers for the expanded basis set
C *
C * Pointer             Length				Memory type
C * nprm_pt		num_bset*natoms*nshells		integer
C * angm_pt		num_bset*natoms*nshells		integer
C * hyb_pt		num_bset*natoms*nshells		integer
C * ctre_pt		num_bset*natoms*nshells		integer
C * pstr_pt		num_bset*natoms*nshells		integer
C * alp_pt		num_bset*natoms*tot_prm		double
C * cc_pt		num_bset*natoms*tot_prm		double
C *
C *Expand the basis set into the integer scratch arrays
C *
        nprm_pt = allocate_memory(num_bset*maxi_shlA,'i')
        angm_pt = allocate_memory(num_bset*maxi_shlA,'i')
        hyb_pt  = allocate_memory(num_bset*maxi_shlA,'i')
        ctre_pt = allocate_memory(num_bset*maxi_shlA,'i')
        pstr_pt = allocate_memory(num_bset*maxi_shlA,'i')
        alp_pt  = allocate_memory(num_bset*maxi_primA,'d')
        cc_pt   = allocate_memory(num_bset*maxi_primA,'d')

        if(debug_sw) then
          write(6,*) 'Memory requirements for expanding basis sets'
c         write(6,*) ' Start of integer memory ',start_ifreemem
          write(6,*) ' Integer memory required (per block)',
     &                 maxi_shlA
c         write(6,*) ' Start of floating point memory ',start_dfreemem
          write(6,*) ' Floating point memory required (per block)',
     &                 maxi_primA
          write(6,*) 'Pointer		Value'
          write(6,*) 'nprm_pt		',nprm_pt
          write(6,*) 'angm_pt		',angm_pt
          write(6,*) 'hyb_pt		',hyb_pt
          write(6,*) 'ctre_pt		',ctre_pt
          write(6,*) 'pstr_pt		',pstr_pt
          write(6,*) 'alp_pt		',alp_pt
          write(6,*) 'cc_pt		',cc_pt
        endif
        call expand_toshells(memory_int(nprm_pt),
     &                       memory_int(angm_pt),
     &                       memory_int(hyb_pt),
     &                       memory_int(ctre_pt),
     &                       memory_int(pstr_pt),
     &                       memory_fp(alp_pt),
     &                       memory_fp(cc_pt))
C *
C *****************************************************************************
      endif
      if(oe1c_int_sw) then
C *****************************************************************************
C *One electron one centre gaussian integral
C *
C * pointers to the arrays used in the primitive integral routines
C *
C * Pointer     Length                                  Memory type
C * expons	a_nprim					double
C * s		a_nprim					double
C * d		a_nprim*la*la				double
C * e		a_nprim*la*la				double
C * f		a_nprim*la*la				double
C * sigma	a_nprim					double
C *
        a_tag=2
        ncentres=1
        call scratch_size(a_tag,max_pr,max_angm)
        lsum=max_angm
        lprod=max_angm
        prim_sum=max_pr
        size_required=prim_sum*lsum*lprod

        s_pt      = allocate_memory(prim_sum,'d')
        expons_pt = allocate_memory(prim_sum,'d')
        d_pt      = allocate_memory(size_required,'d')
        e_pt      = allocate_memory(size_required,'d')
        f_pt      = allocate_memory(size_required,'d')
        sigma_pt  = allocate_memory(prim_sum,'d')
        memory_alloc=memory_alloc+size_required*3+
     &               prim_sum*3

        call oe_integral_1c(memory_int(nprm_pt),
     &                      memory_int(angm_pt),
     &                      memory_int(hyb_pt),
     &                      memory_int(ctre_pt),
     &                      memory_int(pstr_pt),
     &                      memory_fp(alp_pt),
     &                      memory_fp(cc_pt),
     &                      memory_fp(s_pt),
     &                      memory_fp(expons_pt),
     &                      memory_fp(d_pt),
     &                      memory_fp(e_pt),
     &                      memory_fp(f_pt),
     &                      memory_fp(sigma_pt),
     &                      a_tag,matrix_out)
C *
C *Free memory
C *
        call free_memory(sigma_pt,'d')
        call free_memory(f_pt,'d')
        call free_memory(e_pt,'d')
        call free_memory(d_pt,'d')
        call free_memory(expons_pt,'d')
        call free_memory(s_pt,'d')

C *
C *****************************************************************************
      endif
      if(oe2c_ove_sw) then
C *****************************************************************************
C *One electron two centre overlap integrals
C *
C * pointers to the arrays used in the primitive integral routines
C *
C * Pointer     Length                                  Memory type
C * expons      ab_nprim*2				double
C * ss          ab_nprim				double
C * bigQ?P        ab_nprim*3				double
C * d           ab_nprim*lsum*lprod			double
C * e           ab_nprim*lsum*lprod			double
C * f           ab_nprim*lsum*lprod			double
C * sigma	ab_nprim				double
C * smallP	ab_nprim*3				double
C *
        a_tag=2
        b_tag=2
        ncentres=2
        call scratch_size(a_tag,max_angm,max_pr)
        lsum=max_angm
        lprod=max_angm*max_angm
        prim_sum=max_pr+max_pr
        size_required=prim_sum*lsum*lprod
        write(6,*) 'size_required ',size_required

        ss_pt     = allocate_memory(prim_sum,'d')
        expons_pt = allocate_memory(2*prim_sum,'d')
        bigP_pt   = allocate_memory(3*prim_sum,'d')
        d_pt      = allocate_memory(size_required,'d')
        e_pt      = allocate_memory(size_required,'d')
        f_pt      = allocate_memory(size_required,'d')
        sigma_pt  = allocate_memory(prim_sum,'d')
        smallP_pt = allocate_memory(3*prim_sum,'d')

        call oe_overlap_2c(memory_int(nprm_pt),
     &                     memory_int(angm_pt),
     &                     memory_int(hyb_pt),
     &                     memory_int(ctre_pt),
     &                     memory_int(pstr_pt),
     &                     memory_fp(alp_pt),
     &                     memory_fp(cc_pt),
     &                     memory_fp(ss_pt),
     &                     memory_fp(expons_pt),
     &                     memory_fp(bigP_pt),
     &                     memory_fp(d_pt),
     &                     memory_fp(e_pt),
     &                     memory_fp(f_pt),
     &                     memory_fp(sigma_pt),
     &                     memory_fp(smallP_pt),
     &                     a_tag,b_tag,
     &                     matrix_out)
C *
C *Free memory
C *
        call free_memory(smallP_pt,'d')
        call free_memory(sigma_pt,'d')
        call free_memory(f_pt,'d')
        call free_memory(e_pt,'d')
        call free_memory(d_pt,'d')
        call free_memory(bigP_pt,'d')
        call free_memory(expons_pt,'d')
        call free_memory(ss_pt,'d')
C *
C *
C *****************************************************************************
      endif
      if(te2c_rep_sw) then
C *****************************************************************************
C *Two electron two centre repulsion integrals
C *

         call caserr('obsolete intpack  fn called')
ccc        call te_repulsion_2c(memory_fp,
ccc     &                       matrix_out)

C *
C *****************************************************************************
      endif
      if(te3c_rep_sw) then
C *****************************************************************************
C *Two electron three centre repulsion integrals

        call caserr('obsolete intpack  fn called')
c        call te_repulsion_3c(memory_fp,
c     &                       Vform_sw,
c     &                       matrix_scr,
c     &                       matrix_out)
C *
C *****************************************************************************
      endif
C *
C *Free remaining memory
C *
      if(expand_sw) then

        call free_memory(cc_pt,'d')
        call free_memory(alp_pt,'d')
        call free_memory(pstr_pt,'i')
        call free_memory(ctre_pt,'i')
        call free_memory(hyb_pt,'i')
        call free_memory(angm_pt,'i')
        call free_memory(nprm_pt,'i')

      endif

      return
      end
      subroutine scratch_size(tag,max_pr,max_angm)
C *****************************************************************************
C *Description:								      *
C *Basic routine which calculates the largest primitive contraction and the   *
C *greatest L value. Will always give values which are larger than needed,    *
C *hence will be replaced with a more efficient method of determining the     *
C *required memory. Will be obsolete in F90 version of CCP1-DFT	 	      *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations								      *
C *									      *
C *Parameters								      *
INCLUDE(common/dft_parameters)
C *In variables								      *
INCLUDE(common/dft_basis)
      integer tag
C *Out variables							      *
      integer max_pr,max_angm
C *Local variables							      *
      integer lcen,lshl
C *End declarations							      *
C *****************************************************************************
      max_angm=0
      max_pr=0
      do lcen=1,num_types(tag)
        do lshl=1,num_shl(tag,lcen)
          max_angm=max(max_angm,angmom(tag,lcen,lshl))
          max_pr=max(max_pr,nprim(tag,lcen,lshl))
        enddo
      enddo
      return
      end 
      subroutine order_fill
C *****************************************************************************
C *Description:								      *
C *Form array which determines basis function ordering in matrices            *
C *This routine is package dependent. The following packages are supported    *
C * GAMESS-UK							              *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations								      *
C *									      *
C *In variables								      *
INCLUDE(common/dft_order_info)
C *Out variables							      *

C *Local variables							      *
      REAL sqrt3,sqrt5,sqrt7
C *End declarations							      *
C *****************************************************************************
      sqrt3=dsqrt(3.0d0)
      sqrt5=dsqrt(5.0d0)
      sqrt7=dsqrt(7.0d0)
C *****************************************************************************
C *Gamess - UK basis function ordering					      * 
C *
C *S
      bf_order(1,1)  = 0
      bf_order(1,2)  = 0
      bf_order(1,3)  = 0
C *
C *P
C *
      bf_order(2,1)  = 1
      bf_order(2,2)  = 0
      bf_order(2,3)  = 0
      bf_order(3,1)  = 0
      bf_order(3,2)  = 1
      bf_order(3,3)  = 0
      bf_order(4,1)  = 0
      bf_order(4,2)  = 0
      bf_order(4,3)  = 1
C *
C *D
C *
      bf_order(5,1)  = 2
      bf_order(5,2)  = 0
      bf_order(5,3)  = 0
      bf_order(6,1)  = 0
      bf_order(6,2)  = 2
      bf_order(6,3)  = 0
      bf_order(7,1)  = 0
      bf_order(7,2)  = 0
      bf_order(7,3)  = 2
      bf_order(8,1)  = 1
      bf_order(8,2)  = 1
      bf_order(8,3)  = 0
      bf_order(9,1)  = 1
      bf_order(9,2)  = 0
      bf_order(9,3)  = 1
      bf_order(10,1) = 0
      bf_order(10,2) = 1
      bf_order(10,3) = 1
C *
C *F
C *
      bf_order(11,1)  = 3
      bf_order(11,2)  = 0
      bf_order(11,3)  = 0
      bf_order(12,1)  = 0
      bf_order(12,2)  = 3
      bf_order(12,3)  = 0
      bf_order(13,1)  = 0
      bf_order(13,2)  = 0
      bf_order(13,3)  = 3
      bf_order(14,1)  = 2
      bf_order(14,2)  = 1
      bf_order(14,3)  = 0
      bf_order(15,1)  = 2
      bf_order(15,2)  = 0
      bf_order(15,3)  = 1
      bf_order(16,1)  = 1
      bf_order(16,2)  = 2
      bf_order(16,3)  = 0      
      bf_order(17,1)  = 0
      bf_order(17,2)  = 2
      bf_order(17,3)  = 1
      bf_order(18,1)  = 1
      bf_order(18,2)  = 0
      bf_order(18,3)  = 2
      bf_order(19,1)  = 0
      bf_order(19,2)  = 1
      bf_order(19,3)  = 2
      bf_order(20,1)  = 1
      bf_order(20,2)  = 1
      bf_order(20,3)  = 1
c *
c *G
c *
      bf_order(21,1)  = 4
      bf_order(21,2)  = 0
      bf_order(21,3)  = 0

      bf_order(22,1)  = 0
      bf_order(22,2)  = 4
      bf_order(22,3)  = 0

      bf_order(23,1)  = 0
      bf_order(23,2)  = 0
      bf_order(23,3)  = 4

      bf_order(24,1)  = 3
      bf_order(24,2)  = 1
      bf_order(24,3)  = 0

      bf_order(25,1)  = 3
      bf_order(25,2)  = 0
      bf_order(25,3)  = 1

      bf_order(26,1)  = 1
      bf_order(26,2)  = 3
      bf_order(26,3)  = 0

      bf_order(27,1)  = 0
      bf_order(27,2)  = 3
      bf_order(27,3)  = 1

      bf_order(28,1)  = 1
      bf_order(28,2)  = 0
      bf_order(28,3)  = 3

      bf_order(29,1)  = 0
      bf_order(29,2)  = 1
      bf_order(29,3)  = 3

      bf_order(30,1)  = 2
      bf_order(30,2)  = 2
      bf_order(30,3)  = 0

      bf_order(31,1)  = 2
      bf_order(31,2)  = 0
      bf_order(31,3)  = 2

      bf_order(32,1)  = 0
      bf_order(32,2)  = 2
      bf_order(32,3)  = 2

      bf_order(33,1)  = 2
      bf_order(33,2)  = 1
      bf_order(33,3)  = 1

      bf_order(34,1)  = 1
      bf_order(34,2)  = 2
      bf_order(34,3)  = 1

      bf_order(35,1)  = 1
      bf_order(35,2)  = 1
      bf_order(35,3)  = 2

C *End Gamess - UK basis function ordering				      *
C *****************************************************************************
C *****************************************************************************
C *Basis function indices
C *
C *S
C *
      bf_start(1)  = 1
      bf_end(1)    = 1
      bf_num(1)    = 1
C *
C *P
C *
      bf_start(2)  = 2
      bf_end(2)    = 4
      bf_num(2)    = 3
C *
C *D
C *
      bf_start(3)  = 5
      bf_end(3)    = 10
      bf_num(3)    = 6
C *
C *F
C *
      bf_start(4)  = 11
      bf_end(4)    = 20
      bf_num(4)    = 10
C *
C *G
C *
      bf_start(5)  = 21
      bf_end(5)    = 35
      bf_num(5)    = 15

C *End basis function indices
C *****************************************************************************
C *****************************************************************************
C *Basis function extra normalisation
C *
C *S
C *
      bf_norm(1)  = 1.0d0
C *
C *P
C *
      bf_norm(2)  = 1.0d0
      bf_norm(3)  = 1.0d0
      bf_norm(4)  = 1.0d0
C *
C *D
C *
      bf_norm(5)  = 1.0d0
      bf_norm(6)  = 1.0d0
      bf_norm(7)  = 1.0d0
      bf_norm(8)  = sqrt3
      bf_norm(9)  = sqrt3
      bf_norm(10) = sqrt3
C *
C *F
C *
c
c corrected ps 11/2/97 - but doesn't seem to be used ??
c
      bf_norm(11)  = 1.0d0
      bf_norm(12)  = 1.0d0
      bf_norm(13)  = 1.0d0
      bf_norm(14)  = sqrt5
      bf_norm(15)  = sqrt5
      bf_norm(16)  = sqrt5
      bf_norm(17)  = sqrt5
      bf_norm(18)  = sqrt5
      bf_norm(19)  = sqrt5
      bf_norm(20)  = sqrt5*sqrt3
C *
C *G
C *
c
      bf_norm(21)  = 1.0d0
      bf_norm(22)  = 1.0d0
      bf_norm(23)  = 1.0d0
      bf_norm(24)  = sqrt7
      bf_norm(25)  = sqrt7
      bf_norm(26)  = sqrt7
      bf_norm(27)  = sqrt7
      bf_norm(28)  = sqrt7
      bf_norm(29)  = sqrt7
      bf_norm(30)  = sqrt7*sqrt5/sqrt3
      bf_norm(31)  = sqrt7*sqrt5/sqrt3
      bf_norm(32)  = sqrt7*sqrt5/sqrt3
      bf_norm(33)  = sqrt7*sqrt5
      bf_norm(34)  = sqrt7*sqrt5
      bf_norm(35)  = sqrt7*sqrt5

C *
C *End basis function extra normalisation
C *****************************************************************************
      return
      end
      subroutine ang_norm(alpha,la,an,al,am,norml)
      implicit none
INCLUDE(common/dft_physical_constants)
      REAL alpha
      integer la,an,al,am
c     REAL norml(5)
      REAL norml
c     integer lang,c
c     REAL pi2,n14
      REAL p1,p2
      REAL ffact(5)
      ffact(1)=1.0d0
      ffact(2)=1.0d0
c     ffact(3)=0.5773502
c     ffact(4)=0.2581988
c     ffact(5)=0.0975900
      ffact(3)=3.0d0
      ffact(4)=15.00d0
      ffact(5)=105.00d0
c     pi2=2.0d0/pi
c     n14=1.0d0/4.0d0
c     p1=(pi2*alpha)**(n14)
c     c=0
c     do lang=0,la 
c       c=c+1
c       p2=p1*((4.0d0*alpha)**(dble(lang/2.0d0)))
c       norml(c)=p2*ffact(c)
c     enddo

      p1=(pi/(2.0d0*alpha))**(3.0d0/2.0d0)
      p2=(4.0d0*alpha)**(an+al+am)
      norml=(ffact(an+1)*ffact(al+1)*ffact(am+1)*(p1/p2))**(-0.5d0)
      return
      end
      subroutine o3driver(ao_tag,kf_tag,
     &                    memory_int,memory_fp,Ckfit,
     &                    kma,kmb)   
C *****************************************************************************
C *Description:								      *
C *One electron three centre overlap integral driver.			      *
C *****************************************************************************
      implicit none
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis)
      integer max_pshell
      parameter(max_pshell=6)
c
      integer ao_tag,kf_tag
      integer memory_int(*)
      REAL memory_fp(*)
      REAL CKfit(*)
      REAL kma,kmb
c
      integer nprm_pt,angm_pt,hyb_pt,ctre_pt,pstr_pt
      integer alp_pt,cc_pt
      integer size_prm,size_coef
      integer sss_pt,expons_pt,bigQ_pt,d_pt,e_pt,f_pt
      integer sigma_pt,smallq_pt,sabc_pt
c
      integer ao_basfn
      integer nu,mu,cf
c
      integer allocate_memory
C
C Expand the basis set into scratch arrays
      nprm_pt = allocate_memory(num_bset*maxi_shlA,'i')
      angm_pt = allocate_memory(num_bset*maxi_shlA,'i')
      hyb_pt  = allocate_memory(num_bset*maxi_shlA,'i')
      ctre_pt = allocate_memory(num_bset*maxi_shlA,'i')
      pstr_pt = allocate_memory(num_bset*maxi_shlA,'i')
      alp_pt  = allocate_memory(num_bset*maxi_primA,'d')
      cc_pt   = allocate_memory(num_bset*maxi_primA*max_ang,'d')
      call expand_toshells(memory_int(nprm_pt),
     &                     memory_int(angm_pt),
     &                     memory_int(hyb_pt),
     &                     memory_int(ctre_pt),
     &                     memory_int(pstr_pt),
     &                     memory_fp(alp_pt),
     &                     memory_fp(cc_pt))
      write(6,*) 'num_bset:',num_bset
      write(6,*) 'maxi_shlA:',maxi_shlA
      write(6,*) 'maxi_primA:',maxi_primA
      write(6,*) 'Checking guards after expansion'
      size_prm  = max_pshell**3
      sss_pt    = allocate_memory(size_prm,'d')
      expons_pt = allocate_memory(3*size_prm,'d')
      bigQ_pt   = allocate_memory(3*size_prm,'d')
      size_coef = size_prm*(max_pshell*3)*max_ang**3
      d_pt      = allocate_memory(size_coef,'d')
      e_pt      = allocate_memory(size_coef,'d')
      f_pt      = allocate_memory(size_coef,'d')
      sigma_pt  = allocate_memory(size_prm,'d')
      smallQ_pt = allocate_memory(3*size_prm,'d')
      Sabc_pt   = allocate_memory(size_prm,'d')
      ao_basfn=totbfn(ao_tag)
      do nu=1,totshl(ao_tag)
        do mu=1,nu
          do cf=1,totshl(kf_tag)
            call oe_overlap_3c(ao_tag,kf_tag,
     &                         nu,mu,cf,
     &                         memory_int(nprm_pt),
     &                         memory_int(angm_pt),
     &                         memory_int(hyb_pt),
     &                         memory_int(ctre_pt),
     &                         memory_int(pstr_pt),
     &                         memory_fp(alp_pt),
     &                         memory_fp(cc_pt),
     &                         memory_fp(sss_pt),
     &                         memory_fp(expons_pt),
     &                         memory_fp(bigQ_pt),
     &                         memory_fp(d_pt),
     &                         memory_fp(e_pt),
     &                         memory_fp(f_pt),
     &                         memory_fp(sigma_pt),
     &                         memory_fp(smallQ_pt),
     &                         memory_fp(Sabc_pt))
          write(6,*) 'Sabc:',memory_fp(Sabc_pt)
          enddo 
c         call dgemv('N',ao_basfn,ao_basfn,1.0d0,CKfit,ao_basfn,
c    &               memory_fp(Sabc_pt),1,1.0d0,kma,1)
        enddo
      enddo
C *
C *Free memory
C *
      call free_memory(Sabc_pt,'d')
      call free_memory(smallQ_pt,'d')
      call free_memory(sigma_pt,'d')
      call free_memory(f_pt,'d')
      call free_memory(e_pt,'d')
      call free_memory(d_pt,'d')
      call free_memory(bigQ_pt,'d')
      call free_memory(expons_pt,'d')
      call free_memory(sss_pt,'d')

      call free_memory(cc_pt,'d')
      call free_memory(alp_pt,'d')
      call free_memory(pstr_pt,'i')
      call free_memory(ctre_pt,'i')
      call free_memory(hyb_pt,'i')
      call free_memory(angm_pt,'i')
      call free_memory(nprm_pt,'i')
 
      return
      end
      subroutine oe_integral_1c(nprim_e,angmom_e,hybrid_e,
     &                          centre_e,pstart_e,
     &                          alpha_e,cont_coeff_e,
     &                          s,expons,d,e,f,sigma,a_tag,
     &                          matrix)
C *****************************************************************************
C *Description								      *
C *One centre gaussian integral routine					      *
C *									      *
C * Variables								      *
C * --------- 								      *
C * a_nprim 	Total number of primitives for a given shell                  *
C * expons	Array containing exponents for each primitive in shell        * 
C * s		Stores contraction coefficient for given primitive in shell   *
C * d		Coefficients described by McM & D., This is the x set.        *
C * e		Coefficients described by McM & D., This is the y set.        *
C * f		Coefficients described by McM & D., This is the z set.        *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations                                                               *
C *									      *
C *Parameters								      *
INCLUDE(common/dft_parameters)
C *In variables								      *
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_basis)
      integer a_tag

      integer nprim_e(num_bset,*),angmom_e(num_bset,*)
      integer hybrid_e(num_bset,*),centre_e(num_bset,*)
      integer pstart_e(num_bset,*)
      REAL alpha_e(num_bset,*),cont_coeff_e(num_bset,*)

      REAL s(*),expons(*),sigma(*)
      REAL d(*),e(*),f(*)
C *Out variables                                                              *
      REAL matrix(*)
C *Local variables                                                            *
      integer ashell
      integer a_nprim,lsum,la,ah
      integer matpos
C *End declarations                                                           *
C *****************************************************************************
C *****************************************************************************
C *Loop over shells                                                           
c     integer i

c      write(6,*)'OVERLAP'

c      do i=1,10
c        write(6,*) nprim_e(a_tag,i),angmom_e(a_tag,i),
c     &        hybrid_e(a_tag,i),
c     &        alpha_e(a_tag,i),cont_coeff_e(a_tag,i)
c      enddo

      matpos=0
      do ashell=1,totshl(a_tag)
c         write(6,*)'shell',ashell
C *
C *Form the necessary prefactors
C *
         call shl_inf(nprim_e,angmom_e,hybrid_e,centre_e,pstart_e,
     &               alpha_e,cont_coeff_e,
     &               num_bset,a_tag,ashell,la,ah,a_nprim,s,expons)
C *
C *Generate the D coefficients
C *
        call D_form_1c(la,a_nprim,s,expons,d,e,f,sigma)
C *
C *Now form the contracted integral
C *
c        write(6,*)'la,ah',la,ah
c        write(6,*)'d e f',d(1),e(1),f(1)

        call mat_form_1c(a_nprim,expons,s,lsum,la,ah,d,e,f,matrix,
     &                   matpos)

      enddo
C *End loop over shells							      *
C *****************************************************************************
      return
      end
      subroutine shl_inf(nprim_e,angmom_e,hybrid_e,centre_e,pstart_e,
     &                   alpha_e,cont_coeff_e,
     &                   num_bset,a_tag,ashell,la,ah,a_nprim,s,expons)
C *****************************************************************************
C *Description								      *
C *Form shell preliminary factors					      *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations                                                               *
C *									      *
C *Parameters								      *
INCLUDE(common/dft_parameters)
C *In variables								      *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_physical_constants)
INCLUDE(common/dft_numbers)
      integer num_bset
      integer nprim_e(num_bset,*)
      integer angmom_e(num_bset,*)
      integer hybrid_e(num_bset,*)
      integer centre_e(num_bset,*)
      integer pstart_e(num_bset,*)
      REAL alpha_e(num_bset,*),cont_coeff_e(num_bset,*)
      integer a_tag
      integer ashell
C *Out variables                                                              *
      integer la,ah,a_nprim
      REAL expons(*)
      REAL s(*)
C *Local variables							      *
      integer a_start
      integer icount,lprma,lap
      REAL alpha,exp_ab
C *Functions                                                                  *
c     REAL dij2
C *End declarations                                                           *
C *****************************************************************************
C *****************************************************************************
C *Gather shell information						      *

      a_start = pstart_e(a_tag,ashell)
      a_nprim = nprim_e(a_tag,ashell)
      la      = angmom_e(a_tag,ashell)-1
      ah      = hybrid_e(a_tag,ashell)-1

c      write(6,*)'shlinf',pstart_e(a_tag,ashell),nprim_e(a_tag,ashell),la,lh
C *
C *For each shell form s integral and store exponents     
C *
      icount=0 
      do lprma=a_start,(a_start+a_nprim)-1
        icount=icount+1
        expons(icount)=alpha_e(a_tag,lprma)
        s(icount)=cont_coeff_e(a_tag,lprma)
c        write(6,*) 'EXPS:',lprma,expons(icount),s(icount),
c     &       cont_coeff_e(a_tag,lprma)
      enddo
C *
C *Loop over all primitives in shell and form factors
C * 
      icount=0
      do lap=1,a_nprim
        alpha          = expons(lap)
        exp_ab         = sqrt(pi/alpha)
        icount         = icount + 1
        expons(icount) = alpha
        s(icount)      = exp_ab * s(lap)
      enddo
      a_nprim=icount
C *
C *Finish gathering shell information					      *
C *****************************************************************************
      return
      end


      subroutine D_form_1c(la,a_nprim,s,expons,d,e,f,sigma)
C *****************************************************************************
C *Description:								      *
C *Form the d, e and f coefficents described by McMurchie and Davidson.       *
C *L. E. McMurchie and E. R. Davidson J. Comput. Phys. 26 (1978) 218-231      *
C *The d coefficients also contain exponent and contraction information.      *
C *									      *
C *The correspondance between variables and notation in the paper is:         *
C *									      *
C * xL = L = n + l + m where x corresponds to the shell.                      *
C * d = dNnn                                                                  *
C * e = eLll								      *
C * f = fMmm								      *
C * K is a sum of N, L and M                                                  *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations 							      *
C *									      *
C *In variables 							      *
INCLUDE(common/dft_physical_constants)
      integer a_nprim,la
      REAL expons(*)
      REAL s(*)
      REAL sigma(*)
C *Out variables                                                              *
      REAL d(a_nprim,0:la,0:la)
      REAL e(a_nprim,0:la,0:la)
      REAL f(a_nprim,0:la,0:la)
C *Local variables							      *
      integer lap,loa,lK
      REAL alpha,st
      REAL array_size
C *End declarations 							      *
C *****************************************************************************
      array_size=a_nprim*la*la
C * 
C *Set each x coefficient to equal ss
C *
      do lap=1,a_nprim
        alpha = expons(lap)
        st = sqrt(pi/alpha)
        d(lap,0,0)=s(lap)
        e(lap,0,0)=st
        f(lap,0,0)=st
        sigma(lap)=0.5d0/alpha
      enddo
C *****************************************************************************
C *For shell a form the coefficients dn+1nn, el+1ll, fm+1mm	 	      *
      if(la.gt.0) then
        
        do lap=1,a_nprim
          d(lap,0,1)=0.0d0
          e(lap,0,1)=0.0d0
          f(lap,0,1)=0.0d0
          d(lap,1,1)=sigma(lap)*d(lap,0,0)
          e(lap,1,1)=sigma(lap)*e(lap,0,0)
          f(lap,1,1)=sigma(lap)*f(lap,0,0)
        enddo
       do loa=2,la
C *
C *K = 0
C *
          do lap=1,a_nprim
            d(lap,0,loa)  = d(lap,1,loa-1)
            e(lap,0,loa)  = e(lap,1,loa-1)
            f(lap,0,loa)  = f(lap,1,loa-1)
          enddo
C *
C *K ne 0, la
C *
          do lK=1,loa-2
            do lap=1,a_nprim
              d(lap,lK,loa) = sigma(lap)   *   d(lap,lK-1,loa-1)
     &                      + (lK+1)       *   d(lap,lK+1,loa-1)
              e(lap,lK,loa) = sigma(lap)   *   e(lap,lK-1,loa-1)
     &                      + (lK+1)       *   e(lap,lK+1,loa-1)
              f(lap,lK,loa) = sigma(lap)   *   f(lap,lK-1,loa-1)
     &                      + (lK+1)       *   f(lap,lK+1,loa-1)
            enddo
          enddo
C *
C *K = loa
C *
          lK=loa-1
          do lap=1,a_nprim
            d(lap,lK,loa)   = sigma(lap)   *   d(lap,lK-1,loa-1)
            e(lap,lK,loa)   = sigma(lap)   *   e(lap,lK-1,loa-1)
            f(lap,lK,loa)   = sigma(lap)   *   f(lap,lK-1,loa-1)
            d(lap,lK+1,loa) = sigma(lap)   *   d(lap,lK  ,loa-1)
            e(lap,lK+1,loa) = sigma(lap)   *   e(lap,lK  ,loa-1)
            f(lap,lK+1,loa) = sigma(lap)   *   f(lap,lK  ,loa-1)
          enddo
C *
C *End K
C *
        enddo
      endif
C *									      *
C *****************************************************************************
      return
      end
      subroutine mat_form_1c(a_nprim,ex,s,lsum,la,ah,d,e,f,matrix,
     &                       matpos)
C *****************************************************************************
C *Description:								      *
C *Form 1 centre integral matrix					      *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations								      *
C *									      *
C *In variables								      *
INCLUDE(common/dft_order_info)
      integer la,ah
      integer a_nprim,lsum
      REAL ex(*),s(*)
      REAL d(a_nprim,0:la,0:la)
      REAL e(a_nprim,0:la,0:la)
      REAL f(a_nprim,0:la,0:la)
C *Out variables							      *
      REAL matrix(*)
C *Local variables							      *
      integer aijk,anang,ha
      integer lap,matpos
      integer an,al,am
      REAL enorm
C *End declarations 							      *
C *****************************************************************************
C *****************************************************************************
C *Loop over hybrid shells then over all combinations of angular components   *
C *for shell 
      do ha=ah,la
c       anang=(((1+ha)*(2+ha))/2)+1
        anang=ha+1
        do aijk=bf_start(anang),bf_end(anang)
c           write(6,*) 'AIJK:',aijk
          an=bf_order(aijk,1)
          al=bf_order(aijk,2)
          am=bf_order(aijk,3)
          matpos=matpos+1
          matrix(matpos)=0.0d0
          enorm=1.0d0
c          write(6,*)'matpos',matpos
          do lap=1,a_nprim
c            enorm=bf_norm(an+1)*bf_norm(al+1)*bf_norm(am+1)
             matrix(matpos)=matrix(matpos)+ (d(lap,0,an)
     &                                    *  e(lap,0,al)
     &                                    *  f(lap,0,am))*enorm
c             write(6,*) 'Es:',an,d(lap,0,an),al,e(lap,0,al),am,
c     &            f(lap,0,am),d(lap,0,1)
          enddo
        enddo
      enddo
C *End loop								      *
C *****************************************************************************
      return
      end
      subroutine oe_overlap_2c(nprim_e,angmom_e,hybrid_e,
     &                         centre_e,pstart_e,
     &                         alpha_e,cont_coeff_e,
     &                         ss,expons,bigP,
     &                         d,e,f,sigma,smallP,
     &                         a_tag,b_tag,
     &                         matrix)
C *****************************************************************************
C *Description								      *
C *Two centre overlap routine. Uses McMurchie and Davidson method.            *
C *+++Needs tidying up and placed within modulaI formulism                    *
C *									      *
C * Variables								      *
C * --------- 								      *
C * ab_nprim 	Total number of primitives for a given shell doublet          *
C * expons	Array containing two exponents per 2 centre primitive         *
C * ss		Basic primitive integral for given 2 centre primitive         *
C * bigP	Position of 2 centre primitive				      *
C * d		Coefficients described by McM & D., This is the x set.        *
C * e		Coefficients described by McM & D., This is the y set.        *
C * f		Coefficients described by McM & D., This is the z set.        *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations                                                               *
C *									      *
C *Parameters								      *
INCLUDE(common/dft_parameters)
C *In variables								      *
INCLUDE(common/dft_module_comm)
INCLUDE(common/dft_basis)
      integer a_tag,b_tag
      integer nprim_e(num_bset,*),angmom_e(num_bset,*)
      integer hybrid_e(num_bset,*),centre_e(num_bset,*)
      integer pstart_e(num_bset,*)
      REAL alpha_e(*),cont_coeff_e(*)
      REAL ss(*),expons(*),bigP(*)
      REAL d(*),e(*),f(*)
      REAL sigma(*),smallP(*)
C *Out variables                                                              *
      REAL matrix(*)
C *Local variables                                                            *
      integer istart,ashell,bshell
      integer ab_nprim,lsum,la,lb,ah,bh
      integer matpos
      REAL bigA,bigB
C *End declarations                                                           *
C *****************************************************************************
C *****************************************************************************
C *Loop over shells                                                           
      istart=1
      matpos=1
C      if(a_tag.eq.b_tag) triangle_sw = .true.
cc      write(6,*) 'shl',totshl(a_tag),totshl(b_tag)
      do ashell=1,totshl(a_tag)
        if(triangle_sw) istart=ashell
        do bshell=istart,totshl(b_tag)
C *
C *Form the necessary prefactors
C *
          call shlpair_inf(nprim_e,angmom_e,hybrid_e,
     &                     centre_e,pstart_e,
     &                     alpha_e,cont_coeff_e,
     &                     num_bset,
     &                     a_tag,b_tag,
     &                     ashell,bshell,
     &                     bigA,bigB,
     &                     la,lb,ah,bh,ab_nprim,
     &                     ss,expons,bigP)
          lsum=la+lb
C *
C *Generate the D coefficients
C *
          call D_form_2c(la,lb,bigA,bigB,
     &                   ab_nprim,
     &                   ss,expons,bigP,
     &                   d,e,f,lsum,
     &                   sigma,smallP)
C *
C *Now form the contracted 2 centre overlap integral
C *
          call Smat_form_2c(ab_nprim,lsum,la,lb,d,e,f,ah,bh,
     &                      matrix,matpos)
        enddo
      enddo
C *End loop over shells							      *
C *****************************************************************************
      return
      end
      subroutine shlpair_inf(nprim_e,angmom_e,hybrid_e,
     &                       centre_e,pstart_e,
     &                       alpha_e,cont_coeff_e,
     &                       num_bset,
     &                       a_tag,b_tag,
     &                       ashell,bshell,
     &                       bigA,bigB,
     &                       la,lb,ah,bh,ab_nprim,
     &                       ss,expons,bigP)
C *****************************************************************************
C *Description								      *
C *Form shell pair preliminary factors                                        *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations                                                               *
C *									      *
C *Parameters								      *
INCLUDE(common/dft_parameters)
C *In variables								      *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_physical_constants)
INCLUDE(common/dft_numbers)
      integer num_bset,a_tag,b_tag
      integer nprim_e(num_bset,*)
      integer angmom_e(num_bset,*)
      integer hybrid_e(num_bset,*)
      integer centre_e(num_bset,*)
      integer pstart_e(num_bset,*)
      REAL alpha_e(num_bset,*),cont_coeff_e(num_bset,*)
      integer ashell,bshell
C *Out variables                                                              *
      integer la,lb,ah,bh,ab_nprim
      REAL bigA(3),bigB(3)
      REAL expons(2,*)
      REAL ss(*),bigP(3,*)
C *Local variables							      *
      integer a_start,a_nprim,atomA
      integer b_start,b_nprim,atomB
      integer icount,lprma,lprmb,lab
      REAL alpha,beta,sigmaP,isigmaP,exp_ab
C *Functions                                                                  *
c     REAL dij2
C *End declarations                                                           *
C *****************************************************************************
C *****************************************************************************
C *Gather shell information						      *
cc      write(6,*) 'Atag - ',a_tag,ashell
cc      write(6,*) 'Btag - ',b_tag,bshell
      a_start = pstart_e(a_tag,ashell)
      a_nprim = nprim_e(a_tag,ashell)
      la      = angmom_e(a_tag,ashell)
      atomA   = centre_e(a_tag,ashell)
      ah      = hybrid_e(a_tag,ashell)
 
      b_start = pstart_e(b_tag,bshell)
      b_nprim = nprim_e(b_tag,bshell)
      lb      = angmom_e(b_tag,bshell)
      atomB   = centre_e(b_tag,bshell)
      bh      = hybrid_e(b_tag,bshell)

      bigA(1)=atom_c(atomA,1)
      bigA(2)=atom_c(atomA,2)
      bigA(3)=atom_c(atomA,3)
      bigB(1)=atom_c(atomB,1)
      bigB(2)=atom_c(atomB,2)
      bigB(3)=atom_c(atomB,3)
C *
C *For each shell doublet form ss integral and store exponents     
C *
      icount=0 
      ab_nprim=a_nprim*b_nprim
cc      write(6,*) 'ab_nprim ',ab_nprim
cc      write(6,*) 'A',a_start,a_nprim
cc      write(6,*) 'B',b_start,b_nprim
      do lprma=a_start,(a_start+a_nprim)-1
        do lprmb=b_start,(b_start+b_nprim)-1
          icount=icount+1
          expons(1,icount)=alpha_e(a_tag,lprma)
          expons(2,icount)=alpha_e(b_tag,lprmb)
          ss(icount)=cont_coeff_e(a_tag,lprma)*
     &               cont_coeff_e(b_tag,lprmb)
cc          write(6,*) 'ex',expons(1,icount),expons(2,icount),ss(icount)
        enddo
      enddo
C *
C *Loop over all primitives in shell doublet and form factors
C *Remove all primitives whose exponent is less than global accuracy
C * 
      icount=0
      do lab=1,ab_nprim
        alpha   = expons(1,lab)
        beta    = expons(2,lab)
        sigmaP  = alpha+beta
        isigmaP = 1.0d0/sigmaP
        bigP(1,lab)=isigmaP*(alpha*bigA(1)+beta*bigB(1))
        bigP(2,lab)=isigmaP*(alpha*bigA(2)+beta*bigB(2))
        bigP(3,lab)=isigmaP*(alpha*bigA(3)+beta*bigB(3))
        exp_ab=(pi*isigmaP)**n15
        if(exp_ab.gt.global_accuracy) then
          icount=icount+1
          expons(1,icount)=alpha
          expons(2,icount)=beta
          bigP(1,icount)=bigP(1,lab)
          bigP(2,icount)=bigP(2,lab)
          bigP(3,icount)=bigP(3,lab)
          ss(icount)=exp_ab*ss(lab)
        endif
      enddo
      ab_nprim=icount
C *
C *Finish gathering shell information					      *
C *****************************************************************************
      return
      end
      subroutine D_form_2c(la,lb,
     &                     bigA,bigB,
     &                     ab_nprim,
     &                     ss,expons,bigP,
     &                     d,e,f,lsum,
     &                     sigma,smallP)
C *****************************************************************************
C *Description:								      *
C *Form the d, e and f coefficents described by McMurchie and Davidson.       *
C *L. E. McMurchie and E. R. Davidson J. Comput. Phys. 26 (1978) 218-231      *
C *The d coefficients also contain exponent and contraction information.      *
C *									      *
C *The correspondance between variables and notation in the paper is:         *
C *									      *
C * xL = L = n + l + m where x corresponds to the shell.                      *
C * d = dNnn                                                                  *
C * e = eLll								      *
C * f = fMmm								      *
C * K is a sum of N, L and M                                                  *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations 							      *
C *									      *
C *In variables 							      *
      integer ab_nprim,lsum,la,lb
      REAL expons(2,*)
      REAL ss(*),bigP(3,*)
      REAL sigma(*),smallP(3,*)
      REAL bigA(3),bigB(3)
C *Out variables                                                              *
      REAL d(ab_nprim,0:lsum,0:la,0:lb)
      REAL e(ab_nprim,0:lsum,0:la,0:lb)
      REAL f(ab_nprim,0:lsum,0:la,0:lb)
C *Local variables							      *
      integer lab,loa,lob,lK
      REAL alpha,beta
      REAL array_size
C *End declarations 							      *
C *****************************************************************************
cc      write(6,*) 'In dform'
      array_size=ab_nprim*lsum*la*lb
cc      write(6,*) 'array_size ',array_size
      do lab=1,ab_nprim
        do lK=1,lsum
          do loa=1,la
            do lob=1,lb
              d(lab,lK,loa,lob)=0.0d0
              e(lab,lK,loa,lob)=0.0d0
              f(lab,lK,loa,lob)=0.0d0
            enddo
          enddo
        enddo
      enddo
c     call aclear_dp(d,array_size,0.0d0)
c     call aclear_dp(e,array_size,0.0d0)
c     call aclear_dp(f,array_size,0.0d0)
C * 
C *Set each x coefficient to equal ss
C *
cc      write(6,*) 'starting small loop'
      do lab=1,ab_nprim
        d(lab,0,0,0)=ss(lab)
        e(lab,0,0,0)=1.0d0
        f(lab,0,0,0)=1.0d0
        alpha = expons(1,lab)
        beta  = expons(2,lab)
        sigma(lab)=0.5d0/alpha+beta
        smallP(1,lab)=bigP(1,lab)-bigA(3)
        smallP(2,lab)=bigP(2,lab)-bigA(3) 
        smallP(3,lab)=bigP(3,lab)-bigA(3)
      enddo
cc      write(6,*) 'out of small loop'
C *****************************************************************************
C *For shell a form the coefficients dn+1nn, el+1ll, fm+1mm	 	      *
      do loa=1,la
C *
C *K = 0
C *
        do lab=1,ab_nprim
          d(lab,0,loa,0)  = smallP(1,lab) *   d(lab,0,  loa-1,0)
     &                    +                   d(lab,1,  loa-1,0)
          e(lab,0,loa,0)  = smallP(2,lab) *   e(lab,0,  loa-1,0)
     &                    +                   e(lab,1,  loa-1,0)
          f(lab,0,loa,0)  = smallP(3,lab) *   f(lab,0,  loa-1,0)
     &                    +                   f(lab,1,  loa-1,0)
        enddo
C *
C *K ne 0, la
C *
        do lK=1,loa-1
          do lab=1,ab_nprim
            d(lab,lK,loa,0) = sigma(lab)   *   d(lab,lK-1,loa-1,0)
     &                      + smallP(1,lab)*   d(lab,lK,  loa-1,0)
     &                      + (lK+1)       *   d(lab,lK+1,loa-1,0)
            e(lab,lK,loa,0) = sigma(lab)   *   e(lab,lK-1,loa-1,0)
     &                      + smallP(2,lab)*   e(lab,lK,  loa-1,0)
     &                      + (lK+1)       *   e(lab,lK+1,loa-1,0)
            f(lab,lK,loa,0) = sigma(lab)   *   f(lab,lK-1,loa-1,0)
     &                      + smallP(3,lab)*   f(lab,lK,  loa-1,0)
     &                      + (lK+1)       *   f(lab,lK+1,loa-1,0)
          enddo
        enddo
C *
C *K = loa
C *
        do lab=1,ab_nprim
          d(lab,lK,loa,0) = sigma(lab)     *   d(lab,lK-1,loa-1,0)
     &                    + smallP(1,lab)  *   d(lab,lK,  loa-1,0)
          e(lab,lK,loa,0) = sigma(lab)     *   e(lab,lK-1,loa-1,0)
     &                    + smallP(2,lab)  *   e(lab,lK,  loa-1,0)
          f(lab,lK,loa,0) = sigma(lab)     *   f(lab,lK-1,loa-1,0)
     &                    + smallP(3,lab)  *   f(lab,lK,  loa-1,0)
        enddo
C *
C *End K
C *
      enddo
C *									      *
C *****************************************************************************
      do lab=1,ab_nprim
        smallP(1,lab)=bigP(1,lab)-bigB(1)
        smallP(2,lab)=bigP(2,lab)-bigB(2)
        smallP(3,lab)=bigP(3,lab)-bigB(3)
      enddo
C *****************************************************************************
C *For each shell pair a,b form the coefficients			      *
      do lob=1,lb
        do loa=0,la
C *
C *K = 0
C *
          do lab=1,ab_nprim
            d(lab,0,loa,lob) = smallP(1,lab) *   d(lab,0,  loa,lob-1)
     &                       +                   d(lab,1,  loa,lob-1)
            e(lab,0,loa,lob) = smallP(2,lab) *   e(lab,0,  loa,lob-1)
     &                       +                   e(lab,1,  loa,lob-1)
            f(lab,0,loa,lob) = smallP(3,lab) *   f(lab,0,  loa,lob-1)
     &                       +                   f(lab,1,  loa,lob-1)
          enddo
C *
C *K ne 0, la + lb
C *
          do lK=1,la+lb-1
            do lab=1,ab_nprim
              d(lab,lK,loa,lob)= sigma(lab)   *   d(lab,lK-1,loa,lob-1)
     &                         + smallP(1,lab)*   d(lab,lK,  loa,lob-1)
     &                         + (lK+1)       *   d(lab,lK+1,loa,lob-1)
              e(lab,lK,loa,lob)= sigma(lab)   *   e(lab,lK-1,loa,lob-1)
     &                         + smallP(2,lab)*   e(lab,lK,  loa,lob-1)
     &                         + (lK+1)       *   e(lab,lK+1,loa,lob-1)
              f(lab,lK,loa,lob)= sigma(lab)   *   f(lab,lK-1,loa,lob-1)
     &                         + smallP(3,lab)*   f(lab,lK,  loa,lob-1)
     &                         + (lK+1)       *   f(lab,lK+1,loa,lob-1)
            enddo
          enddo
C *
C *K = la + lb
C *
          do lab=1,ab_nprim
            d(lab,lK,loa,lob) = sigma(lab)    *   d(lab,lK-1,loa,lob-1)
     &                        + smallP(1,lab) *   d(lab,lK,  loa,lob-1)
            e(lab,lK,loa,lob) = sigma(lab)    *   e(lab,lK-1,loa,lob-1)
     &                        + smallP(2,lab) *   e(lab,lK,  loa,lob-1)
            f(lab,lK,loa,lob) = sigma(lab)    *   f(lab,lK-1,loa,lob-1)
     &                        + smallP(3,lab) *   f(lab,lK,  loa,lob-1)
          enddo
C *
C *End K
C *
        enddo 
      enddo
cc      write(6,*) 'Out of dform'
C *									      *
C *****************************************************************************
      return
      end
      subroutine Smat_form_2c(ab_nprim,la,lb,lsum,d,e,f,
     &                        ah,bh,
     &                        matrix,matpos)
C *****************************************************************************
C *Description:								      *
C *Form 2 centre overlap matrix elements for basis functions in shell doublet *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations								      *
C *									      *
C *In variables								      *
INCLUDE(common/dft_order_info)
      integer la,lb,ah,bh
      integer ab_nprim,lsum
      REAL d(ab_nprim,0:lsum,0:la,0:lb)
      REAL e(ab_nprim,0:lsum,0:la,0:lb)
      REAL f(ab_nprim,0:lsum,0:la,0:lb)
C *Out variables							      *
      REAL matrix(*)
C *Local variables							      *
      integer aijk,bijk,anang,bnang,ha,hb
      integer lab,matpos
      integer an,al,am,bn,bl,bm
C *End declarations 							      *
C *****************************************************************************
C *****************************************************************************
C *Loop over hybrid shells then over all combinations of angular components   *
C *for shell doublet         						      *
c      write(6,*) ab_nprim,ah,la,bh,lb
      do ha=ah,la
        do hb=bh,lb
          anang=((1+ha)*(2+ha))/2
          bnang=((1+bh)*(2+bh))/2
          do aijk=1,anang
c            write(6,*) 'aijk ',aijk
            an=bf_order(aijk,1)
            al=bf_order(aijk,2)
            am=bf_order(aijk,3)
c            write(6,*) 'aijk ',aijk,an,al,am
            do bijk=1,bnang
              bn=bf_order(bijk,1)
              bl=bf_order(bijk,2)
              bm=bf_order(bijk,3)
c              write(6,*) 'bijk ',bijk,bn,bl,bm
              matpos=matpos+1
              matrix(matpos)=0.0d0
              do lab=1,ab_nprim
                matrix(matpos)=matrix(matpos)+ d(lab,0,an,bn)
     &                                       + e(lab,0,al,bl)
     &                                       + f(lab,0,am,bm)
              enddo
c              write(6,*) 'Matrix ',matrix(matpos),matpos
            enddo
          enddo
        enddo
      enddo
C *End loop								      *
C *****************************************************************************
      return
      end
C *****************************************************************************
C *Description								      *
C *Collection of routines which form 3 centre overlap integrals. 	      *
C *The module uses the McMurchie and Davidson method to form the integrals.   *
C *See L. E. McMurchie and E. R. Davidson J. Comput. Phys. 26 (1978) 218-231  *
C * and V. R. Saunders in Methods in Computational Physics, NATO ASI 1982     * 
C *									      *
C *Notes								      *
C *This is son of oe_overlap_2c						      *
C *									      *
C *Variables and what they mean						      *
C *abc_nprim 	Total number of primitives within a given shell triplet       *
C *expons	Array containing three exponents per 3 centre primitive       *
C *sss		Basic primitive integral for given 3 centre primitive         *
C *bigQ		Position of 3 centre primitive				      *
C *d		Coefficients described by McM & D., This is the x set.        *
C *e		Coefficients described by McM & D., This is the y set.        *
C *f		Coefficients described by McM & D., This is the z set.        *
C *****************************************************************************
      subroutine oe_overlap_3c(ao_tag,kf_tag,
     &                         nu,mu,cf,
     &                         nprim_e,angmom_e,hybrid_e,
     &                         centre_e,pstart_e,
     &                         alpha_e,cont_coeff_e,
     &                         sss,expons,bigQ,
     &                         d,e,f,sigma,smallQ,
     &                         Sabc)
C *****************************************************************************
C *Description:								      *
C *Low level driver							      *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations								      *
C *                                                                           *
C *Parameters                                                                 *
INCLUDE(common/dft_parameters)
C *In variables                                                               *
INCLUDE(common/dft_basis)
INCLUDE(common/dft_module_comm)
      integer ao_tag,kf_tag
      integer nprim_e(*),angmom_e(*),hybrid_e(*),centre_e(*)
      integer pstart_e(*)
c
      REAL alpha_e(*),cont_coeff_e(*)
      REAL sss(*),expons(*),bigQ(*)
      REAL d(*),e(*),f(*)
      REAL sigma(*),smallQ(*)
c
      integer nu,mu,cf
C *Out variables                                                              *
      REAL Sabc(*)
C *Local variables                                                            *
      integer abc_nprim,lsum
      integer matpos
      integer la,lb,lc,ah,bh,ch
      REAL bigA(3),bigB(3),bigC(3)
C *End declarations							      *
C *****************************************************************************
      call shltriple_inf(nprim_e,angmom_e,hybrid_e,
     &                   centre_e,pstart_e,
     &                   alpha_e,cont_coeff_e,
     &                   num_bset,
     &                   ao_tag,kf_tag,
     &                   nu,mu,cf,
     &                   bigA,bigB,bigC,
     &                   la,lb,lc,ah,bh,ch,
     &                   abc_nprim,
     &                   sss,expons,bigQ)
      lsum=la+lb+lc
C *
C *Generate the D coefficients
C *
      call D_form_3c(la,lb,lc,
     &               bigA,bigB,bigC,
     &               abc_nprim,
     &               sss,expons,bigQ,
     &               d,e,f,lsum,
     &               sigma,smallQ)
C *
C *Now form the contracted 3 centre overlap integral
C *
      call Smat_form_3c(abc_nprim,lsum,la,lb,lc,d,e,f,
     &                  ah,bh,ch,Sabc,matpos)
      return
      end
      subroutine shltriple_inf(nprim_e,angmom_e,hybrid_e,
     &                         centre_e,pstart_e,
     &                         alpha_e,cont_coeff_e,
     &                         num_bset,
     &                         ao_tag,kf_tag,
     &                         nu,mu,cf,
     &                         bigA,bigB,bigC,
     &                         la,lb,lc,ah,bh,ch,
     &                         abc_nprim,
     &                         sss,expons,bigQ)
C *****************************************************************************
C *Description								      *
C *Form shell triples preliminary factors                                     *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations                                                               *
C *									      *
C *Parameters								      *
INCLUDE(common/dft_parameters)
C *In variables								      *
INCLUDE(common/dft_mol_info)
INCLUDE(common/dft_physical_constants)
INCLUDE(common/dft_numbers)
      integer num_bset
      integer nprim_e(num_bset,*)
      integer angmom_e(num_bset,*)
      integer hybrid_e(num_bset,*)
      integer centre_e(num_bset,*)
      integer pstart_e(num_bset,*)
      REAL alpha_e(num_bset,*),cont_coeff_e(num_bset,*)
      integer ao_tag,kf_tag
      integer nu,mu,cf
C *Out variables                                                              *
      integer la,lb,lc,ah,bh,ch,abc_nprim
      REAL bigA(3),bigB(3),bigC(3)
      REAL expons(3,*)
      REAL sss(*)
      REAL bigQ(3,*)
C *Local variables							      *
      integer a_start,a_nprim,atomA
      integer b_start,b_nprim,atomB
      integer c_start,c_nprim,atomC
      integer icount,lprma,lprmb,lprmc,labc
      REAL alpha,beta,gamma,zeta
      REAL sigmaP,isigmaP,exp_ab
      REAL isigmaP_c,fact,exp_abc
      REAL rpc2,top_of(3),rab2
      REAL bigP(3)
C *Functions                                                                  *
      REAL dij2
C *End declarations                                                           *
C *****************************************************************************
C *****************************************************************************
C *Gather shell information						      *
      a_start = pstart_e(ao_tag,nu)
      a_nprim = nprim_e(ao_tag,nu)
      la      = angmom_e(ao_tag,nu)-1
      atomA   = centre_e(ao_tag,nu)
      ah      = hybrid_e(ao_tag,nu)
      
      b_start = pstart_e(ao_tag,mu)
      b_nprim = nprim_e(ao_tag,mu)
      lb      = angmom_e(ao_tag,mu)-1
      atomB   = centre_e(ao_tag,mu)
      bh      = hybrid_e(ao_tag,mu)

      c_start = pstart_e(kf_tag,cf)
      c_nprim = nprim_e(kf_tag,cf)
      lc      = angmom_e(kf_tag,cf)-1
      atomC   = centre_e(kf_tag,cf)
      ch      = hybrid_e(kf_tag,cf)
 
      bigA(1)=atom_c(atomA,1)
      bigA(2)=atom_c(atomA,2)
      bigA(3)=atom_c(atomA,3)
      bigB(1)=atom_c(atomB,1)
      bigB(2)=atom_c(atomB,2)
      bigB(3)=atom_c(atomB,3)
      bigC(1)=atom_c(atomC,1)
      bigC(2)=atom_c(atomC,2)
      bigC(3)=atom_c(atomC,3)
C *
C *For given shell triplet load contraction coefficients into sss array 
C *and store exponents 
C *
      icount=0 
      abc_nprim=a_nprim*b_nprim*c_nprim 
      do lprmc=c_start,(c_start+c_nprim)-1
      do lprma=a_start,(a_start+a_nprim)-1
        do lprmb=b_start,(b_start+b_nprim)-1
c         do lprmc=c_start,(c_start+c_nprim)-1
            icount=icount+1
c           expons(1,icount)=alpha_e(ao_tag,lprma)
c           expons(2,icount)=alpha_e(ao_tag,lprmb)  
c           expons(3,icount)=alpha_e(kf_tag,lprmc)
            expons(2,icount)=alpha_e(ao_tag,lprma)
            expons(3,icount)=alpha_e(ao_tag,lprmb)
            expons(1,icount)=alpha_e(kf_tag,lprmc)
c           sss(icount)=cont_coeff_e(ao_tag,lprma)*
c    &                  cont_coeff_e(ao_tag,lprmb)*
c    &                  1.0d0
c    &                  cont_coeff_e(kf_tag,lprmc)
            sss(icount)=1.0d0
          enddo
        enddo
      enddo
      write(6,*) 'ABC_nprim:',abc_nprim,icount
C *
C *Loop over all primitives in shell triplet and form factors
C *Remove all primitives whose exponent is less than global accuracy
C * 
      icount=0
      RAB2=dij2(bigA,bigB)
      do labc=1,abc_nprim
        alpha     = expons(1,labc)
        beta      = expons(2,labc)
        gamma     = expons(3,labc)
        sigmaP    = alpha+beta
        isigmaP   = 1.0d0/sigmaP
        zeta      = alpha*beta*isigmaP
        exp_ab    = exp(-zeta*RAB2)
        top_of(1) = (alpha*bigA(1)+beta*bigB(1))
        top_of(2) = (alpha*bigA(2)+beta*bigB(2))
        top_of(3) = (alpha*bigA(3)+beta*bigB(3))
        bigP(1)   = isigmaP*top_of(1)
        bigP(2)   = isigmaP*top_of(2)
        bigP(3)   = isigmaP*top_of(3)
        isigmaP_c = 1.0d0/(sigmaP+gamma)
        RPC2      = dij2(bigP,bigC)
        exp_abc   = exp_ab*exp(-sigmaP*gamma*isigmaP_c*RPC2)
        fact      = sqrt(pi*isigmaP_c)
        exp_abc   = exp_abc*fact*fact*fact
        write(6,*) 'EXP_ABC:',exp_abc,labc,alpha,beta,gamma
c       if(exp_abc.gt.global_accuracy) then
          icount=icount+1
          expons(1,icount)=alpha
          expons(2,icount)=beta
          expons(3,icount)=gamma
          bigQ(1,icount)=isigmaP_c*(top_of(1)+gamma*bigC(1))
          bigQ(2,icount)=isigmaP_c*(top_of(2)+gamma*bigC(2))
          bigQ(3,icount)=isigmaP_c*(top_of(3)+gamma*bigC(3))
          sss(icount)=exp_abc*sss(labc)
c       endif
      enddo
      abc_nprim=icount
C *
C *Finish gathering shell information					      *
C *****************************************************************************
      return
      end
      subroutine D_form_3c(la,lb,lc,
     &                     bigA,bigB,bigC,
     &                     abc_nprim,
     &                     sss,expons,bigQ,
     &                     d,e,f,lsum,
     &                     sigma,smallQ)
C *****************************************************************************
C *Description:								      *
C *Form the D coefficents described by McMurchie and Davidson.                *
C *L. E. McMurchie and E. R. Davidson J. Comput. Phys. 26 (1978) 218-231      *
C *This is a three centre version. Each coefficient also contains exponent    *
C *and contraction information.                                               *
C *									      *
C *The correspondance between variables and notation in the paper is:         *
C *									      *
C * xL = L = n + l + m where x corresponds to the shell.                      *
C * d = dNnnn                                                                 *
C * e = eLlll								      *
C * f = fMmmm								      *
C * K is a sum of N, L and M                                                  *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations 							      *
C *									      *
C *In variables 							      *
      integer abc_nprim,lsum,la,lb,lc
      REAL expons(3,*)
      REAL sss(*)
      REAL bigQ(3,*)
      REAL bigA(3),bigB(3),bigC(3)
C *Scratch space and pointers						      *
      REAL sigma(*)
      REAL smallQ(3,*)
C *Out variables                                                              *
      REAL d(abc_nprim,0:lsum,0:la,0:lb,0:lc)
      REAL e(abc_nprim,0:lsum,0:la,0:lb,0:lc)
      REAL f(abc_nprim,0:lsum,0:la,0:lb,0:lc)
C *Local variables							      *
      integer labc,loa,lob,loc,lK
      REAL alpha,beta,gamma
      REAL array_size
C *End declarations 							      *
C *****************************************************************************
      array_size=abc_nprim*lsum*la*lb*lc
      do lK=0,lsum
        do loa=0,la
          do lob=0,lb
            do loc=0,lc
              d(1,lK,loa,lob,loc)=0.0d0
              e(1,lK,loa,lob,loc)=0.0d0
              f(1,lK,loa,lob,loc)=0.0d0
            enddo
          enddo
        enddo
      enddo
C * 
C *Set each x coefficient to equal sss
C *
      do labc=1,abc_nprim
        d(labc,0,0,0,0)=sss(labc)
        e(labc,0,0,0,0)=1.0d0
        f(labc,0,0,0,0)=1.0d0
        alpha = expons(1,labc)
        beta  = expons(2,labc)
        gamma = expons(3,labc)
        sigma(labc)=0.5d0/alpha+beta+gamma
        sigma(labc+abc_nprim)=sigma(labc)
        sigma(labc+2*abc_nprim)=sigma(labc)
        smallQ(1,labc)=bigQ(1,labc)-bigA(3)
        smallQ(2,labc)=bigQ(2,labc)-bigA(3) 
        smallQ(3,labc)=bigQ(3,labc)-bigA(3)
      enddo
C *****************************************************************************
C *For shell a form the coefficients dn+1nn, el+1ll, fm+1mm	 	      *
      do loa=1,la
C *
C *K = 0
C *
        do labc=1,abc_nprim
      d(labc,0,loa,0,0) = smallQ(1,labc)*   d(labc,0,loa-1,0,0)
     &                  +                   d(labc,1,loa-1,0,0)
      e(labc,0,loa,0,0) = smallQ(2,labc)*   e(labc,0,loa-1,0,0)
     &                  +                   e(labc,1,loa-1,0,0)
      f(labc,0,loa,0,0) = smallQ(3,labc)*   f(labc,0,loa-1,0,0)
     &                  +                   f(labc,1,loa-1,0,0)
        enddo
C *
C *K ne 0, la
C *
        do lK=1,loa-1
          do labc=1,abc_nprim
      d(labc,lK,loa,0,0) = sigma(labc)   *   d(labc,lK-1,loa-1,0,0)
     &                   + smallQ(1,labc)*   d(labc,lK,  loa-1,0,0)
     &                   + (lK+1)        *   d(labc,lK+1,loa-1,0,0)
      e(labc,lK,loa,0,0) = sigma(labc)   *   e(labc,lK-1,loa-1,0,0)
     &                   + smallQ(2,labc)*   e(labc,lK,  loa-1,0,0)
     &                   + (lK+1)        *   e(labc,lK+1,loa-1,0,0)
      f(labc,lK,loa,0,0) = sigma(labc)   *   f(labc,lK-1,loa-1,0,0)
     &                   + smallQ(3,labc)*   f(labc,lK,  loa-1,0,0)
     &                   + (lK+1)        *   f(labc,lK+1,loa-1,0,0)
          enddo
        enddo
C *
C *K = a
C *
        do labc=1,abc_nprim
      d(labc,lK,loa,0,0) = sigma(labc)   *   d(labc,lK-1,loa-1,0,0)
     &                   + smallQ(1,labc)*   d(labc,lK,  loa-1,0,0)
      e(labc,lK,loa,0,0) = sigma(labc)   *   e(labc,lK-1,loa-1,0,0)
     &                   + smallQ(2,labc)*   e(labc,lK,  loa-1,0,0)
      f(labc,lK,loa,0,0) = sigma(labc)   *   f(labc,lK-1,loa-1,0,0)
     &                   + smallQ(3,labc)*   f(labc,lK,  loa-1,0,0)
        enddo
C *
C *End K
C *
      enddo
C *									      *
C *****************************************************************************
      do labc=1,abc_nprim
        smallQ(1,labc)=bigQ(1,labc)-bigB(1)
        smallQ(2,labc)=bigQ(2,labc)-bigB(2)
        smallQ(3,labc)=bigQ(3,labc)-bigB(3)
      enddo
C *****************************************************************************
C *For each shell pair a,b form the coefficients			      *
      do lob=1,lb
        do loa=0,la
C *
C *K = 0
C *
          do labc=1,abc_nprim
      d(labc,0,loa,lob,0) = smallQ(1,labc)*   d(labc,0,  loa,lob-1,0)
     &                    +                   d(labc,1,  loa,lob-1,0)
      e(labc,0,loa,lob,0) = smallQ(2,labc)*   e(labc,0,  loa,lob-1,0)
     &                    +                   e(labc,1,  loa,lob-1,0)
      f(labc,0,loa,lob,0) = smallQ(3,labc)*   f(labc,0,  loa,lob-1,0)
     &                    +                   f(labc,1,  loa,lob-1,0)
          enddo
C *
C *K ne 0, la + lb 
C *
          do lK=1,loa+lob-1
            do labc=1,abc_nprim
      d(labc,lK,loa,lob,0)= sigma(labc)   *   d(labc,lK-1,loa,lob-1,0)
     &                    + smallQ(1,labc)*   d(labc,lK,  loa,lob-1,0)
     &                    + (lK+1)        *   d(labc,lK+1,loa,lob-1,0)
      e(labc,lK,loa,lob,0)= sigma(labc)   *   e(labc,lK-1,loa,lob-1,0)
     &                    + smallQ(2,labc)*   e(labc,lK,  loa,lob-1,0)
     &                    + (lK+1)        *   e(labc,lK+1,loa,lob-1,0)
      f(labc,lK,loa,lob,0)= sigma(labc)   *   f(labc,lK-1,loa,lob-1,0)
     &                    + smallQ(3,labc)*   f(labc,lK,  loa,lob-1,0)
     &                    + (lK+1)        *   f(labc,lK+1,loa,lob-1,0)
            enddo
          enddo
C *
C *K = la + lb 
C *
          do labc=1,abc_nprim
      d(labc,lK,loa,lob,0) = sigma(labc)   *  d(labc,lK-1,loa,lob-1,0)
     &                     + smallQ(1,labc)*  d(labc,lK,  loa,lob-1,0)
      e(labc,lK,loa,lob,0) = sigma(labc)   *  e(labc,lK-1,loa,lob-1,0)
     &                     + smallQ(2,labc)*  e(labc,lK,  loa,lob-1,0)
      f(labc,lK,loa,lob,0) = sigma(labc)   *  f(labc,lK-1,loa,lob-1,0)
     &                     + smallQ(3,labc)*  f(labc,lK,  loa,lob-1,0)
      
          enddo
C *
C *End K
C *
        enddo 
      enddo
C *									      *
C *****************************************************************************
      do labc=1,abc_nprim
        smallQ(1,labc)=bigQ(1,labc)-bigC(1)
        smallQ(2,labc)=bigQ(2,labc)-bigC(2)
        smallQ(3,labc)=bigQ(3,labc)-bigC(3)
      enddo
C *****************************************************************************
C *For each shell triple a,b,c form the coefficients			      *
      do loc=1,lc
        do lob=0,lb  
          do loa=0,la 
C *
C *K = 0
C *
            do labc=1,abc_nprim
      d(labc,0,loa,lob,0) = smallQ(1,labc)*   d(labc,0,  loa,lob,loc-1)
     &                    +                   d(labc,1,  loa,lob,loc-1)
      e(labc,0,loa,lob,0) = smallQ(2,labc)*   e(labc,0,  loa,lob,loc-1)
     &                    +                   e(labc,1,  loa,lob,loc-1)
      f(labc,0,loa,lob,0) = smallQ(3,labc)*   f(labc,0,  loa,lob,loc-1)
     &                    +                   f(labc,1,  loa,lob,loc-1)
            enddo
C *
C *K ne 0, loa + lob + loc
C *
            do lK=1,loa+lob+loc-1
              do labc=1,abc_nprim
      d(labc,lK,loa,lob,0) = sigma(labc)*   d(labc,lK-1,loa,lob,loc-1)
     &                  + smallQ(1,labc)*   d(labc,lK,  loa,lob,loc-1)
     &                  + (lK+1)        *   d(labc,lK+1,loa,lob,loc-1)
      e(labc,lK,loa,lob,0) = sigma(labc)*   e(labc,lK-1,loa,lob,loc-1)
     &                  + smallQ(2,labc)*   e(labc,lK,  loa,lob,loc-1)
     &                  + (lK+1)        *   e(labc,lK+1,loa,lob,loc-1)
      f(labc,lK,loa,lob,0) = sigma(labc)*   f(labc,lK-1,loa,lob,loc-1)
     &                  + smallQ(3,labc)*   f(labc,lK,  loa,lob,loc-1)
     &                  + (lK+1)        *   f(labc,lK+1,loa,lob,loc-1)

              enddo
            enddo
C *
C * K = la + lb + lc
C *
            do labc=1,abc_nprim
      d(labc,lK,loa,lob,0) = sigma(labc)*   d(labc,lK-1,loa,lob,loc-1)
     &                  + smallQ(1,labc)*   d(labc,lK,  loa,lob,loc-1)
      e(labc,lK,loa,lob,0) = sigma(labc)*   e(labc,lK-1,loa,lob,loc-1)
     &                  + smallQ(2,labc)*   e(labc,lK,  loa,lob,loc-1)
      f(labc,lK,loa,lob,0) = sigma(labc)*   f(labc,lK-1,loa,lob,loc-1)
     &                  + smallQ(3,labc)*   f(labc,lK,  loa,lob,loc-1)
            enddo
C *
C *End of K
C *
          enddo 
        enddo
      enddo
C *									      *
C *****************************************************************************
      return
      end
      subroutine Smat_form_3c(abc_nprim,la,lb,lc,lsum,d,e,f,
     &                        ah,bh,ch,
     &                        Sabc,matpos)
C *****************************************************************************
C *Description:								      *
C *Form 3 centre overlap matrix elements for basis functions in shell triple  *
C *****************************************************************************
        implicit none
C *****************************************************************************
C *Declarations								      *
C *									      *
C *In variables								      *
INCLUDE(common/dft_order_info)
      integer abc_nprim,lsum,la,lb,lc
      REAL d(abc_nprim,0:lsum,0:la,0:lb,0:lc)
      REAL e(abc_nprim,0:lsum,0:la,0:lb,0:lc)
      REAL f(abc_nprim,0:lsum,0:la,0:lb,0:lc)
      integer ah,bh,ch
C *Out variables							      *
      REAL Sabc(*)
C *Local variables							      *
      integer matpos
      integer aijk,bijk,cijk,anang,bnang,cnang
      integer an,al,am
      integer bn,bl,bm
      integer cn,cl,cm
      integer labc
C *End declarations 							      *
C *****************************************************************************
C *****************************************************************************
C *Loop over hybrid shells then over all combinations of angular components   *
C *for shell triplet                                                          *
      matpos=0
      anang=((1+la)*(2+la))/2
      bnang=((1+lb)*(2+lb))/2
      cnang=((1+lc)*(2+lc))/2
      write(6,*)'loops:',anang,bnang,cnang,abc_nprim
      do aijk=1,anang
        an=bf_order(aijk,1)
        al=bf_order(aijk,2)
        am=bf_order(aijk,3)
        do bijk=1,bnang
          bn=bf_order(bijk,1)
          bl=bf_order(bijk,2)
          bm=bf_order(bijk,3)
          do cijk=1,cnang
            cn=bf_order(cijk,1)
            cl=bf_order(cijk,2)
            cm=bf_order(cijk,3)
            matpos=matpos+1
            Sabc(matpos)=0.0d0
            do labc=1,abc_nprim
              Sabc(matpos)=Sabc(matpos)
     &                    +d(labc,0,an,bn,cn)
     &                    *e(labc,0,al,bl,cl)
     &                    *f(labc,0,am,bm,cm)
      write(6,*) 'COEFS:',d(labc,0,an,bn,cn),e(labc,0,al,bl,cl),f(labc,0
     +,am,bm,cm)
            enddo
          enddo
        enddo
      enddo
      write(6,*) 'matpos:',matpos
      return
      end

      subroutine te_norm_2c(memory_fp,matrix_out)
C *****************************************************************************
C *Description:                                                               *
C *Two centre 2 electron repulsion integral routine                           *
C *****************************************************************************
      implicit none
C *****************************************************************************
C *Declarations
C *Parameters
INCLUDE(common/dft_parameters)
INCLUDE(common/dft_basis)
INCLUDE(common/dft_memory_info)
C *In variables
      REAL memory_fp(*)
C *Out variables
      REAL matrix_out(*)
C *Local variables
      integer gout_pt,cd_tag
C *Functions
      integer allocate_memory
C *End declarations                                                           *
C *****************************************************************************
C *****************************************************************************
C *Pointers for jkintx
C *
C * Pointer			Length			Type
C * gout_pt			50625			double
      cd_tag=2
      gout_pt=allocate_memory(50625,'d')

      call te2c_rep_norm(cd_tag,memory_fp(gout_pt),
     &                   matrix_out)

      call free_memory(gout_pt,'d')
      return
      end
      subroutine ver_dft_intpack(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/dft/intpack.m,v $
     +     "/
      data revision /
     +     "$Revision: 5774 $"
     +      /
      data date /
     +     "$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $"
     +     /
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end
