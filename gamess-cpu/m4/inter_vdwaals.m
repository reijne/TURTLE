c=======================================================================
c
c     begin van der Waals correction stuff
c
      subroutine vdwaals_input
c
_IFN(f77)
      use vdwaals_corr
      implicit none
c
      character*4  ytest         ! keywords read in
      character*4  yelm          ! element name
      character*4  ymodel        ! c6 pair coefficient model name
      character*10 dft_xc        ! the DFT functional
      integer      ielm          ! atomic number
      integer      ierr          ! error code
      REAL         val           ! either c6 coefficient or radius
      logical :: odone = .false. ! are we done yet
c
      integer  isubst
      external isubst
c
c     First extract the functional if any is available
c
      call vdw_set_c6pair("geometric",ierr)
      call vdw_get_functional(dft_xc)
      call vdw_set_defaults
      call vdw_switch_on
      call vdw_set_functional(dft_xc)
      do while (.not.odone)
        call input
        call inpa4(ytest)
        select case (ytest)
        case ('on')                      ! on
          call vdw_switch_on 
        case ('off')                     ! off
          call vdw_switch_off 
        case ('scal')                    ! scale <scale factor>
          call inpf(val)
          if (val.lt.0.0d0) then
            call caserr("negative scale factor is not allowed")
          endif
          call vdw_set_scale(val)
        case ('alph')                    ! alpha <exponent>
          call inpf(val)
          if (val.lt.0.0d0) then
            call caserr("negative exponent is not allowed")
          endif
          call vdw_set_alpha(val)
        case ('rad','radi')              ! radius <element> <radius>
          call inpa4(yelm)
          ielm = isubst(yelm)
          call inpf(val)
          call vdw_set_radius(ielm,val)
        case ('c6')                      ! c6 <element> <c6 coefficient>
          call inpa4(yelm)
          ielm = isubst(yelm)
          call inpf(val)
          call vdw_set_c6coef(ielm,val)
        case ('c6mo')                    ! c6model <model>
          call inpa4(ymodel)
          ierr = 0
          select case (ymodel)
            case ("aver")
              call vdw_set_c6pair("average",ierr)
            case ("geom")
              call vdw_set_c6pair("geometric",ierr)
            case default
              call caserr("invalid c6 pair coefficient model")
          end select
          if (ierr.gt.0) then
            call caserr("invalid c6 pair coefficient model")
          endif
        case ('end')
          odone = .true.                 ! end
        case default
          call caserr('unknown directive in van der Waals input block')
        end select 
      enddo
_ELSE
      call caserr("no F90, so no van der Waals corrections")
_ENDIF
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine vdwaals_basic_control(ytest)
c
c     A very simple routine to control the basic options for using the
c     DFT-D dispersion corrections.
c     The only valid options are:
c     - on   - switch dispersion corrections on with default settings
c     - off  - switch dispersion corrections off
c     - aver - switch dispersion corrections on with the "average" 
c              c6 coefficient model
c     - geom - switch dispersion corrections on with the 
c              "geometric mean" c6 coefficient model
c
_IFN(f77)
      use vdwaals_corr
      implicit none
c
      character*4 ytest          ! keywords read in
      character*10 dft_xc        ! density functional
      integer :: ierr
c
c     Need to select the right pair model first to guarantee that
c     the correct functional is returned by vdw_get_functional.
c
      select case (ytest)
        case ('aver')                    ! c6model average
          call vdw_set_c6pair("average",ierr)
        case ('geom')                    ! c6model geometric mean
          call vdw_set_c6pair("geometric",ierr)
      end select 
      call vdw_get_functional(dft_xc)
      call vdw_set_defaults
      call vdw_switch_on
      select case (ytest)
        case ('on')                      ! on
          call vdw_switch_on 
        case ('off')                     ! off
          call vdw_switch_off 
        case ('aver')                    ! c6model average
          call vdw_set_c6pair("average",ierr)
        case ('geom')                    ! c6model geometric mean
          call vdw_set_c6pair("geometric",ierr)
        case default
          call caserr('unknown option van der Waals control')
      end select 
      call vdw_set_functional(dft_xc)
_ELSE
      call caserr("no F90, so no van der Waals corrections")
_ENDIF
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine vdwaals_functional(dft)
_IFN(f77)
      use vdwaals_corr
      implicit none
      character*(*) dft
c
      call vdw_set_functional(dft)
_ENDIF
      end
c
c-----------------------------------------------------------------------
c
      subroutine vdwaals_print
c
_IFN(f77)
      use vdwaals_corr
      implicit none
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
INCLUDE(common/iofile)
c
      integer maxelm
      parameter (maxelm = 118)
      logical oelm(maxelm) ! .true. if element present in molecule
      integer ii           ! counter
      character*90   ref   ! reference
      character*4    yelm  ! chemical symbol
      character*10   dft   ! density functional
      REAL c6, rad, scale, expo
c
      integer  isubst
      external isubst
      logical  opg_root
      external opg_root
c
      if (vdw_is_on().and.opg_root()) then
c
        call vdw_get_reference(ref)
        write(iwr,*)" ===============================================",
     +              "====================================="
        write(iwr,*)" Fudged van der Waals correction for effective ",
     +              "1-e models"
        write(iwr,*)" ",ref
        write(iwr,*)" -----------------------------------------------",
     +              "-------------------------------------"

        call vdw_get_functional(dft)
        call vdw_get_scale(scale)
        call vdw_get_alpha(expo)
        if (dft.eq."") then
          write(iwr,*)"  Scale       Exponent"
          write(iwr,'(1x,f7.2,f15.2)')scale,expo
        else
          write(iwr,*)"  Scale       Exponent   Functional"
          write(iwr,'(1x,f7.2,f15.2,3x,a)')scale,expo,dft
        endif
        write(iwr,*)" -----------------------------------------------",
     +              "-------------------------------------"
        do ii = 1, maxelm
          oelm(ii) = .false.
        enddo
        do ii = 1, nat
          oelm(isubst(zaname(ii))) = .true.
        enddo
        write(iwr,*)" Element    C6 (Hartree*Bohr^6)           R (Bohr)"
        do ii = 1, maxelm
          call vdw_get_c6coef(ii,c6)
          call vdw_get_radius(ii,rad)
          call ztoname(ii,yelm)
          if (oelm(ii).and.c6.gt.0.0d0) then
            write(iwr,'("  ",a4,8x,f18.2,f18.2)')yelm,c6,rad
          endif
        enddo
        write(iwr,*)" ===============================================",
     +              "====================================="
c
      endif
_ENDIF
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine ver_vdwaals(s,r,d)
c
_IFN(f77)
      use vdwaals_corr
      implicit none
c
      character*(*) s,r,d
c
      call vdw_get_version(s,r,d)
_ENDIF
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine vdwaals_energy(etot)
c
_IFN(f77)
      use vdwaals_corr
      implicit none
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
c
      REAL  etot 
      integer ielm(maxat) ! the atomic number of each atom
      integer ii          ! counter
c
      integer  isubst
      external isubst
c
      do ii = 1, nat
        ielm(ii) = isubst(zaname(ii))
      enddo
      call vdw_energy(nat,ielm,c,etot)
_ENDIF
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine vdwaals_gradient(grad)
_IFN(f77)
c
      use vdwaals_corr
      implicit none
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
c
      REAL  grad(3,maxat) ! the gradient 
      integer ielm(maxat) ! the atomic number of each atom
      integer ii          ! counter
c
      integer  isubst
      external isubst
c
      do ii = 1, nat
        ielm(ii) = isubst(zaname(ii))
      enddo
      call vdw_gradient(nat,ielm,c,grad)
_ENDIF
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine vdwaals_hessian(hess,ld)
c
_IFN(f77)
      use vdwaals_corr
      implicit none
c
INCLUDE(common/sizes)
INCLUDE(common/infoa)
INCLUDE(common/runlab)
c
      integer ld          ! leading dimension
      REAL  hess(ld,ld)   ! the hessian 
      integer ielm(maxat) ! the atomic number of each atom
      integer ii          ! counter
c
      integer  isubst
      external isubst
c
      do ii = 1, nat
        ielm(ii) = isubst(zaname(ii))
      enddo
      call vdw_hessian(nat,ielm,c,hess)
_ENDIF
c
      end
c
c     end van der Waals correction stuff
c
c=======================================================================
c
      subroutine ver_inter_vdwaals(s,r,d)
      character*80 source
      character*30 revision
      character*60 date
      character s*(*), r*(*), d*(*)
      data source /
     +     "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/inter_vdwaals.m,v $
     +     "/
      data revision /"$Revision: 5774 $"/
      data date /"$Date: 2008-12-05 00:26:07 +0100 (Fri, 05 Dec 2008) $
     +     "/
      s=source(9:)
      r=revision(11:)
      d=date(7:)
      return
      end

