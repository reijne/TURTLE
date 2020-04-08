!=======================================================================
!
!  Empirical dispersion/van der Waals corrections for effective
!  one electron models
!
!  [G]  "Accurate Description of van der Waals Complexes by Density
!        Functional Theory Including Empirical Corrections"
!       Stefan Grimme
!       Journal of Computational Chemistry 25 (2004) 1463-1473.
!
!  An alternative model for the pair C6 coefficients was suggested in
!
!  [AG] "Density functional theory including dispersion corrections for
!        intermolecular interactions in a large benchmark set of
!        biologically relevant molecules"
!       Jens Antony, Stefan Grimme
!       Physical Chemistry Chemical Physics 8 (2006) 5287-5293.
!
!=======================================================================
!
!  Written by Huub van Dam, STFC, May 2007.
!  Extended by Huub van Dam, STFC, May 2008:
!  - Added geometric mean C6 coefficient pair model.
!
!  Due to a user request to include empirical correction for the
!  dispersion energy in Hartree-Fock and DFT calculations this was
!  implemented. As requested we follow the paper by Grimme but any
!  user should realize that this is very much an approach to fudge the
!  results. The approach is implemented for the nuclear energy, the
!  gradient and the hessian.
!
!
!  The dispersion energy is approximated as (equation numbers starting
!  with G refer to equations in Grimme paper, equation numbers starting
!  with AG refer to the Antony-Grimme paper):
!
!     E    = - s6 sum(i=1:N,j=i+1:N) E_ij
!
!
!            C6_ij
!     E_ij = ------ * damp(R_ij)                                (eq. G2)
!            R_ij^6
!
!  where:
!     - N:    the number of atoms
!     - s6:   an overall scale factor (functional dependent)
!     - C6:   the C6 coefficient for an atom pair
!     - R:    the inter-atomic distance
!     - damp: a damping to suppress singularities at R=0
!
!
!  The damping function is given by:
!
!     damp(R) = [1+exp(-a{R/R0-1})]^(-1)                        (eq. G3)
!
!  where:
!     - a:    an exponent (value 23, see page 1465)
!     - R0:   the sum of atomic van der Waals radii
!
!
!  The C6 coefficient for a pair of atoms is approximated as either
!  an "average":
!
!                 C6_i * C6_j
!     C6_ij = 2 * -----------                                   (eq. G4)
!                 C6_i + C6_j
!
!  or a "geometric" mean:
!
!     C6_ij = sqrt(C6_i * C6_j)                              (eq. AG2.3)
!
!
!  The gradient of the above expression can be derived as follows:
!
!     dE_ij      C6_ij                 C6_ij    d damp(R_ij)
!     ----- = -6 ------ * damp(R_ij) + ------ * ------------     (eq. 5)
!     dR_ij      R_ij^7                R_ij^6   dR_ij
!
!
!  The second derivative becomes:
!
!     d2E_ij       C6_ij
!     ------- = 42 ------ * damp(R_ij)
!     dR_ij^2      R_ij^8
!
!                  C6_ij    d damp(R_ij)
!             - 12 ------ * ------------
!                  R_ij^7   dR_ij
!
!
!                  C6_ij    d2 damp(R_ij)
!             +    ------ * -------------                        (eq. 6)
!                  R_ij^6   dR_ij^2
!
!
!  The 1st derivative of the damping functions becomes:
!
!    d damp(R)    a * exp(-a{R/R0-1})
!    --------- =  --------------------------                     (eq. 7)
!    d R          R0 * [1+exp(-a{R/R0-1})]^2
!
!
!  The 2nd derivative becomes:
!
!    d^2 damp(R)    2 * a^2 * exp(-a{R/R0-1})^2
!    ----------- =  ----------------------------
!    d R^2          R0^2 * [1+exp(-a{R/R0-1})]^3
!
!                   a^2 * exp(-a{R/R0-1})
!                -  ----------------------------                 (eq. 8)
!                   R0^2 * [1+exp(-a{R/R0-1})]^2
!
!
!  To complete the derivative expressions the derivative of the inter-
!  atomic distances is needed as well. Obviously the inter-atomic
!  distance is defined as
!
!    R_ij = sqrt((x_i-x_j)^2+(y_i-y_j)^2+(y_i-y_j)^2)            (eq. 9)
!
!  Hence the gradients become:
!
!    d R_ij   (x_i-x_j)     d R_ij
!    ------ = --------- = - ------                              (eq. 10)
!    d x_i     R_ij         d x_j
!
!  The second derivative becomes
!
!    d^2 R_ij      d^2 R_ij        d^2 R_ij        d^2 R_ij
!    ----------- = ----------- = - ----------- = - ----------- =
!    d x_i^2       d x_j^2         d x_i d x_j     d x_j d x_i
!
!                    (x_i-x_j)^2    1
!                = - ----------- + ----                         (eq. 11)
!                    R_ij^3        R_ij
!
!    d^2 R_ij      d^2 R_ij        d^2 R_ij        d^2 R_ij
!    ----------- = ----------- = - ----------- = - ----------- =
!    d x_i d y_i   d x_j d y_j     d x_i d y_j     d x_j d y_i
!
!                    (x_i-x_j)(y_i-y_j)
!                = - ------------------                         (eq. 12)
!                    R_ij^3
!
!  Now the total gradient becomes
!
!    d E_ij   d E_ij   d R_ij
!    ------ = ------ * ------
!    d x_i    d R_ij   d x_i
!
!  The Hessian becomes
!
!    d2 E_ij        d2 E_ij    d R_ij   d R_ij   d E_ij   d2 R_ij
!    ----------- =  -------- * ------ * ------ + ------ * -----------
!    d x_i d y_i    d R_ij^2   d x_i    d y_i    d R_ij   d x_i d y_i
!
!  In the implementation it is recommended to keep the damping functions
!  as well as the calculation of the C6 coefficients separate as these
!  expressions vary greatly between different authors. 
!
!  Furthermore a modular design is followed for the implementation to
!  isolate this code as much as possible.
!
!  It is assumed that any code using this will handle all the input
!  processing as well as any printing. The reason is that most codes
!  have their own conventions in this respect.
!
module vdwaals_corr

  implicit none
  private 

  public :: vdw_set_defaults  ! set default atomic C6 coefficients and 
                              ! van der Waals radii

  public :: vdw_switch_on   ! turn van der Waals terms on
  public :: vdw_switch_off  ! turn van der Waals terms off
  public :: vdw_is_on       ! are the van der Waals terms turned on?

  public :: vdw_set_radius  ! set the van der Waals radius for an element
  public :: vdw_set_c6coef  ! set the C6 coefficient for an element
  public :: vdw_set_alpha   ! set the exponent in the damping function
  public :: vdw_set_scale   ! set the overall scale factor
  public :: vdw_set_c6pair  ! set the C6 pair coefficient model

  public :: vdw_get_radius  ! get the van der Waals radius of an element
  public :: vdw_get_c6coef  ! get the C6 coefficient of an element
  public :: vdw_get_alpha   ! get the exponent in the damping function
  public :: vdw_get_scale   ! get the overall scale factor
  public :: vdw_get_c6pair  ! get the C6 pair coefficient model

  public :: vdw_set_functional ! set the DFT functional (also sets scale)
  public :: vdw_get_functional ! get the DFT functional

  public :: vdw_get_reference ! return the literature reference for the
                              ! implementation of the van der Waals terms
  public :: vdw_get_version   ! return the version information

  public :: vdw_energy      ! compute the van der Waals energy correction
  public :: vdw_gradient    ! compute the van der Waals gradient correction
  public :: vdw_hessian     ! compute the van der Waals hessian correction
  
  integer, parameter :: ilr = selected_real_kind( 13 )

  integer, parameter :: c6pair_average   = 1
  integer, parameter :: c6pair_geometric = 2

  integer :: c6pair_model = c6pair_average

  logical :: initialised = .false.
  logical :: active      = .false.

  ! Has the vdw_set_scale routine been used? (false=no,true=yes)
  logical, dimension (1:2) :: set_scale = .false.

  ! The functional for each pair model
  character(len=10), dimension (1:2) :: functional = ""

  ! atomic C6 coefficients for each pair model
  real (ilr), dimension (-1:118,1:2) :: atom_c6

  ! atomic van der Waals radii for each pair model
  real (ilr), dimension (-1:118,1:2) :: vdw_rad

  ! damping function exponent for each pair model
  real (ilr), dimension (1:2) :: alpha

  ! overall scale factor for each pair model
  ! note: this should be dependent on the functional as well as the
  !       pair model.
  real (ilr), dimension (1:2) :: scale

contains


  pure real (ilr) function pair_c6( c6i, c6j)

    ! This function implements eq. G4. or AG2.3.

    real (ilr), intent( in ) :: c6i, c6j

    if (c6i == 0.0_ilr .or. c6j == 0.0_ilr) then

      ! The c6 coefficient for at least one atom is zero so suppress
      ! this pair contribution.

      pair_c6 = 0.0_ilr

    else

      ! Something sensible can be done

      select case (c6pair_model)

        case (c6pair_average)

          ! Use the normal average as in G4.

          pair_c6 = 2*c6i*c6j/(c6i+c6j)

        case (c6pair_geometric)

          ! Use the geometric mean as in AG2.3.

          pair_c6 = sqrt(c6i*c6j)

        case default

          ! The code has become very broken, ensure that clearly wrong
          ! answers result. Using the stop statement here could be
          ! problematic in parallel contexts, so do not do that.

          pair_c6 = 1.0e300_ilr

      end select

    endif

    return

  end function pair_c6


  pure real (ilr) function damp_exp( rr, aa, R0m1 )

    ! Utility function to evaluate the exponent part of eq. G3

    real (ilr), intent ( in ) :: rr   ! the inter atomic distance
    real (ilr), intent ( in ) :: aa   ! the scale factor
    real (ilr), intent ( in ) :: R0m1 ! the inverse of the sum of van 
                                      ! der Waals radii

    damp_exp = exp(-aa*(rr*R0m1-1.0_ilr))

    return

  end function damp_exp


  pure real (ilr) function damp( rr, aa, R0m1 )

    ! Evaluate the damping function of eq. G3

    real (ilr), intent ( in ) :: rr   ! the inter atomic distance
    real (ilr), intent ( in ) :: aa   ! the scale factor
    real (ilr), intent ( in ) :: R0m1 ! the inverse of the sum of van 
                                      ! der Waals radii

    damp = 1.0_ilr/(1.0_ilr+damp_exp( rr, aa, R0m1 ))

    return

  end function damp


  pure real (ilr) function d_damp( rr, aa, R0m1 )

    ! Evaluate the 1st derivative of damping function, i.e. eq. 7.

    real (ilr), intent ( in ) :: rr   ! the inter atomic distance
    real (ilr), intent ( in ) :: aa   ! the scale factor
    real (ilr), intent ( in ) :: R0m1 ! the inverse of the sum of van 
                                      ! der Waals radii

    d_damp = aa*R0m1*damp_exp(rr, aa, R0m1)*(damp(rr, aa, R0m1)**2)

    return

  end function d_damp


  pure real (ilr) function d2_damp( rr, aa, R0m1 )

    ! Evaluate the 2nd derivative of damping function, i.e. eq. 8.

    real (ilr), intent ( in ) :: rr   ! the inter atomic distance
    real (ilr), intent ( in ) :: aa   ! the scale factor
    real (ilr), intent ( in ) :: R0m1 ! the inverse of the sum of van 
                                      ! der Waals radii

    d2_damp = aa*R0m1*d_damp(rr, aa, R0m1) * &
              (2*damp_exp(rr, aa, R0m1)*damp(rr, aa, R0m1)-1.0_ilr)

    return

  end function d2_damp


  subroutine vdw_energy( natm, iatomno, coord, Etot )

    ! Evaluate van der Waals energy contribution according to eq. G2.
    ! This contribution is added onto Etot.

    integer,    intent ( in ) :: natm             ! the number of atoms
    integer,    intent ( in ) :: iatomno( natm )  ! the atomic numbers
                                                  ! of atoms in the
                                                  ! molecule
    real (ilr), intent ( in ) :: coord( 3, natm ) ! the atomic
                                                  ! coordinates

    real (ilr), intent ( inout ) :: Etot ! the total energy

    ! local variables:

    real (ilr) :: Evdw ! the van der Waals energy contribution
    real (ilr) :: c6ij ! the pair C6 coefficient
    real (ilr) :: rij  ! the inter atomic distance
    real (ilr) :: rij2 ! the inter atomic distance squared
    real (ilr) :: R0m1 ! the inverse of the sum of van der Waals radii
    real (ilr) :: Ep   ! the energy of a single pair

    integer    :: ii, jj, kk ! counters

    if (.not.active) return

    Evdw = 0.0_ilr

    do ii = 1, natm
      do jj = ii+1, natm
        c6ij = pair_c6( atom_c6(iatomno(ii),c6pair_model), &
                        atom_c6(iatomno(jj),c6pair_model) )
        if (c6ij > 0.0_ilr) then
          R0m1 = 1.0_ilr/(vdw_rad(iatomno(ii),c6pair_model)+ &
                          vdw_rad(iatomno(jj),c6pair_model))
          rij2 = 0.0_ilr
          do kk = 1, 3
            rij2 = rij2 + (coord(kk,ii)-coord(kk,jj))**2
          enddo
          rij = sqrt(rij2)
          Ep = c6ij*damp(rij,alpha(c6pair_model),R0m1)/(rij**6)
          Evdw = Evdw + Ep
        endif
      enddo
    enddo

    Etot = Etot - scale(c6pair_model) * Evdw

  end subroutine vdw_energy


  subroutine vdw_gradient( natm, iatomno, coord, grad )

    ! Evaluate van der Waals gradient contribution according to eq. 5
    ! and 10. This contribution is added onto grad.

    integer,    intent ( in ) :: natm             ! the number of atoms
    integer,    intent ( in ) :: iatomno( natm )  ! the atomic numbers
                                                  ! of atoms in the
                                                  ! molecule
    real (ilr), intent ( in ) :: coord( 3, natm ) ! the atomic
                                                  ! coordinates

    real (ilr), intent ( inout ) :: grad( 3, natm ) ! the gradient

    ! local variables:

    real (ilr) :: c6ij   ! the pair C6 coefficient
    real (ilr) :: rij    ! the inter atomic distance
    real (ilr) :: rij2   ! the inter atomic distance squared
    real (ilr) :: R0m1   ! the inverse of the sum of van der Waals radii
    real (ilr) :: xij(3) ! the coordinate differences
    real (ilr) :: dEi(3) ! the gradient with respect to the coordinates
                         ! of atom i.

    integer    :: ii, jj, kk ! counters

    if (.not.active) return

    do ii = 1, natm
      do jj = ii+1, natm
        c6ij = pair_c6( atom_c6(iatomno(ii),c6pair_model), &
                        atom_c6(iatomno(jj),c6pair_model) )
        if (c6ij > 0.0_ilr) then
          R0m1 = 1.0_ilr/(vdw_rad(iatomno(ii),c6pair_model)+ &
                          vdw_rad(iatomno(jj),c6pair_model))
          rij2 = 0.0_ilr
          do kk = 1, 3
            xij(kk) = coord(kk,ii)-coord(kk,jj)
            rij2 = rij2 + xij(kk)**2
          enddo
          rij = sqrt(rij2)
          do kk = 1, 3
            dEi(kk) = scale(c6pair_model) * (xij(kk)/rij) * &
            (-6.0_ilr*c6ij*damp(rij,alpha(c6pair_model),R0m1)/(rij**7) &
             + c6ij*d_damp(rij,alpha(c6pair_model),R0m1)/(rij**6))
          enddo
          do kk = 1, 3
            grad(kk,ii) = grad(kk,ii) - dEi(kk)
            grad(kk,jj) = grad(kk,jj) + dEi(kk)
          enddo
        endif
      enddo
    enddo

  end subroutine vdw_gradient


  subroutine vdw_hessian( natm, iatomno, coord, hess )

    ! Evaluate van der Waals Hessian contribution according to eq. 6,
    ! 11 and 12. This contribution is added onto hess.

    integer,    intent ( in ) :: natm             ! the number of atoms
    integer,    intent ( in ) :: iatomno( natm )  ! the atomic numbers
                                                  ! of atoms in the
                                                  ! molecule
    real (ilr), intent ( in ) :: coord( 3, natm ) ! the atomic
                                                  ! coordinates

    real (ilr), intent ( inout ) :: hess( 3*natm, 3*natm ) ! hessian

    ! local variables:

    real (ilr) :: c6ij        ! the pair C6 coefficient
    real (ilr) :: rij         ! the inter atomic distance
    real (ilr) :: rij2        ! the inter atomic distance squared
    real (ilr) :: R0m1        ! the inverse of the sum of van der Waals
                              ! radii
    real (ilr) :: xij(3)      ! the coordinate differences
    real (ilr) :: dEij        ! the 1st derivative wrt rij (not xi) 
                              ! (eq. 5)
    real (ilr) :: d2Eij       ! the 2nd derivative wrt rij (not xi) 
                              ! (eq. 6)
    real (ilr) :: d2Exii(3,3) ! as d2Eij but now wrt xi (eq. 11, 12)

    integer    :: ii, jj, kk, ll ! counters
    integer    :: ic, jc         ! scratch for coordinates

    if (.not.active) return

    do ii = 1, natm
      do jj = ii+1, natm
        c6ij = pair_c6( atom_c6(iatomno(ii),c6pair_model), &
                        atom_c6(iatomno(jj),c6pair_model) )
        if (c6ij > 0.0_ilr) then
          R0m1 = 1.0_ilr/(vdw_rad(iatomno(ii),c6pair_model)+ &
                          vdw_rad(iatomno(jj),c6pair_model))
          rij2 = 0.0_ilr
          do kk = 1, 3
            xij(kk) = coord(kk,ii)-coord(kk,jj)
            rij2 = rij2 + xij(kk)**2
          enddo
          rij = sqrt(rij2)
          d2Eij = scale(c6pair_model) * c6ij / rij**6 * &
                  (  42 * damp(rij,alpha(c6pair_model),R0m1)/rij2 &
                   - 12 * d_damp(rij,alpha(c6pair_model),R0m1)/rij &
                   +      d2_damp(rij,alpha(c6pair_model),R0m1))
          dEij  = scale(c6pair_model) * c6ij / rij**6 * &
                  (  -6 * damp(rij,alpha(c6pair_model),R0m1)/rij &
                    +     d_damp(rij,alpha(c6pair_model),R0m1))
          do kk = 1, 3
            d2Exii(kk,kk) = dEij*((xij(kk)**2)/rij2+1.0_ilr)/rij &
                          + d2Eij*(xij(kk)**2)/rij2
          enddo
          do kk = 1, 3
            do ll = kk+1, 3
              d2Exii(kk,ll) = dEij*xij(kk)*xij(ll)/(rij2*rij) &
                            + d2Eij*xij(kk)*xij(ll)/rij2
              d2Exii(ll,kk) = d2Exii(kk,ll)
            enddo
          enddo
          ic = (ii-1)*3
          jc = (jj-1)*3
          do kk = 1, 3
            do ll = 1, 3
              hess(ic+kk,ic+ll) = hess(ic+kk,ic+ll) - d2Exii(kk,ll)
              hess(jc+kk,jc+ll) = hess(jc+kk,jc+ll) - d2Exii(kk,ll)
              hess(ic+kk,jc+ll) = hess(ic+kk,jc+ll) + d2Exii(kk,ll)
              hess(jc+kk,ic+ll) = hess(jc+kk,ic+ll) + d2Exii(kk,ll)
            enddo
          enddo
        endif
      enddo
    enddo

  end subroutine vdw_hessian


  subroutine vdw_set_radius( ielement, radius )

    ! Modify the default van der Waals radius for atoms with atomic
    ! number ielement.

    integer,    intent( in ) :: ielement 
    real (ilr), intent( in ) :: radius

    vdw_rad(ielement,c6pair_model) = radius

  end subroutine vdw_set_radius


  subroutine vdw_get_radius( ielement, radius )

    ! Return the current van der Waals radius for atoms with atomic
    ! number ielement in radius.

    integer,    intent( in  ) :: ielement 
    real (ilr), intent( out ) :: radius

    radius = vdw_rad(ielement,c6pair_model)

  end subroutine vdw_get_radius


  subroutine vdw_set_c6coef( ielement, c6 )

    ! Modify the default C6 coefficient for atoms with atomic
    ! number ielement.

    integer,    intent( in ) :: ielement 
    real (ilr), intent( in ) :: c6

    atom_c6(ielement,c6pair_model) = c6

  end subroutine vdw_set_c6coef


  subroutine vdw_get_c6coef( ielement, c6 )

    ! Return the current C6 coefficient for atoms with atomic
    ! number ielement in c6.

    integer,    intent( in  ) :: ielement 
    real (ilr), intent( out ) :: c6

    c6 = atom_c6(ielement,c6pair_model)

  end subroutine vdw_get_c6coef


  subroutine vdw_set_alpha( aa )

    ! Modify the default damping exponent.

    real (ilr), intent( in ) :: aa

    alpha(c6pair_model) = aa

  end subroutine vdw_set_alpha


  subroutine vdw_get_alpha( aa )

    ! Return the current damping exponent in aa.

    real (ilr), intent( out ) :: aa

    aa = alpha(c6pair_model)

  end subroutine vdw_get_alpha


  subroutine vdw_set_scale( s6 )

    ! Modify the default overall scaling factor.

    real (ilr), intent( in ) :: s6

    scale(c6pair_model) = s6
    set_scale(c6pair_model) = .true.
    functional(c6pair_model) = ""

  end subroutine vdw_set_scale


  subroutine vdw_get_scale( s6 )

    ! Return the current overall scaling factor in s6.

    real (ilr), intent( out ) :: s6

    s6 = scale(c6pair_model)

  end subroutine vdw_get_scale


  subroutine vdw_set_functional( dft )

    ! Set the chosen functional and the corresponding scale factors.
    ! This subroutine is just meant to be an aid for selecting a set
    ! up that is fully consistent with the literature. To avoid clashes
    ! with explicitly specified scale factors the values set with
    ! vdw_set_scale take precedence.

    character(len=*), intent( in ) :: dft

    if (.not.set_scale(c6pair_average)) then

      functional(c6pair_average)   = dft

      select case (dft)

        case ("pbe")

          scale(c6pair_average)        = 0.70_ilr

        case ("blyp")
 
          scale(c6pair_average)        = 1.40_ilr

        case ("b-p86","bp86")

          scale(c6pair_average)        = 1.30_ilr

        case default

          functional(c6pair_average)   = ""
          scale(c6pair_average)        = 1.00_ilr

     end select

    endif

    if (.not.set_scale(c6pair_geometric)) then

      functional(c6pair_geometric) = dft

      select case (dft)

        case ("pbe")

          scale(c6pair_geometric)      = 0.75_ilr

        case ("blyp")
 
          scale(c6pair_geometric)      = 1.20_ilr

        case ("b-p86","bp86")

          scale(c6pair_geometric)      = 1.05_ilr

        case ("tpss")

          scale(c6pair_geometric)      = 1.00_ilr

        case ("b3lyp")

          scale(c6pair_geometric)      = 1.05_ilr

        case ("b97-d")

          scale(c6pair_geometric)      = 1.25_ilr

        case default

          functional(c6pair_geometric) = ""
          scale(c6pair_geometric)      = 1.00_ilr

     end select

    endif

  end subroutine vdw_set_functional


  subroutine vdw_get_functional( dft )

    ! Return the chosen functional for the current active pair model.
    ! "" is returned if no functional with known parameters has been
    !    selected or if no functional was selected.

    character(len=*), intent( out ) :: dft

    dft = functional(c6pair_model)

  end subroutine vdw_get_functional


  subroutine vdw_set_c6pair( c6model, ierror )

    ! Modify the default C6 pair coefficient model,
    ! ierror is incremented if an invalid model is specified.

    character(len=*), intent( in )    :: c6model 
    integer,          intent( inout ) :: ierror

    select case (c6model)

      case ("average")
   
        c6pair_model = c6pair_average

      case ("geometric") 
   
        c6pair_model = c6pair_geometric

      case default

        ierror = ierror + 1

    end select 

  end subroutine vdw_set_c6pair


  subroutine vdw_get_c6pair( c6model )

    ! Return the current C6 pair coefficient model.

    character(len=*), intent( out ) :: c6model 

    select case (c6pair_model)

      case (c6pair_average)
   
        c6model = "average"

      case (c6pair_geometric) 
   
        c6model = "geometric"

      case default

        c6model = "ERROR in VdW"

    end select 

  end subroutine vdw_get_c6pair


  subroutine vdw_switch_on

    ! Turn the van der Waals contributions on

    active = .true.

  end subroutine vdw_switch_on


  subroutine vdw_switch_off

    ! Turn the van der Waals contributions off

    active = .false.

  end subroutine vdw_switch_off


  logical function vdw_is_on()

    ! Return .true. if van der Waals contributions are turned on
    ! Return .false. otherwise

    vdw_is_on = active.and.initialised

  end function vdw_is_on


  subroutine vdw_get_reference( reference )

    ! Return the literature reference in a short hand notation

    character*(*), intent ( out ) :: reference

    character(len=*), parameter :: ref_G = &
    "S. Grimme, J.Comput.Chem. 25 (2004) 1463-1473,&
    & doi:10.1002/jcc.20078"

    character(len=*), parameter :: ref_AG = &
    "J. Antony, S. Grimme, Phys.Chem.Chem.Phys. 8 (2006) 5287-5293,&
    & doi:10.1039/b612585a"

    integer :: lref_G, lref_AG, lreference ! length of character buffers

    lref_G     = len(ref_G)
    lref_AG    = len(ref_AG)
    lreference = len(reference)

    select case (c6pair_model)

      case (c6pair_average)

        if (lreference < lref_G) then
          reference = ref_G(1:lreference)
        else
          reference = ref_G
        endif

      case (c6pair_geometric)

        if (lreference < lref_AG) then
          reference = ref_AG(1:lreference)
        else
          reference = ref_AG
        endif

      case default

        reference = "ERROR in VdW"

    end select

  end subroutine vdw_get_reference


  subroutine vdw_get_version(s,r,d)

    ! Return the version information

    character*(*), intent ( out ) :: s, r, d

    character(len=*), parameter :: source = &
    "$Source: /c/qcg/cvs/psh/GAMESS-UK/m4/vdwaals_corr.f90,v $"
    character(len=*), parameter :: revision = &
    "$Revision: 1.2 $"
    character(len=*), parameter :: date = &
    "$Date: 2007-05-24 14:13:42 $"

    s=source(9:)
    r=revision(11:)
    d=date(7:)

  end subroutine vdw_get_version


  subroutine vdw_set_defaults

    ! Set the default parameters for the 
    ! 1 Atomic C6 coefficients
    ! 2 van der Waals radii
    ! 3 damping function exponent
    ! 4 overall scaling function
    ! 5 the C6 pair coefficient model
    ! Also tick that the stuff has been initialised. These setting may
    ! subsequently be overridden by the user.

    real (ilr), parameter :: ang_to_bohr = 1.0_ilr/0.529177249_ilr
    real (ilr), parameter :: nm_to_bohr = 10.0_ilr*ang_to_bohr
    real (ilr), parameter :: jmol_to_hartree = 1.0_ilr/2625.50e3_ilr

    integer :: i

    ! Set parameters up for the "average" pair model
    ! following Grimme [G]

    i = c6pair_average

    alpha(i) = 23.0_ilr ! see Grimme [G]
    scale(i) =  1.0_ilr ! for now

    ! van der Waals radii in au (also from Grimme [G])

    vdw_rad( -1,i) = 0.0_ilr
    vdw_rad(  0,i) = 0.0_ilr
    vdw_rad(  1,i) = 1.11_ilr * ang_to_bohr
    vdw_rad(  2,i) = 0.0_ilr
    vdw_rad(  3,i) = 0.0_ilr
    vdw_rad(  4,i) = 0.0_ilr
    vdw_rad(  5,i) = 0.0_ilr
    vdw_rad(  6,i) = 1.61_ilr * ang_to_bohr
    vdw_rad(  7,i) = 1.55_ilr * ang_to_bohr
    vdw_rad(  8,i) = 1.49_ilr * ang_to_bohr
    vdw_rad(  9,i) = 1.43_ilr * ang_to_bohr
    vdw_rad( 10,i) = 1.38_ilr * ang_to_bohr
    vdw_rad( 11,i) = 0.0_ilr
    vdw_rad( 12,i) = 0.0_ilr
    vdw_rad( 13,i) = 0.0_ilr
    vdw_rad( 14,i) = 0.0_ilr
    vdw_rad( 15,i) = 0.0_ilr
    vdw_rad( 16,i) = 0.0_ilr
    vdw_rad( 17,i) = 0.0_ilr
    vdw_rad( 18,i) = 0.0_ilr
    vdw_rad( 19,i) = 0.0_ilr
    vdw_rad( 20,i) = 0.0_ilr
    vdw_rad( 21,i) = 0.0_ilr
    vdw_rad( 22,i) = 0.0_ilr
    vdw_rad( 23,i) = 0.0_ilr
    vdw_rad( 24,i) = 0.0_ilr
    vdw_rad( 25,i) = 0.0_ilr
    vdw_rad( 26,i) = 0.0_ilr
    vdw_rad( 27,i) = 0.0_ilr
    vdw_rad( 28,i) = 0.0_ilr
    vdw_rad( 29,i) = 0.0_ilr
    vdw_rad( 30,i) = 0.0_ilr
    vdw_rad( 31,i) = 0.0_ilr
    vdw_rad( 32,i) = 0.0_ilr
    vdw_rad( 33,i) = 0.0_ilr
    vdw_rad( 34,i) = 0.0_ilr
    vdw_rad( 35,i) = 0.0_ilr
    vdw_rad( 36,i) = 0.0_ilr
    vdw_rad( 37,i) = 0.0_ilr
    vdw_rad( 38,i) = 0.0_ilr
    vdw_rad( 39,i) = 0.0_ilr
    vdw_rad( 40,i) = 0.0_ilr
    vdw_rad( 41,i) = 0.0_ilr
    vdw_rad( 42,i) = 0.0_ilr
    vdw_rad( 43,i) = 0.0_ilr
    vdw_rad( 44,i) = 0.0_ilr
    vdw_rad( 45,i) = 0.0_ilr
    vdw_rad( 46,i) = 0.0_ilr
    vdw_rad( 47,i) = 0.0_ilr
    vdw_rad( 48,i) = 0.0_ilr
    vdw_rad( 49,i) = 0.0_ilr
    vdw_rad( 50,i) = 0.0_ilr
    vdw_rad( 51,i) = 0.0_ilr
    vdw_rad( 52,i) = 0.0_ilr
    vdw_rad( 53,i) = 0.0_ilr
    vdw_rad( 54,i) = 0.0_ilr
    vdw_rad( 55,i) = 0.0_ilr
    vdw_rad( 56,i) = 0.0_ilr
    vdw_rad( 57,i) = 0.0_ilr
    vdw_rad( 58,i) = 0.0_ilr
    vdw_rad( 59,i) = 0.0_ilr
    vdw_rad( 60,i) = 0.0_ilr
    vdw_rad( 61,i) = 0.0_ilr
    vdw_rad( 62,i) = 0.0_ilr
    vdw_rad( 63,i) = 0.0_ilr
    vdw_rad( 64,i) = 0.0_ilr
    vdw_rad( 65,i) = 0.0_ilr
    vdw_rad( 66,i) = 0.0_ilr
    vdw_rad( 67,i) = 0.0_ilr
    vdw_rad( 68,i) = 0.0_ilr
    vdw_rad( 69,i) = 0.0_ilr
    vdw_rad( 70,i) = 0.0_ilr
    vdw_rad( 71,i) = 0.0_ilr
    vdw_rad( 72,i) = 0.0_ilr
    vdw_rad( 73,i) = 0.0_ilr
    vdw_rad( 74,i) = 0.0_ilr
    vdw_rad( 75,i) = 0.0_ilr
    vdw_rad( 76,i) = 0.0_ilr
    vdw_rad( 77,i) = 0.0_ilr
    vdw_rad( 78,i) = 0.0_ilr
    vdw_rad( 79,i) = 0.0_ilr
    vdw_rad( 80,i) = 0.0_ilr
    vdw_rad( 81,i) = 0.0_ilr
    vdw_rad( 82,i) = 0.0_ilr
    vdw_rad( 83,i) = 0.0_ilr
    vdw_rad( 84,i) = 0.0_ilr
    vdw_rad( 85,i) = 0.0_ilr
    vdw_rad( 86,i) = 0.0_ilr
    vdw_rad( 87,i) = 0.0_ilr
    vdw_rad( 88,i) = 0.0_ilr
    vdw_rad( 89,i) = 0.0_ilr
    vdw_rad( 90,i) = 0.0_ilr
    vdw_rad( 91,i) = 0.0_ilr
    vdw_rad( 92,i) = 0.0_ilr
    vdw_rad( 93,i) = 0.0_ilr
    vdw_rad( 94,i) = 0.0_ilr
    vdw_rad( 95,i) = 0.0_ilr
    vdw_rad( 96,i) = 0.0_ilr
    vdw_rad( 97,i) = 0.0_ilr
    vdw_rad( 98,i) = 0.0_ilr
    vdw_rad( 99,i) = 0.0_ilr
    vdw_rad(100,i) = 0.0_ilr
    vdw_rad(101,i) = 0.0_ilr
    vdw_rad(102,i) = 0.0_ilr
    vdw_rad(103,i) = 0.0_ilr
    vdw_rad(104,i) = 0.0_ilr
    vdw_rad(105,i) = 0.0_ilr
    vdw_rad(106,i) = 0.0_ilr
    vdw_rad(107,i) = 0.0_ilr
    vdw_rad(108,i) = 0.0_ilr
    vdw_rad(109,i) = 0.0_ilr
    vdw_rad(110,i) = 0.0_ilr
    vdw_rad(111,i) = 0.0_ilr
    vdw_rad(112,i) = 0.0_ilr
    vdw_rad(113,i) = 0.0_ilr
    vdw_rad(114,i) = 0.0_ilr
    vdw_rad(115,i) = 0.0_ilr
    vdw_rad(116,i) = 0.0_ilr
    vdw_rad(117,i) = 0.0_ilr
    vdw_rad(118,i) = 0.0_ilr

    ! van der Waals radii in au (also from Grimme [G])

    atom_c6( -1,i) = 0.0_ilr
    atom_c6(  0,i) = 0.0_ilr
    atom_c6(  1,i) = 0.16_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  2,i) = 0.0_ilr
    atom_c6(  3,i) = 0.0_ilr
    atom_c6(  4,i) = 0.0_ilr
    atom_c6(  5,i) = 0.0_ilr
    atom_c6(  6,i) = 1.65_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  7,i) = 1.11_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  8,i) = 0.70_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  9,i) = 0.57_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 10,i) = 0.45_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 11,i) = 0.0_ilr
    atom_c6( 12,i) = 0.0_ilr
    atom_c6( 13,i) = 0.0_ilr
    atom_c6( 14,i) = 0.0_ilr
    atom_c6( 15,i) = 0.0_ilr
    atom_c6( 16,i) = 0.0_ilr
    atom_c6( 17,i) = 0.0_ilr
    atom_c6( 18,i) = 0.0_ilr
    atom_c6( 19,i) = 0.0_ilr
    atom_c6( 20,i) = 0.0_ilr
    atom_c6( 21,i) = 0.0_ilr
    atom_c6( 22,i) = 0.0_ilr
    atom_c6( 23,i) = 0.0_ilr
    atom_c6( 24,i) = 0.0_ilr
    atom_c6( 25,i) = 0.0_ilr
    atom_c6( 26,i) = 0.0_ilr
    atom_c6( 27,i) = 0.0_ilr
    atom_c6( 28,i) = 0.0_ilr
    atom_c6( 29,i) = 0.0_ilr
    atom_c6( 30,i) = 0.0_ilr
    atom_c6( 31,i) = 0.0_ilr
    atom_c6( 32,i) = 0.0_ilr
    atom_c6( 33,i) = 0.0_ilr
    atom_c6( 34,i) = 0.0_ilr
    atom_c6( 35,i) = 0.0_ilr
    atom_c6( 36,i) = 0.0_ilr
    atom_c6( 37,i) = 0.0_ilr
    atom_c6( 38,i) = 0.0_ilr
    atom_c6( 39,i) = 0.0_ilr
    atom_c6( 40,i) = 0.0_ilr
    atom_c6( 41,i) = 0.0_ilr
    atom_c6( 42,i) = 0.0_ilr
    atom_c6( 43,i) = 0.0_ilr
    atom_c6( 44,i) = 0.0_ilr
    atom_c6( 45,i) = 0.0_ilr
    atom_c6( 46,i) = 0.0_ilr
    atom_c6( 47,i) = 0.0_ilr
    atom_c6( 48,i) = 0.0_ilr
    atom_c6( 49,i) = 0.0_ilr
    atom_c6( 50,i) = 0.0_ilr
    atom_c6( 51,i) = 0.0_ilr
    atom_c6( 52,i) = 0.0_ilr
    atom_c6( 53,i) = 0.0_ilr
    atom_c6( 54,i) = 0.0_ilr
    atom_c6( 55,i) = 0.0_ilr
    atom_c6( 56,i) = 0.0_ilr
    atom_c6( 57,i) = 0.0_ilr
    atom_c6( 58,i) = 0.0_ilr
    atom_c6( 59,i) = 0.0_ilr
    atom_c6( 60,i) = 0.0_ilr
    atom_c6( 61,i) = 0.0_ilr
    atom_c6( 62,i) = 0.0_ilr
    atom_c6( 63,i) = 0.0_ilr
    atom_c6( 64,i) = 0.0_ilr
    atom_c6( 65,i) = 0.0_ilr
    atom_c6( 66,i) = 0.0_ilr
    atom_c6( 67,i) = 0.0_ilr
    atom_c6( 68,i) = 0.0_ilr
    atom_c6( 69,i) = 0.0_ilr
    atom_c6( 70,i) = 0.0_ilr
    atom_c6( 71,i) = 0.0_ilr
    atom_c6( 72,i) = 0.0_ilr
    atom_c6( 73,i) = 0.0_ilr
    atom_c6( 74,i) = 0.0_ilr
    atom_c6( 75,i) = 0.0_ilr
    atom_c6( 76,i) = 0.0_ilr
    atom_c6( 77,i) = 0.0_ilr
    atom_c6( 78,i) = 0.0_ilr
    atom_c6( 79,i) = 0.0_ilr
    atom_c6( 80,i) = 0.0_ilr
    atom_c6( 81,i) = 0.0_ilr
    atom_c6( 82,i) = 0.0_ilr
    atom_c6( 83,i) = 0.0_ilr
    atom_c6( 84,i) = 0.0_ilr
    atom_c6( 85,i) = 0.0_ilr
    atom_c6( 86,i) = 0.0_ilr
    atom_c6( 87,i) = 0.0_ilr
    atom_c6( 88,i) = 0.0_ilr
    atom_c6( 89,i) = 0.0_ilr
    atom_c6( 90,i) = 0.0_ilr
    atom_c6( 91,i) = 0.0_ilr
    atom_c6( 92,i) = 0.0_ilr
    atom_c6( 93,i) = 0.0_ilr
    atom_c6( 94,i) = 0.0_ilr
    atom_c6( 95,i) = 0.0_ilr
    atom_c6( 96,i) = 0.0_ilr
    atom_c6( 97,i) = 0.0_ilr
    atom_c6( 98,i) = 0.0_ilr
    atom_c6( 99,i) = 0.0_ilr
    atom_c6(100,i) = 0.0_ilr
    atom_c6(101,i) = 0.0_ilr
    atom_c6(102,i) = 0.0_ilr
    atom_c6(103,i) = 0.0_ilr
    atom_c6(104,i) = 0.0_ilr
    atom_c6(105,i) = 0.0_ilr
    atom_c6(106,i) = 0.0_ilr
    atom_c6(107,i) = 0.0_ilr
    atom_c6(108,i) = 0.0_ilr
    atom_c6(109,i) = 0.0_ilr
    atom_c6(110,i) = 0.0_ilr
    atom_c6(111,i) = 0.0_ilr
    atom_c6(112,i) = 0.0_ilr
    atom_c6(113,i) = 0.0_ilr
    atom_c6(114,i) = 0.0_ilr
    atom_c6(115,i) = 0.0_ilr
    atom_c6(116,i) = 0.0_ilr
    atom_c6(117,i) = 0.0_ilr
    atom_c6(118,i) = 0.0_ilr


    ! Set parameters up for the "average" pair model
    ! following Antony-Grimme [AG]

    i = c6pair_geometric

    alpha(i) = 20.0_ilr ! see Antony-Grimme [AG]
    scale(i) =  1.0_ilr ! for now

    ! van der Waals radii in au (from Antony-Grimme [AG])

    vdw_rad( -1,i) = 0.0_ilr
    vdw_rad(  0,i) = 0.0_ilr
    vdw_rad(  1,i) = 1.001_ilr * ang_to_bohr
    vdw_rad(  2,i) = 1.012_ilr * ang_to_bohr
    vdw_rad(  3,i) = 0.825_ilr * ang_to_bohr
    vdw_rad(  4,i) = 1.408_ilr * ang_to_bohr
    vdw_rad(  5,i) = 1.485_ilr * ang_to_bohr
    vdw_rad(  6,i) = 1.452_ilr * ang_to_bohr
    vdw_rad(  7,i) = 1.397_ilr * ang_to_bohr
    vdw_rad(  8,i) = 1.342_ilr * ang_to_bohr
    vdw_rad(  9,i) = 1.287_ilr * ang_to_bohr
    vdw_rad( 10,i) = 1.243_ilr * ang_to_bohr
    vdw_rad( 11,i) = 1.144_ilr * ang_to_bohr
    vdw_rad( 12,i) = 1.364_ilr * ang_to_bohr
    vdw_rad( 13,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 14,i) = 1.716_ilr * ang_to_bohr
    vdw_rad( 15,i) = 1.705_ilr * ang_to_bohr
    vdw_rad( 16,i) = 1.683_ilr * ang_to_bohr
    vdw_rad( 17,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 18,i) = 1.595_ilr * ang_to_bohr
    vdw_rad( 19,i) = 1.485_ilr * ang_to_bohr
    vdw_rad( 20,i) = 1.474_ilr * ang_to_bohr
    vdw_rad( 21,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 22,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 23,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 24,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 25,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 26,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 27,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 28,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 29,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 30,i) = 1.562_ilr * ang_to_bohr
    vdw_rad( 31,i) = 1.650_ilr * ang_to_bohr
    vdw_rad( 32,i) = 1.727_ilr * ang_to_bohr
    vdw_rad( 33,i) = 1.760_ilr * ang_to_bohr
    vdw_rad( 34,i) = 1.771_ilr * ang_to_bohr
    vdw_rad( 35,i) = 1.749_ilr * ang_to_bohr
    vdw_rad( 36,i) = 1.727_ilr * ang_to_bohr
    vdw_rad( 37,i) = 1.628_ilr * ang_to_bohr
    vdw_rad( 38,i) = 1.606_ilr * ang_to_bohr
    vdw_rad( 39,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 40,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 41,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 42,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 43,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 44,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 45,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 46,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 47,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 48,i) = 1.639_ilr * ang_to_bohr
    vdw_rad( 49,i) = 1.672_ilr * ang_to_bohr
    vdw_rad( 50,i) = 1.804_ilr * ang_to_bohr
    vdw_rad( 51,i) = 1.881_ilr * ang_to_bohr
    vdw_rad( 52,i) = 1.892_ilr * ang_to_bohr
    vdw_rad( 53,i) = 1.892_ilr * ang_to_bohr
    vdw_rad( 54,i) = 1.881_ilr * ang_to_bohr
    vdw_rad( 55,i) = 0.0_ilr
    vdw_rad( 56,i) = 0.0_ilr
    vdw_rad( 57,i) = 0.0_ilr
    vdw_rad( 58,i) = 0.0_ilr
    vdw_rad( 59,i) = 0.0_ilr
    vdw_rad( 60,i) = 0.0_ilr
    vdw_rad( 61,i) = 0.0_ilr
    vdw_rad( 62,i) = 0.0_ilr
    vdw_rad( 63,i) = 0.0_ilr
    vdw_rad( 64,i) = 0.0_ilr
    vdw_rad( 65,i) = 0.0_ilr
    vdw_rad( 66,i) = 0.0_ilr
    vdw_rad( 67,i) = 0.0_ilr
    vdw_rad( 68,i) = 0.0_ilr
    vdw_rad( 69,i) = 0.0_ilr
    vdw_rad( 70,i) = 0.0_ilr
    vdw_rad( 71,i) = 0.0_ilr
    vdw_rad( 72,i) = 0.0_ilr
    vdw_rad( 73,i) = 0.0_ilr
    vdw_rad( 74,i) = 0.0_ilr
    vdw_rad( 75,i) = 0.0_ilr
    vdw_rad( 76,i) = 0.0_ilr
    vdw_rad( 77,i) = 0.0_ilr
    vdw_rad( 78,i) = 0.0_ilr
    vdw_rad( 79,i) = 0.0_ilr
    vdw_rad( 80,i) = 0.0_ilr
    vdw_rad( 81,i) = 0.0_ilr
    vdw_rad( 82,i) = 0.0_ilr
    vdw_rad( 83,i) = 0.0_ilr
    vdw_rad( 84,i) = 0.0_ilr
    vdw_rad( 85,i) = 0.0_ilr
    vdw_rad( 86,i) = 0.0_ilr
    vdw_rad( 87,i) = 0.0_ilr
    vdw_rad( 88,i) = 0.0_ilr
    vdw_rad( 89,i) = 0.0_ilr
    vdw_rad( 90,i) = 0.0_ilr
    vdw_rad( 91,i) = 0.0_ilr
    vdw_rad( 92,i) = 0.0_ilr
    vdw_rad( 93,i) = 0.0_ilr
    vdw_rad( 94,i) = 0.0_ilr
    vdw_rad( 95,i) = 0.0_ilr
    vdw_rad( 96,i) = 0.0_ilr
    vdw_rad( 97,i) = 0.0_ilr
    vdw_rad( 98,i) = 0.0_ilr
    vdw_rad( 99,i) = 0.0_ilr
    vdw_rad(100,i) = 0.0_ilr
    vdw_rad(101,i) = 0.0_ilr
    vdw_rad(102,i) = 0.0_ilr
    vdw_rad(103,i) = 0.0_ilr
    vdw_rad(104,i) = 0.0_ilr
    vdw_rad(105,i) = 0.0_ilr
    vdw_rad(106,i) = 0.0_ilr
    vdw_rad(107,i) = 0.0_ilr
    vdw_rad(108,i) = 0.0_ilr
    vdw_rad(109,i) = 0.0_ilr
    vdw_rad(110,i) = 0.0_ilr
    vdw_rad(111,i) = 0.0_ilr
    vdw_rad(112,i) = 0.0_ilr
    vdw_rad(113,i) = 0.0_ilr
    vdw_rad(114,i) = 0.0_ilr
    vdw_rad(115,i) = 0.0_ilr
    vdw_rad(116,i) = 0.0_ilr
    vdw_rad(117,i) = 0.0_ilr
    vdw_rad(118,i) = 0.0_ilr

    ! van der Waals radii in au (also from Antony-Grimme [AG])

    atom_c6( -1,i) =  0.0_ilr
    atom_c6(  0,i) =  0.0_ilr
    atom_c6(  1,i) =  0.14_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  2,i) =  0.08_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  3,i) =  1.61_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  4,i) =  1.61_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  5,i) =  3.13_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  6,i) =  1.75_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  7,i) =  1.23_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  8,i) =  0.70_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6(  9,i) =  0.75_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 10,i) =  0.63_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 11,i) =  5.71_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 12,i) =  5.71_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 13,i) = 10.79_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 14,i) =  9.23_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 15,i) =  7.84_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 16,i) =  5.57_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 17,i) =  5.07_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 18,i) =  4.61_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 19,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 20,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 21,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 22,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 23,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 24,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 25,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 26,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 27,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 28,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 29,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 30,i) = 10.80_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 31,i) = 16.99_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 32,i) = 17.10_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 33,i) = 16.37_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 34,i) = 12.64_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 35,i) = 12.47_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 36,i) = 12.01_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 37,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 38,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 39,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 40,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 41,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 42,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 43,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 44,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 45,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 46,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 47,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 48,i) = 24.67_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 49,i) = 37.32_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 50,i) = 38.71_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 51,i) = 38.44_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 52,i) = 31.74_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 53,i) = 31.50_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 54,i) = 29.99_ilr * jmol_to_hartree * ( nm_to_bohr ** 6 )
    atom_c6( 55,i) =  0.0_ilr
    atom_c6( 56,i) =  0.0_ilr
    atom_c6( 57,i) =  0.0_ilr
    atom_c6( 58,i) =  0.0_ilr
    atom_c6( 59,i) =  0.0_ilr
    atom_c6( 60,i) =  0.0_ilr
    atom_c6( 61,i) =  0.0_ilr
    atom_c6( 62,i) =  0.0_ilr
    atom_c6( 63,i) =  0.0_ilr
    atom_c6( 64,i) =  0.0_ilr
    atom_c6( 65,i) =  0.0_ilr
    atom_c6( 66,i) =  0.0_ilr
    atom_c6( 67,i) =  0.0_ilr
    atom_c6( 68,i) =  0.0_ilr
    atom_c6( 69,i) =  0.0_ilr
    atom_c6( 70,i) =  0.0_ilr
    atom_c6( 71,i) =  0.0_ilr
    atom_c6( 72,i) =  0.0_ilr
    atom_c6( 73,i) =  0.0_ilr
    atom_c6( 74,i) =  0.0_ilr
    atom_c6( 75,i) =  0.0_ilr
    atom_c6( 76,i) =  0.0_ilr
    atom_c6( 77,i) =  0.0_ilr
    atom_c6( 78,i) =  0.0_ilr
    atom_c6( 79,i) =  0.0_ilr
    atom_c6( 80,i) =  0.0_ilr
    atom_c6( 81,i) =  0.0_ilr
    atom_c6( 82,i) =  0.0_ilr
    atom_c6( 83,i) =  0.0_ilr
    atom_c6( 84,i) =  0.0_ilr
    atom_c6( 85,i) =  0.0_ilr
    atom_c6( 86,i) =  0.0_ilr
    atom_c6( 87,i) =  0.0_ilr
    atom_c6( 88,i) =  0.0_ilr
    atom_c6( 89,i) =  0.0_ilr
    atom_c6( 90,i) =  0.0_ilr
    atom_c6( 91,i) =  0.0_ilr
    atom_c6( 92,i) =  0.0_ilr
    atom_c6( 93,i) =  0.0_ilr
    atom_c6( 94,i) =  0.0_ilr
    atom_c6( 95,i) =  0.0_ilr
    atom_c6( 96,i) =  0.0_ilr
    atom_c6( 97,i) =  0.0_ilr
    atom_c6( 98,i) =  0.0_ilr
    atom_c6( 99,i) =  0.0_ilr
    atom_c6(100,i) =  0.0_ilr
    atom_c6(101,i) =  0.0_ilr
    atom_c6(102,i) =  0.0_ilr
    atom_c6(103,i) =  0.0_ilr
    atom_c6(104,i) =  0.0_ilr
    atom_c6(105,i) =  0.0_ilr
    atom_c6(106,i) =  0.0_ilr
    atom_c6(107,i) =  0.0_ilr
    atom_c6(108,i) =  0.0_ilr
    atom_c6(109,i) =  0.0_ilr
    atom_c6(110,i) =  0.0_ilr
    atom_c6(111,i) =  0.0_ilr
    atom_c6(112,i) =  0.0_ilr
    atom_c6(113,i) =  0.0_ilr
    atom_c6(114,i) =  0.0_ilr
    atom_c6(115,i) =  0.0_ilr
    atom_c6(116,i) =  0.0_ilr
    atom_c6(117,i) =  0.0_ilr
    atom_c6(118,i) =  0.0_ilr

    c6pair_model = c6pair_average

    initialised = .true.
   
  end subroutine vdw_set_defaults

end module vdwaals_corr
