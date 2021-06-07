!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module theorywaves

!BOP
!\newpage
! !MODULE: theorywaves
!
! !AUTHOR:
!  
!
! !DESCRIPTION:
!\\
!\\
!  References:\\
!\\
!\\

! !USES:

  use cvmix_kinds_and_types, only : cvmix_r8,                                 &
                                    cvmix_strlen,                             &
                                    cvmix_zero,                               &
                                    cvmix_one,                                &
                                    cvmix_PI,                                 &
                                    cvmix_data_type,                          &
                                    cvmix_global_params_type,                 &
                                    CVMIX_OVERWRITE_OLD_VAL,                  &
                                    CVMIX_SUM_OLD_AND_NEW_VALS,               &
                                    CVMIX_MAX_OLD_AND_NEW_VALS
  use cvmix_math, only :            CVMIX_MATH_INTERP_LINEAR,                 &
                                    CVMIX_MATH_INTERP_QUAD,                   &
                                    CVMIX_MATH_INTERP_CUBE_SPLINE,            &
                                    cvmix_math_poly_interp,                   &
                                    cvmix_math_cubic_root_find,               &
                                    cvmix_math_evaluate_cubic
  use cvmix_put_get,         only : cvmix_put
  use cvmix_utils,           only : cvmix_update_wrap

!EOP

  implicit none
  private
  save

!BOP

! !DEFINED PARAMETERS:
  integer, parameter :: CVMIX_KPP_INTERP_LMD94       = -1
  integer, parameter :: CVMIX_KPP_MATCH_BOTH         = 1
  integer, parameter :: CVMIX_KPP_MATCH_GRADIENT     = 2
  integer, parameter :: CVMIX_KPP_SIMPLE_SHAPES      = 3
  integer, parameter :: CVMIX_KPP_PARABOLIC_NONLOCAL = 4
  integer, parameter :: NO_LANGMUIR_MIXING           = -1
  integer, parameter :: LANGMUIR_MIXING_LWF16        = 1
  integer, parameter :: LANGMUIR_MIXING_RWHGK16      = 2
  integer, parameter :: NO_LANGMUIR_ENTRAINMENT      = -1
  integer, parameter :: LANGMUIR_ENTRAINMENT_LWF16   = 1
  integer, parameter :: LANGMUIR_ENTRAINMENT_LF17    = 2
  integer, parameter :: LANGMUIR_ENTRAINMENT_RWHGK16 = 3

! !PUBLIC MEMBER FUNCTIONS:

  public :: cvmix_init_kpp
  ! Note: cvmix_kpp_compute_OBL_depth would be part of cvmix_coeffs_kpp but
  !       CVMix can not smooth the boundary layer depth or correct the
  !       buoyancy flux term
  ! These are public for testing, may end up private later
  public :: EFactor_model
  public :: ustokes_SL_model

! !PUBLIC TYPES:

  ! cvmix_kpp_params_type contains the necessary parameters for KPP mixing
  type, public :: cvmix_kpp_params_type
    private
      real(cvmix_r8) :: Ri_crit        ! Critical Richardson number
                                       ! (OBL_depth = where bulk Ri = Ri_crit)

      real(cvmix_r8) :: minOBLdepth    ! Minimum allowable OBL depth
                                       ! (Default is 0 m => no minimum)
      real(cvmix_r8) :: maxOBLdepth    ! Maximum allowable OBL depth
                                       ! (Default is 0 m => no maximum)
      real(cvmix_r8) :: minVtsqr       ! Minimum allowable unresolved shear
                                       ! (Default is 1e-10 m^2/s^2)

      real(cvmix_r8) :: vonkarman      ! von Karman constant

      real(cvmix_r8) :: Cstar          ! coefficient for nonlinear transport
      real(cvmix_r8) :: nonlocal_coeff ! Cs from Eq (20) in LMD94
                                       ! Default value comes from paper, but
                                       ! some users may set it = 1.

      ! For velocity scale function, _m => momentum and _s => scalar (tracer)
      real(cvmix_r8) :: zeta_m         ! parameter for computing vel scale func
      real(cvmix_r8) :: zeta_s         ! parameter for computing vel scale func
      real(cvmix_r8) :: a_m            ! parameter for computing vel scale func
      real(cvmix_r8) :: c_m            ! parameter for computing vel scale func
      real(cvmix_r8) :: a_s            ! parameter for computing vel scale func
      real(cvmix_r8) :: c_s            ! parameter for computing vel scale func

      real(cvmix_r8) :: surf_layer_ext ! nondimensional extent of surface layer
                                       ! (expressed in sigma-coordinates)

      integer        :: interp_type    ! interpolation type used to interpolate
                                       ! bulk Richardson number
      integer        :: interp_type2   ! interpolation type used to interpolate
                                       ! diff and visc at OBL_depth

      ! Cv is a parameter used to compute the unresolved shear. By default, the
      ! formula from Eq. (A3) of Danabasoglu et al. is used, but a single
      ! scalar value can be set instead.
      real(cvmix_r8) :: Cv

      ! MatchTechnique is set by a string of the same name as an argument in
      ! cvmix_init_kpp. It determines how matching between the boundary layer
      ! and ocean interior is handled at the interface. Note that this also
      ! controls whether the shape function used to compute the coefficient in
      ! front of the nonlocal term is the same as that used to compute the
      ! gradient term.
      ! Options (for cvmix_init_kpp) are
      ! (i) SimpleShapes => Shape functions for both the gradient and nonlocal
      !                     terms vanish at interface
      ! (ii) MatchGradient => Shape function for nonlocal term vanishes at
      !                       interface, but gradient term matches interior
      !                       values.
      ! (iii) MatchBoth => Shape functions for both the gradient and nonlocal
      !                    term match interior values at interface
      ! (iv) ParabolicNonLocal => Shape function for the nonlocal term is
      !                         (1-sigma)^2, gradient term is sigma*(1-sigma)^2
      integer :: MatchTechnique

      ! Flag for what to do with old values of CVmix_vars%[MTS]diff
      integer :: handle_old_vals

      ! Logic flags to dictate if / how various terms are computed
      logical        :: lscalar_Cv     ! True => use the scalar Cv value
      logical        :: lEkman         ! True => compute Ekman depth limit
      logical        :: lMonOb         ! True => compute Monin-Obukhov limit
      logical        :: lnoDGat1       ! True => G'(1) = 0 (shape function)
                                       ! False => compute G'(1) as in LMD94
      logical        :: lenhanced_diff ! True => enhance diffusivity at OBL
      integer        :: Langmuir_Mixing_Opt
                                       ! Option of Langmuir enhanced mixing
                                       ! - apply an enhancement factor to the
                                       ! turbulent velocity scale
      integer        :: Langmuir_Entrainment_Opt
                                       ! Option of Langmuir turbulence enhanced
                                       ! entrainment - modify the unresolved shear
      logical        :: l_LMD_ws       ! flag to use original Large et al. (1994)
                                       ! equations for computing turbulent scales
                                       ! rather than the updated methodology in
                                       ! Danabasoglu et al. (2006). The latter
                                       ! limits sigma to be < surf_layer_extent
                                       ! when computing turbulent scales while
                                       ! the former only imposes this restriction
                                       ! in unstable regimes.
      real(cvmix_r8) :: c_LT, c_ST, c_CT  ! Empirical constants in the scaling of the
                                          ! entrainment buoyancy flux
                                          ! (20) in Li and Fox-Kemper, 2017, JPO
      real(cvmix_r8) :: p_LT              ! Power of Langmuir number in the above
                                          ! scaling
      !BGR
      real(cvmix_r8) :: RWHGK_ENTR_COEF,& ! Coefficient and exponent from
                        RWHGK_ENTR_EXP    ! RWHGK16 Langmuir parameterization

  end type cvmix_kpp_params_type

!EOP

type(cvmix_kpp_params_type), target :: CVmix_kpp_params_saved

contains

!EOC

  function EFactor_model(u10, ustar, hbl, CVmix_params_in)

! This function returns the enhancement factor, given the 10-meter
! wind (m/s), friction velocity (m/s) and the boundary layer depth (m).
!
! Qing Li, 160606

! Input
    real(cvmix_r8), intent(in) :: &
        ! 10 meter wind (m/s)
        u10, &
        ! water-side surface friction velocity (m/s)
        ustar, &
        ! boundary layer depth (m)
        hbl
    type(cvmix_global_params_type), intent(in) :: CVmix_params_in

! Local variables
    real(cvmix_r8) :: us_sl, lasl_sqr_i
    real(cvmix_r8) :: EFactor_model

    if (u10 .gt. cvmix_zero .and. ustar .gt. cvmix_zero) then
      ! surface layer averaged Stokes drift
      us_sl = ustokes_SL_model(u10, hbl, CVmix_params_in)
      !
      ! LaSL^{-2}
      lasl_sqr_i = us_sl/ustar
      !
      ! enhancement factor (Li et al., 2016)
      EFactor_model = sqrt(cvmix_one &
                 +cvmix_one/1.5_cvmix_r8**2*lasl_sqr_i &
                 +cvmix_one/5.4_cvmix_r8**4*lasl_sqr_i**2)
    else
      ! otherwise set to one
      EFactor_model = cvmix_one
    endif

  end function EFactor_model

  function ustokes_SL_model(u10, hbl, CVmix_params_in)

! This function returns the surface layer averaged Stokes drift, given
! the 10-meter wind (m/s) and the boundary layer depth (m).
!
! Qing Li, 20180130

! Input
    real(cvmix_r8), intent(in) :: &
        ! 10 meter wind (m/s)
        u10, &
        ! boundary layer depth (m)
        hbl
    type(cvmix_global_params_type), intent(in) :: CVmix_params_in
! Local variables
    ! parameters
    real(cvmix_r8), parameter :: &
        ! ratio of U19.5 to U10 (Holthuijsen, 2007)
        u19p5_to_u10 = 1.075_cvmix_r8, &
        ! ratio of mean frequency to peak frequency for
        ! Pierson-Moskowitz spectrum (Webb, 2011)
        fm_to_fp = 1.296_cvmix_r8, &
        ! ratio of surface Stokes drift to U10
        us_to_u10 = 0.0162_cvmix_r8, &
        ! loss ratio of Stokes transport
        r_loss = 0.667_cvmix_r8

    real(cvmix_r8) :: us, hm0, fm, fp, vstokes, kphil, kstar
    real(cvmix_r8) :: z0, z0i, r1, r2, r3, r4, tmp
    real(cvmix_r8) :: ustokes_SL_model

    if (u10 .gt. cvmix_zero) then
      ! surface Stokes drift
      us = us_to_u10*u10
      !
      ! significant wave height from Pierson-Moskowitz
      ! spectrum (Bouws, 1998)
      hm0 = 0.0246_cvmix_r8*u10**2
      !
      ! peak frequency (PM, Bouws, 1998)
      tmp = 2.0_cvmix_r8*cvmix_PI*u19p5_to_u10*u10
      fp = 0.877_cvmix_r8*CVmix_params_in%Gravity/tmp
      !
      ! mean frequency
      fm = fm_to_fp*fp
      !
      ! total Stokes transport (a factor r_loss is applied to account
      !  for the effect of directional spreading, multidirectional waves
      !  and the use of PM peak frequency and PM significant wave height
      !  on estimating the Stokes transport)
      vstokes = 0.125_cvmix_r8*cvmix_PI*r_loss*fm*hm0**2
      !
      ! the general peak wavenumber for Phillips' spectrum
      ! (Breivik et al., 2016) with correction of directional spreading
      kphil = 0.176_cvmix_r8*us/vstokes
      !
      ! surface layer averaged Stokes dirft with Stokes drift profile
      ! estimated from Phillips' spectrum (Breivik et al., 2016)
      ! the directional spreading effect from Webb and Fox-Kemper, 2015
      ! is also included
      kstar = kphil*2.56_cvmix_r8
      ! surface layer
      z0 = 0.2_cvmix_r8*abs(hbl)
      z0i = cvmix_one/z0
      ! term 1 to 4
      r1 = (0.151_cvmix_r8/kphil*z0i-0.84_cvmix_r8) &
            *(cvmix_one-exp(-2.0_cvmix_r8*kphil*z0))
      r2 = -(0.84_cvmix_r8+0.0591_cvmix_r8/kphil*z0i) &
             *sqrt(2.0_cvmix_r8*cvmix_PI*kphil*z0) &
             *erfc(sqrt(2.0_cvmix_r8*kphil*z0))
      r3 = (0.0632_cvmix_r8/kstar*z0i+0.125_cvmix_r8) &
            *(cvmix_one-exp(-2.0_cvmix_r8*kstar*z0))
      r4 = (0.125_cvmix_r8+0.0946_cvmix_r8/kstar*z0i) &
             *sqrt(2.0_cvmix_r8*cvmix_PI*kstar*z0) &
             *erfc(sqrt(2.0_cvmix_r8*kstar*z0))
      ustokes_SL_model = us*(0.715_cvmix_r8+r1+r2+r3+r4)
    else
      ustokes_SL_model = cvmix_zero
    endif

    end function ustokes_SL_model

end module theorywaves
