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
!  uses no other modules
!EOP

  implicit none
  private
  save

!BOP

! !DEFINED PARAMETERS:

  ! Kind Types:
  ! Use double precision for floating point computations.
  integer, parameter, public :: cvmix_r8       = selected_real_kind(15, 307)

  ! Global parameters:
  ! The constant 1 is used repeatedly. 
  ! The value for pi is needed.
  real(cvmix_r8), parameter, public :: cvmix_zero = real(0,cvmix_r8),         &
                                       cvmix_one  = real(1,cvmix_r8)
  real(cvmix_r8), parameter, public :: cvmix_PI   = &
                                       3.14159265358979323846_cvmix_r8

! !PUBLIC MEMBER FUNCTIONS:

  public :: EFactor_model
  public :: ustokes_SL_model

! !PUBLIC TYPES:

  ! cvmix_global_params_type contains global parameters used by multiple
  ! mixing methods.
  type, public :: cvmix_global_params_type
    ! maximum number of levels for any column
    integer :: max_nlev
             ! units: unitless

    real(cvmix_r8) :: Gravity = 9.80616_cvmix_r8

    ! Prandtl number
    real(cvmix_r8) :: prandtl
                    ! units: unitless

    ! Fresh water and salt water densities
    real(cvmix_r8) :: FreshWaterDensity
    real(cvmix_r8) :: SaltWaterDensity
                    ! units: kg m^-3

  end type cvmix_global_params_type

!EOP

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
