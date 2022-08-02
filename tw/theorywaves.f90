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
! CONSTANTS
!EOP
      USE CONSTANTS, ONLY: GRAV, PI, ZERO, ONE

  implicit none
!  private
!  save

!BOP

! !DEFINED PARAMETERS:

  ! Kind Types:
  ! Use double precision for floating point computations.
!  integer, parameter :: tw_r8       = selected_real_kind(15, 307)

  ! Global parameters:
  ! The constant 1 is used repeatedly. 
  ! The value for pi is needed.
!  real(tw_r8), parameter :: tw_zero = real(0,tw_r8),         &
!                            tw_one  = real(1,tw_r8)
!  real(tw_r8), parameter :: PI      = &
!                               3.14159265358979323846_tw_r8
!  real(tw_r8), parameter :: Gravity = &
!                               9.80616_tw_r8


! !PUBLIC MEMBER FUNCTIONS:

  public :: EFactor_model
  public :: ustokes_SL_model

! !PUBLIC TYPES:

!EOP

contains

!EOC

  function EFactor_model(u10, ustar, hbl)

! This function returns the enhancement factor, given the 10-meter
! wind (m/s), friction velocity (m/s) and the boundary layer depth (m).
!
! Qing Li, 160606

! Input
!    real(tw_r8), intent(in) :: &
!        ! 10 meter wind (m/s)
!        u10, &
!        ! water-side surface friction velocity (m/s)
!        ustar, &
!        ! boundary layer depth (m)
!        hbl
    REAL, intent(in) :: u10, ustar, hbl

! Local variables
    REAL :: us_sl, lasl_sqr_i
    REAL :: EFactor_model

    if (u10 .gt. ZERO .and. ustar .gt. ZERO) then
      ! surface layer averaged Stokes drift
      us_sl = ustokes_SL_model(u10, hbl)
      !
      ! LaSL^{-2}
      lasl_sqr_i = us_sl/ustar
      !
      ! enhancement factor (Li et al., 2016)
      EFactor_model = sqrt(ONE &
                 +ONE/1.5**2*lasl_sqr_i &
                 +ONE/5.4**4*lasl_sqr_i**2)
    else
      ! otherwise set to one
      EFactor_model = ONE
    endif

  end function EFactor_model

  function ustokes_SL_model(u10, hbl)

! This function returns the surface layer averaged Stokes drift, given
! the 10-meter wind (m/s) and the boundary layer depth (m).
!
! Qing Li, 20180130

! Input
    real(tw_r8), intent(in) :: &
        ! 10 meter wind (m/s)
        u10, &
        ! boundary layer depth (m)
        hbl
! Local variables
    ! parameters
    real(tw_r8), parameter :: &
        ! ratio of U19.5 to U10 (Holthuijsen, 2007)
        u19p5_to_u10 = 1.075, &
        ! ratio of mean frequency to peak frequency for
        ! Pierson-Moskowitz spectrum (Webb, 2011)
        fm_to_fp = 1.296, &
        ! ratio of surface Stokes drift to U10
        us_to_u10 = 0.0162, &
        ! loss ratio of Stokes transport
        r_loss = 0.667

    real(tw_r8) :: us, hm0, fm, fp, vstokes, kphil, kstar
    real(tw_r8) :: z0, z0i, r1, r2, r3, r4, tmp
    real(tw_r8) :: ustokes_SL_model

    if (u10 .gt. ZERO) then
      ! surface Stokes drift
      us = us_to_u10*u10
      !
      ! significant wave height from Pierson-Moskowitz
      ! spectrum (Bouws, 1998)
      hm0 = 0.0246*u10**2
      !
      ! peak frequency (PM, Bouws, 1998)
      tmp = 2.0*PI*u19p5_to_u10*u10
      fp = 0.877*GRAV/tmp
      !
      ! mean frequency
      fm = fm_to_fp*fp
      !
      ! total Stokes transport (a factor r_loss is applied to account
      !  for the effect of directional spreading, multidirectional waves
      !  and the use of PM peak frequency and PM significant wave height
      !  on estimating the Stokes transport)
      vstokes = 0.125*PI*r_loss*fm*hm0**2
      !
      ! the general peak wavenumber for Phillips' spectrum
      ! (Breivik et al., 2016) with correction of directional spreading
      kphil = 0.176*us/vstokes
      !
      ! surface layer averaged Stokes dirft with Stokes drift profile
      ! estimated from Phillips' spectrum (Breivik et al., 2016)
      ! the directional spreading effect from Webb and Fox-Kemper, 2015
      ! is also included
      kstar = kphil*2.56
      ! surface layer
      z0 = 0.2*abs(hbl)
      z0i = ONE/z0
      ! term 1 to 4
      r1 = (0.151/kphil*z0i-0.84) &
            *(ONE-exp(-2.0*kphil*z0))
      r2 = -(0.84+0.0591/kphil*z0i) &
             *sqrt(2.0*PI*kphil*z0) &
             *erfc(sqrt(2.0*kphil*z0))
      r3 = (0.0632/kstar*z0i+0.125) &
            *(ONE-exp(-2.0*kstar*z0))
      r4 = (0.125+0.0946/kstar*z0i) &
             *sqrt(2.0*PI*kstar*z0) &
             *erfc(sqrt(2.0*kstar*z0))
      ustokes_SL_model = us*(0.715+r1+r2+r3+r4)
    else
      ustokes_SL_model = ZERO
    endif

    end function ustokes_SL_model

end module theorywaves
