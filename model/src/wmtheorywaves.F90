!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 MODULE WMTHEORYWAVES

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
!  private
!  save

!BOP

! !DEFINED PARAMETERS:

  ! Kind Types:
  ! Use double precision for floating point computations.
  integer, parameter :: tw_r8       = selected_real_kind(15, 307)

  ! Global parameters:
  ! The constant 1 is used repeatedly. 
  ! The value for pi is needed.
  real(tw_r8), parameter :: tw_zero = real(0,tw_r8),         &
                            tw_one  = real(1,tw_r8)
  real(tw_r8), parameter :: PI      = &
                               3.14159265358979323846_tw_r8
  real(tw_r8), parameter :: Gravity = &
                               9.80616_tw_r8

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
    real(tw_r8), intent(in) :: &
        ! 10 meter wind (m/s)
        u10, &
        ! water-side surface friction velocity (m/s)
        ustar, &
        ! boundary layer depth (m)
        hbl

! Local variables
    real(tw_r8) :: us_sl, lasl_sqr_i
    real(tw_r8) :: EFactor_model

    if (u10 .gt. tw_zero .and. ustar .gt. tw_zero) then
      ! surface layer averaged Stokes drift
      us_sl = ustokes_SL_model(u10, hbl)
      !
      ! LaSL^{-2}
      lasl_sqr_i = us_sl/ustar
      !
      ! enhancement factor (Li et al., 2016)
      EFactor_model = sqrt(tw_one &
                 +tw_one/1.5_tw_r8**2*lasl_sqr_i &
                 +tw_one/5.4_tw_r8**4*lasl_sqr_i**2)
    else
      ! otherwise set to one
      EFactor_model = tw_one
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
        u19p5_to_u10 = 1.075_tw_r8, &
        ! ratio of mean frequency to peak frequency for
        ! Pierson-Moskowitz spectrum (Webb, 2011)
        fm_to_fp = 1.296_tw_r8, &
        ! ratio of surface Stokes drift to U10
        us_to_u10 = 0.0162_tw_r8, &
        ! loss ratio of Stokes transport
        r_loss = 0.667_tw_r8

    real(tw_r8) :: us, hm0, fm, fp, vstokes, kphil, kstar
    real(tw_r8) :: z0, z0i, r1, r2, r3, r4, tmp
    real(tw_r8) :: ustokes_SL_model

    if (u10 .gt. tw_zero) then
      ! surface Stokes drift
      us = us_to_u10*u10
      !
      ! significant wave height from Pierson-Moskowitz
      ! spectrum (Bouws, 1998)
      hm0 = 0.0246_tw_r8*u10**2
      !
      ! peak frequency (PM, Bouws, 1998)
      tmp = 2.0_tw_r8*PI*u19p5_to_u10*u10
      fp = 0.877_tw_r8*Gravity/tmp
      !
      ! mean frequency
      fm = fm_to_fp*fp
      !
      ! total Stokes transport (a factor r_loss is applied to account
      !  for the effect of directional spreading, multidirectional waves
      !  and the use of PM peak frequency and PM significant wave height
      !  on estimating the Stokes transport)
      vstokes = 0.125_tw_r8*PI*r_loss*fm*hm0**2
      !
      ! the general peak wavenumber for Phillips' spectrum
      ! (Breivik et al., 2016) with correction of directional spreading
      kphil = 0.176_tw_r8*us/vstokes
      !
      ! surface layer averaged Stokes dirft with Stokes drift profile
      ! estimated from Phillips' spectrum (Breivik et al., 2016)
      ! the directional spreading effect from Webb and Fox-Kemper, 2015
      ! is also included
      kstar = kphil*2.56_tw_r8
      ! surface layer
      z0 = 0.2_tw_r8*abs(hbl)
      z0i = tw_one/z0
      ! term 1 to 4
      r1 = (0.151_tw_r8/kphil*z0i-0.84_tw_r8) &
            *(tw_one-exp(-2.0_tw_r8*kphil*z0))
      r2 = -(0.84_tw_r8+0.0591_tw_r8/kphil*z0i) &
             *sqrt(2.0_tw_r8*PI*kphil*z0) &
             *erfc(sqrt(2.0_tw_r8*kphil*z0))
      r3 = (0.0632_tw_r8/kstar*z0i+0.125_tw_r8) &
            *(tw_one-exp(-2.0_tw_r8*kstar*z0))
      r4 = (0.125_tw_r8+0.0946_tw_r8/kstar*z0i) &
             *sqrt(2.0_tw_r8*PI*kstar*z0) &
             *erfc(sqrt(2.0_tw_r8*kstar*z0))
      ustokes_SL_model = us*(0.715_tw_r8+r1+r2+r3+r4)
    else
      ustokes_SL_model = tw_zero
    endif

    end function ustokes_SL_model

END MODULE WMTHEORYWAVES
