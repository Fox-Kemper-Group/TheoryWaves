#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SWLDMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  +-----------------------------------+
!/
!/    21-Nov-2011 : Origination.                        ( version 4.07 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Source term module for swell dissipation based on different
!     physics that can be independently selected form the input
!     and whitecapping dissipation terms in the model setup.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SWL4    Subr. Public   Ardhuin et al (2010+) swell dissipation
!      W3SWL6    Subr. Public   Babanin (2011) swell dissipation 
!
!      IRANGE    Func. Private  Generate a sequence of integer values
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!
!  6. Switches :
!
!     !/S  Enable subroutine tracing.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC  :: W3SWL4, W3SWL6
      PRIVATE :: IRANGE
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SWL4 (A, CG, WN, DAIR, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Aug-2021 |
!/                  +-----------------------------------+
!/
!/    29-May-2009 : Origination (w3srcxmd.ftn)          ( version 3.14 )
!/    06-Jan-2012 : Implementation                         (S. Zieger)
!/    13-Aug-2021 : Consider DAIR a variable           ( version x.xx )
!/
!  1. Purpose :
!
!     FIXME 
!
!  2. Method :
!
!  3. Parameters :
!
!      Parameter list
!     ----------------------------------------------------------------
!      A¹      R.A. I  Action density spectrum
!      CG      R.A. I  Group velocities
!      WN      R.A. I  Wavenumbers
!      DAIR    R.A. I   Air density
!      S¹      R.A. O  Source term 
!      D¹      R.A. O  Diagonal term of derivative
!      ¹ Stored as 1-D array with dimension NTH*NK (column by column).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      IRANGE    Func. W3SWLDMD 
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SRCE    Subr. W3SRCEMD Source term integration.
!      W3EXPO    Subr.   N/A    Point output post-processor.
!      GXEXPO    Subr.   N/A    GrADS point output post-processor.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See comments in source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: GRAV, DWAT
      USE W3GDATMD,  ONLY: NK, NTH, NSPEC, SIG2, DDEN, FTE, SWL6B1
#ifdef W3_S
      USE W3SERVMD, ONLY: STRACE
#endif
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
      REAL, INTENT(IN)  :: A(NSPEC), CG(NK), WN(NK), DAIR
      REAL, INTENT(OUT) :: S(NSPEC), D(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
#ifdef W3_S
      INTEGER, SAVE     :: IENT = 0
#endif
      INTEGER           :: IKN(NK), ITH
      REAL, PARAMETER   :: VA = 1.4E-5 ! Air kinematic viscosity (used in WAM).
      REAL              :: EB(NK), WN2(NSPEC), EMEAN
      REAL              :: FE, AORB, RE, RECRIT, UOSIG, CDSV
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'W3SWL4')
#endif
!
      IKN = IRANGE(1,NSPEC,NTH) 
      D   = 0.
      WN2 = 0.
!
      DO ITH = 1, NTH
         WN2(IKN+(ITH-1)) = WN        ! Wavenumbers to all directions.
      END DO
!
      EB     = SUM(RESHAPE(A,(/ NTH,NK /)),1) * DDEN(1:NK) / CG
      EMEAN  = SUM(EB) + (EB(NK) / DDEN(NK)) * FTE
!
      AORB   = 2.0*SQRT(EMEAN)
!     
      EB     = SUM(RESHAPE(A*SIG2**2,(/ NTH,NK /)),1) * DDEN(1:NK) / CG
      UOSIG  = 2.0*SQRT(SUM(EB))

      FE     = SWL6B1    ! (from NAMELIST)
!     FE     = 0.001     ! (from NAMELIST)
!/             0.001 - 0.019 with median value 0.007 (Ardhuin et al 2009, Babanin 2011)
      CDSV   = 1.2
!
      RECRIT = 1.0E5
      RE     = 4.0 * UOSIG * AORB / VA
!
      IF (RE .GT. RECRIT) THEN
         D = -(16.0/GRAV) * (DAIR/DWAT) * FE * (SIG2**2) *UOSIG
      ELSE
         D = -2.0 * (DAIR/DWAT) * CDSV * WN2 * SQRT(2.0 * VA * SIG2)
      END IF
!
      S = D * A
!
!  WRITE(*,*) ' FE       =',FE
!  WRITE(*,*) ' HS       =',4.*SQRT(EMEAN)
!  WRITE(*,*) ' UOSIG    =',UOSIG
!  WRITE(*,*) ' AORB     =',AORB
!  WRITE(*,*) ' RE/RECRIT=',RE/RECRIT
!  WRITE(*,*) ' SWL4_tot =',SUM(SUM(RESHAPE(S,(/ NTH,NK /)),1)*DDEN/CG)
!/
!/ End of W3SWL4 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SWL4
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SWL6 (A, CG, WN, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         16-Feb-2012 |
!/                  +-----------------------------------+
!/
!/    29-May-2009 : Origination (w3srcxmd.ftn)          ( version 3.14 )
!/    16-Feb-2012 : Implementation                      ( version 4.07 )
!/                                                         (S. Zieger)
!/
!  1. Purpose :
!
!     Turbulent dissipation of narrow-banded swell as described in
!     Babanin (2011, Section 7.5). 
!
!     Babanin 2011: Cambridge Press, 295-321, 463pp.
!
!  2. Method :
!
!     S = D * A
!
!  3. Parameters :
!
!      Parameter list
!     ----------------------------------------------------------------
!      A¹      R.A. I  Action density spectrum
!      CG      R.A. I  Group velocities
!      WN      R.A. I  Wavenumbers
!      S¹      R.A. O  Source term 
!      D¹      R.A. O  Diagonal term of derivative
!      ¹ Stored as 1-D array with dimension NTH*NK (column by column).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      IRANGE    Func. W3SWLDMD 
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SRCE    Subr. W3SRCEMD Source term integration.
!      W3EXPO    Subr.   N/A    Point output post-processor.
!      GXEXPO    Subr.   N/A    GrADS point output post-processor.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See comments in source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS, ONLY: GRAV
      USE W3GDATMD,  ONLY: NK, NTH, NSPEC, SIG, DDEN, DTH
      USE W3GDATMD,  ONLY: SWL6CSTB1, SWL6B1, FTE, FTWN
#ifdef W3_T6
      USE W3ODATMD,  ONLY: NDST
#endif
#ifdef W3_S
      USE W3SERVMD, ONLY: STRACE
#endif
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
      REAL, INTENT(IN)    :: A(NSPEC), CG(NK), WN(NK)
      REAL, INTENT(OUT)   :: S(NSPEC), D(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
#ifdef W3_S
      INTEGER, SAVE     :: IENT = 0
#endif
      INTEGER             :: IK, ITH, IKN(NK)
      REAL, DIMENSION(NK) :: ABAND, KMAX, ANAR, BN, AORB, DDIS
      REAL                :: K(NTH,NK), B1
!/
!/ ------------------------------------------------------------------- /
!/
#ifdef W3_S
      CALL STRACE (IENT, 'W3SWL6')
#endif
!
!/ 0) --- Initialize parameters -------------------------------------- /
      IKN   = IRANGE(1,NSPEC,NTH)            ! Index vector for array access, e.g.  
!                                            ! in form of WN(1:NK) == WN2(IKN).
      ABAND = SUM(RESHAPE(A,(/ NTH,NK /)),1) ! action density as function of wavenumber
      DDIS  = 0.
      D     = 0.
      B1    = SWL6B1                         ! empirical constant from NAMELIST
!
!/ 1) --- Choose calculation of steepness a*k ------------------------ /
!/        Replace the measure of steepness with the spectral
!         saturation after Banner et al. (2002) ---------------------- /
      K     = RESHAPE(A,(/ NTH,NK /))
      KMAX  = MAXVAL(K,1)
      DO IK = 1,NK
         IF (KMAX(IK).LT.1.0E-34) THEN
            K(1:NTH,IK) = 1.
         ELSE
            K(1:NTH,IK) = K(1:NTH,IK)/KMAX(IK)
         END IF
      END DO
      ANAR  = 1.0/( SUM(K,1) * DTH )
      BN    = ANAR * ( ABAND * SIG(1:NK) * DTH ) * WN**3
!
      IF (.NOT.SWL6CSTB1) THEN
!
!/    --- A constant value for B1 attenuates swell too strong in the
!/        western central Pacific (i.e. cross swell less than 1.0m).
!/        Workaround is to scale B1 with steepness a*kp, where kp is
!/        the peak wavenumber. SWL6B1 remains a scaling constant, but
!/        with different magnitude.  --------------------------------- /
          IK    = MAXLOC(ABAND,1)         ! Index for peak
!         EMEAN = SUM(ABAND * DDEN / CG)  ! Total sea surface variance
          B1    = SWL6B1 * ( 2. * SQRT(SUM(ABAND*DDEN/CG)) * WN(IK) )
!
      END IF
!
!/ 2) --- Calculate the derivative term only (in units of 1/s) ------- /
      DO IK = 1,NK
         IF (ABAND(IK) .GT. 1.E-30) THEN
            DDIS(IK) = -(2./3.) * B1 * SIG(IK) * SQRT(BN(IK))
         END IF
      END DO
!
!/ 3) --- Apply dissipation term of derivative to all directions ----- /
      DO ITH = 1, NTH
         D(IKN+(ITH-1)) = DDIS 
      END DO
!
      S = D * A
!
!  WRITE(*,*) ' B1       =',B1
!  WRITE(*,*) ' DDIS_tot =',SUM(DDIS*ABAND*DDEN/CG)
!  WRITE(*,*) ' EDENS_tot=',sum(aband*dden/cg)
!  WRITE(*,*) ' EDENS_tot=',sum(aband*sig*dth*dsii/cg)
!  WRITE(*,*) ' '
!  WRITE(*,*) ' SWL6_tot =',sum(SUM(RESHAPE(S,(/ NTH,NK /)),1)*DDEN/CG)
!
!/
!/ End of W3SWL6 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SWL6
!/ ------------------------------------------------------------------- /
!/
      FUNCTION IRANGE(X0,X1,DX) RESULT(IX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |           S. Zieger               |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Feb-2011 |
!/                  +-----------------------------------+
!/
!/    15-Feb-2011 : Origination from W3SRC6MD          ( version 4.07 )
!/                                                        (S. Zieger)
!/
!  1. Purpose :
!         Generate a linear-spaced sequence of integer
!         numbers. Used for array addressing (indexing).
!
!/
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: X0, X1, DX
      INTEGER, ALLOCATABLE :: IX(:)
      INTEGER              :: N
      INTEGER              :: I
!
      N = INT(REAL(X1-X0)/REAL(DX))+1
      ALLOCATE(IX(N))
      DO I = 1, N
         IX(I) = X0+ (I-1)*DX
      END DO
!/
      END FUNCTION IRANGE
!/ ------------------------------------------------------------------- /
!/
      END MODULE W3SWLDMD
