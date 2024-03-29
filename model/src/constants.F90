#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE CONSTANTS
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Jun-2018 |
!/                  +-----------------------------------+
!/
!/    11-Nov-1999 : Fortran 90 version.                 ( version 2.00 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    25-Jun-2011 : Adding Kelvin functions.            ( version 4.05 )
!/    03-Sep-2012 : Adding TSTOUT flag.                 ( version 4.10 )
!/    28-Feb-2013 : Adding cap at 0.5 in FWTABLE        ( version 4.08 )
!/    20-Jan-2017 : Add parameters for ESMF             ( version 6.02 )
!/    01-Mar-2018 : Add UNDEF parameter                 ( version 6.02 )
!/    05-Jun-2018 : Add PDLIB parameters                ( version 6.04 )
!/
!/    Copyright 2009-2012 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Define some much-used constants for global use (all defined
!     as PARAMETER).
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      GRAV      Real  Global   Acc. of gravity                 (m/s2)
!      DWAT      Real  Global   Density of water               (kg/m3)
!      DAIR      Real  Global   Density of air                 (kg/m3)
!      NU_AIR    Real  Global   Kinematic viscosity of air      (m2/s)
!      NU_WATER  Real  Global   Kinematic viscosity of water    (m2/s)
!      SED_SG    Real  Global   Specific gravity of sediments   (N.D.)
!      KAPPA     Real  Global   von Karman's constant           (N.D.)
!      PI        Real  Global   pi.
!      TPI       Real  Global   2pi.
!      HPI       Real  Global   0.5pi.
!      TPIINV    Real  Global   1/2pi.
!      HPIINV    Real  Global   2/pi.
!      RADE      Real  Global   Conv. factor from radians to degrees.
!      DERA      Real  Global   Conv. factor from degrees to radians.
!      RADIUS    Real  Global   Radius of the earth.             (m)
!      TSTOUT    Log.  Global   Flag for generation of test files.
!      UNDEF     Real  Global   Value for undefined variable in output
!     ----------------------------------------------------------------
!
!  5. Remarks
!
!      - The flag for generating test output files is included here as
!        it is needed in both ww3_shel and ww3_multi at the same time.
!        Make sure that this flag is true if you want to write to the 
!        test output file !
!
!/ ------------------------------------------------------------------- /
!/
      LOGICAL, PARAMETER      :: TSTOUT = .FALSE.
!
      REAL, PARAMETER         :: GRAV   =    9.806
      REAL, PARAMETER         :: DWAT   = 1000.
      REAL, PARAMETER         :: DAIR   =    1.225
      REAL, PARAMETER         :: nu_air  = 1.4E-5        
!mdo  *** Changing nu_water to be consistent with DWAT=1000 (assumes 10degC)
      REAL, PARAMETER         :: nu_water  = 1.31E-6    !mdo   WAS: 3.E-6        
      REAL, PARAMETER         :: sed_sg  = 2.65        
      REAL, PARAMETER         :: KAPPA = 0.40       !Von Karman's constant
!
      REAL, PARAMETER         :: PI     = 3.141592653589793
      REAL, PARAMETER         :: TPI    = 2.0 * PI
      REAL, PARAMETER         :: HPI    = 0.5 * PI
      REAL, PARAMETER         :: TPIINV = 1. / TPI
      REAL, PARAMETER         :: HPIINV = 1. / HPI
      REAL, PARAMETER         :: RADE   = 180. / PI
      REAL, PARAMETER         :: DERA   = PI / 180.
!
      REAL, PARAMETER         :: RADIUS = 4.E7 * TPIINV
!
      REAL, PARAMETER         :: G2PI3I = 1. / ( GRAV**2 * TPI**3 )
      REAL, PARAMETER         :: G1PI1I = 1. / ( GRAV * TPI )
!
      REAL                    :: UNDEF = -999.9
!PSH TheoryWaves begin
      REAL, PARAMETER         :: ZERO = 0.
      REAL, PARAMETER         :: ONE = 1. 
!PSH TheoryWaves end

!
! Parameters for friction factor table 
!
      INTEGER, PARAMETER       :: SIZEFWTABLE=300  
      REAL                     :: FWTABLE(0:SIZEFWTABLE)
      REAL                     :: DELAB
      REAL,    PARAMETER       :: ABMIN = -1. 
      REAL, PRIVATE, PARAMETER :: ABMAX = 8.
      INTEGER, PARAMETER       :: srce_direct = 0
      INTEGER, PARAMETER       :: srce_imp_post = 1
      INTEGER, PARAMETER       :: srce_imp_pre = 2
      INTEGER, PARAMETER       :: DEBUG_NODE = 1014
      INTEGER, PARAMETER       :: DEBUG_ELEMENT = 50
      LOGICAL                  :: LPDLIB = .FALSE. 
      LOGICAL                  :: LSETUP
!
! Parameters in support of running as ESMF component
!
! --- Flag indicating whether or not the model has been invoked as an
!     ESMF Component.  This flag is set to true in the WMESMFMD ESMF
!     module during initialization.
      LOGICAL :: IS_ESMF_COMPONENT = .FALSE.
!
      CONTAINS
! ----------------------------------------------------------------------
      SUBROUTINE TABU_FW
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Feb-2013 |
!/                  +-----------------------------------+
!/
!/    19-Oct-2007 : Origination.                        ( version 3.13 )
!/    28-Feb-2013 : Caps the friction factor to 0.5     ( version 4.08 )
!/
!  1. Purpose :
!     TO estimate friction coefficients in oscillatory boundary layers
!     METHOD.
!      tabulation on Kelvin functions
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WW3_GRID  Prog. WW3_GRID Model grid initialization
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
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
      INTEGER, PARAMETER      :: NITER=100
      REAL   , PARAMETER      :: XM=0.50, EPS1=0.00001
!     VARIABLE.   TYPE.     PURPOSE.
!      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *XNU*       REAL      KINEMATIC VISCOSITY OF AIR.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
!      *EPS1*      REAL      SMALL NUMBER TO MAKE SURE THAT A SOLUTION
!                            IS OBTAINED IN ITERATION WITH TAU>TAUW.
! ----------------------------------------------------------------------
      INTEGER I,ITER
      REAL KER, KEI
      REAL ABR,ABRLOG,L10,FACT,FSUBW,FSUBWMEMO,dzeta0,dzeta0memo
!
!
!
      DELAB   = (ABMAX-ABMIN)/REAL(SIZEFWTABLE)
      L10=ALOG(10.)
      DO I=0,SIZEFWTABLE
!
!  index I in this table corresponds to a normalized roughness z0/ABR = 10^ABMIN+REAL(I)*DELAB
!
         ABRLOG=ABMIN+REAL(I)*DELAB
         ABR=EXP(ABRLOG*L10)
         FACT=1/ABR/(21.2*KAPPA)
         FSUBW=0.05
         dzeta0=0.
         DO ITER=1,NITER
            fsubwmemo=fsubw
            dzeta0memo=dzeta0
            dzeta0=fact*fsubw**(-.5)
            CALL KERKEI(2.*SQRT(dzeta0),ker,kei)
            fsubw=.08/(ker**2+kei**2)
            fsubw=.5*(fsubwmemo+fsubw)
            dzeta0=.5*(dzeta0memo+dzeta0)
            END DO   
!
! Maximum value of 0.5 for fe is based on field 
! and lab experiment by Lowe et al. JGR 2005, 2007 
! 
            FWTABLE(I)  = MIN(fsubw,0.5) 
!           WRITE(994,*) 'Friction factor:',I,ABR,FWTABLE(I)
         END DO
      RETURN
      END SUBROUTINE TABU_FW

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KZEONE(X, Y, RE0, IM0, RE1, IM1)                      
!  June 1999 adaptation to CRESTb, all tests on range of (x,y) have been
!  bypassed, we implicitly expect X to be positive or |x,y| non zero
! 
! This subroutine is copyright by ACM
! see http://www.acm.org/pubs/copyright_policy/softwareCRnotice.html
! ACM declines any responsibility of any kind
! 
! THE VARIABLES X AND Y ARE THE REAL AND IMAGINARY PARTS OF
! THE ARGUMENT OF THE FIRST TWO MODIFIED BESSEL FUNCTIONS
! OF THE SECOND KIND,K0 AND K1.  RE0,IM0,RE1 AND IM1 GIVE
! THE REAL AND IMAGINARY PARTS OF EXP(X)*K0 AND EXP(X)*K1,
! RESPECTIVELY.  ALTHOUGH THE REAL NOTATION USED IN THIS
! SUBROUTINE MAY SEEM INELEGANT WHEN COMPARED WITH THE
! COMPLEX NOTATION THAT FORTRAN ALLOWS, THIS VERSION RUNS
! ABOUT 30 PERCENT FASTER THAN ONE WRITTEN USING COMPLEX
! VARIABLES.
! ACM Libraries
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   IMPLICIT NONE
   DOUBLE PRECISION X, Y, X2, Y2, RE0, IM0, RE1, IM1, &
      R1, R2, T1, T2, P1, P2, RTERM, ITERM, L
   DOUBLE PRECISION , PARAMETER, DIMENSION(8) :: EXSQ = &
         (/ 0.5641003087264D0,0.4120286874989D0,0.1584889157959D0, & 
            0.3078003387255D-1,0.2778068842913D-2,0.1000044412325D-3, &
            0.1059115547711D-5,0.1522475804254D-8 /)
   DOUBLE PRECISION , PARAMETER, DIMENSION(8) :: TSQ = &
         (/ 0.0D0,3.19303633920635D-1,1.29075862295915D0, &
            2.95837445869665D0,5.40903159724444D0,8.80407957805676D0, &
            1.34685357432515D1,2.02499163658709D1 /)
   INTEGER N,M,K
! THE ARRAYS TSQ AND EXSQ CONTAIN THE SQUARE OF THE
! ABSCISSAS AND THE WEIGHT FACTORS USED IN THE GAUSS-
! HERMITE QUADRATURE.
      R2 = X*X + Y*Y
      IF (R2.GE.1.96D2) GO TO 50
      IF (R2.GE.1.849D1) GO TO 30
! THIS SECTION CALCULATES THE FUNCTIONS USING THE SERIES
! EXPANSIONS
      X2 = X/2.0D0
      Y2 = Y/2.0D0
      P1 = X2*X2
      P2 = Y2*Y2
      T1 = -(DLOG(P1+P2)/2.0D0+0.5772156649015329D0)
! THE CONSTANT IN THE PRECEDING STATEMENT IS EULER*S
! CONSTANT
      T2 = -DATAN2(Y,X)
      X2 = P1 - P2
      Y2 = X*Y2
      RTERM = 1.0D0
      ITERM = 0.0D0
      RE0 = T1
      IM0 = T2
      T1 = T1 + 0.5D0
      RE1 = T1
      IM1 = T2
      P2 = DSQRT(R2)
      L = 2.106D0*P2 + 4.4D0
      IF (P2.LT.8.0D-1) L = 2.129D0*P2 + 4.0D0
      DO 20 N=1,INT(L)
        P1 = N
        P2 = N*N
        R1 = RTERM
        RTERM = (R1*X2-ITERM*Y2)/P2
        ITERM = (R1*Y2+ITERM*X2)/P2
        T1 = T1 + 0.5D0/P1
        RE0 = RE0 + T1*RTERM - T2*ITERM
        IM0 = IM0 + T1*ITERM + T2*RTERM
        P1 = P1 + 1.0D0
        T1 = T1 + 0.5D0/P1
        RE1 = RE1 + (T1*RTERM-T2*ITERM)/P1
        IM1 = IM1 + (T1*ITERM+T2*RTERM)/P1
   20 CONTINUE
      R1 = X/R2 - 0.5D0*(X*RE1-Y*IM1)
      R2 = -Y/R2 - 0.5D0*(X*IM1+Y*RE1)
      P1 = DEXP(X)
      RE0 = P1*RE0
      IM0 = P1*IM0
      RE1 = P1*R1
      IM1 = P1*R2
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE INTEGRAL
! REPRESENTATION, EQN 3, EVALUATED WITH 15 POINT GAUSS-
! HERMITE QUADRATURE
   30 X2 = 2.0D0*X
      Y2 = 2.0D0*Y
      R1 = Y2*Y2
      P1 = DSQRT(X2*X2+R1)
      P2 = DSQRT(P1+X2)
      T1 = EXSQ(1)/(2.0D0*P1)
      RE0 = T1*P2
      IM0 = T1/P2
      RE1 = 0.0D0
      IM1 = 0.0D0
      DO 40 N=2,8
        T2 = X2 + TSQ(N)
        P1 = DSQRT(T2*T2+R1)
        P2 = DSQRT(P1+T2)
        T1 = EXSQ(N)/P1
        RE0 = RE0 + T1*P2
        IM0 = IM0 + T1/P2
        T1 = EXSQ(N)*TSQ(N)
        RE1 = RE1 + T1*P2
        IM1 = IM1 + T1/P2
   40 CONTINUE
      T2 = -Y2*IM0
      RE1 = RE1/R2
      R2 = Y2*IM1/R2
      RTERM = 1.41421356237309D0*DCOS(Y)
      ITERM = -1.41421356237309D0*DSIN(Y)
! THE CONSTANT IN THE PREVIOUS STATEMENTS IS,OF COURSE,
! SQRT(2.0).
      IM0 = RE0*ITERM + T2*RTERM
      RE0 = RE0*RTERM - T2*ITERM
      T1 = RE1*RTERM - R2*ITERM
      T2 = RE1*ITERM + R2*RTERM
      RE1 = T1*X + T2*Y
      IM1 = -T1*Y + T2*X
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE
! ASYMPTOTIC EXPANSIONS
   50 RTERM = 1.0D0
      ITERM = 0.0D0
      RE0 = 1.0D0
      IM0 = 0.0D0
      RE1 = 1.0D0
      IM1 = 0.0D0
      P1 = 8.0D0*R2
      P2 = DSQRT(R2)
      L = 3.91D0+8.12D1/P2
      R1 = 1.0D0
      R2 = 1.0D0
      M = -8
      K = 3
      DO 60 N=1,INT(L)
        M = M + 8
        K = K - M
        R1 = FLOAT(K-4)*R1
        R2 = FLOAT(K)*R2
        T1 = FLOAT(N)*P1
        T2 = RTERM
        RTERM = (T2*X+ITERM*Y)/T1
        ITERM = (-T2*Y+ITERM*X)/T1
        RE0 = RE0 + R1*RTERM
        IM0 = IM0 + R1*ITERM
        RE1 = RE1 + R2*RTERM
        IM1 = IM1 + R2*ITERM
   60 CONTINUE
      T1 = DSQRT(P2+X)
      T2 = -Y/T1
      P1 = 8.86226925452758D-1/P2
! THIS CONSTANT IS SQRT(PI)/2.0, WITH PI=3.14159...
      RTERM = P1*DCOS(Y)
      ITERM = -P1*DSIN(Y)
      R1 = RE0*RTERM - IM0*ITERM
      R2 = RE0*ITERM + IM0*RTERM
      RE0 = T1*R1 - T2*R2
      IM0 = T1*R2 + T2*R1
      R1 = RE1*RTERM - IM1*ITERM
      R2 = RE1*ITERM + IM1*RTERM
      RE1 = T1*R1 - T2*R2
      IM1 = T1*R2 + T2*R1
      RETURN
      END SUBROUTINE KZEONE
      
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KERKEI(X,KER,KEI)
!**********************************************************************
! Computes the values of the zeroth order Kelvin function Ker and Kei
! These functions are used to determine the friction factor fw as a 
! function of the bottom roughness length assuming a linear profile
! of eddy viscosity (See Grant and Madsen, 1979)
!**********************************************************************
   IMPLICIT NONE
   
   DOUBLE PRECISION ZR,ZI,CYR,CYI,CYR1,CYI1
   REAL X,KER,KEI
   
   ZR=X*.50D0*SQRT(2.0D0)
   ZI=ZR
   CALL KZEONE(ZR, ZI, CYR, CYI,CYR1,CYI1)
   KER=CYR/EXP(ZR)
   KEI=CYI/EXP(ZR)
END SUBROUTINE KERKEI
!/
!/ End of module CONSTANTS ------------------------------------------- /
!/
      END MODULE CONSTANTS
