!                *********************************
                 SUBROUTINE SUSPENSION_COMPUTE_CAE
!                *********************************
!
     &(TAUP,HN,FDM,FD90,NPOIN,CHARR,XMVE,XMVS,VCE,GRAV,HMIN,XWC,
     & ZERO,ZREF,AC,CSTAEQ,QSC,ICQ,U2D,V2D,CSRATIO,T14,DEBUG)
!
!***********************************************************************
! SISYPHE   V7P3                                             28/03/2017
!***********************************************************************
!
!brief    COMPUTES EQUILIBRIUM CONCENTRATION;
!+
!
!history  R. WALTHER (ARTELIA), J. FONTAINE (EDF-LNHE)
!+        28/03/2017
!+        V7P3
!+  Creation of the subroutine.
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AC             |-->| CRITICAL SHIELDS PARAMETER
!| ACLADM         |-->| MEAN DIAMETER OF SEDIMENT
!| CHARR          |-->| LOGICAL, IF BEDLOAD OR NOT
!| CSTAEQ         |<->| EQUILIBRIUM CONCENTRATION
!| CSRATIO        |-->| EQUILIBRIUM CONCENTRATION FOR SOULBY-VAN RIJN EQ.
!| DEBUG          |-->| FLAG FOR DEBUGGING
!| FLUER          |<->| EROSION FLUX
!| GRAV           |-->| ACCELERATION OF GRAVITY
!| HMIN           |-->| MINIMUM VALUE OF WATER DEPTH
!| HN             |-->| WATER DEPTH
!| ICQ            |-->| REFERENCE CONCENTRATION FORMULA
!| NPOIN          |-->| NUMBER OF POINTS
!| QSC            |-->| BED LOAD TRANSPORT RATE
!| TAUP           |-->| SKIN FRICTION
!| VCE            |-->| FLOW VISCOSITY
!| XMVE           |-->| FLUID DENSITY
!| XMVS           |-->| SEDIMENT DENSITY
!| XWC            |-->| SETTLING VELOCITIES
!| ZERO           |-->| ZERO
!| ZREF           |-->| REFERENCE ELEVATION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE INTERFACE_SISYPHE
      USE BIEF
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE

      ! GLOBAL VARIABLES
      ! -------------------
      TYPE (BIEF_OBJ),  INTENT(IN)    :: TAUP,HN,ZREF,QSC
      TYPE (BIEF_OBJ),  INTENT(IN)    :: U2D,V2D,CSRATIO,T14
      INTEGER,          INTENT(IN)    :: NPOIN,DEBUG,ICQ
      LOGICAL,          INTENT(IN)    :: CHARR
      DOUBLE PRECISION, INTENT(IN)    :: XMVE,XMVS,GRAV,HMIN,XWC,VCE
      DOUBLE PRECISION, INTENT(IN)    :: ZERO,FDM,FD90
      DOUBLE PRECISION, INTENT(IN)    :: AC
      TYPE (BIEF_OBJ),  INTENT(INOUT) :: CSTAEQ

      ! LOCAL VARIABLES
      ! -------------------

      INTEGER I
!
!======================================================================!
!======================================================================!
!                               PROGRAM                                !
!======================================================================!
!======================================================================!
!
!  COMPUTES THE NEAR BED EQUILIBRIUM CONCENTRATION --> CSTAEQ (MEAN DIAMETER)
!
      IF(ICQ.EQ.1) THEN
!
        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_FREDSOE'
        CALL SUSPENSION_FREDSOE(FDM,TAUP,NPOIN,
     &                          GRAV,XMVE,XMVS,ZERO,AC,CSTAEQ)
        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_FREDSOE'
!
      ELSEIF(ICQ.EQ.2) THEN
        WRITE(LU,*)'BIJKER TO REPROGRAM IN CASE OF SUSPENSION'
!        TO REPROGRAM
!       BECAUSE MULTIPLICATION BY AVA TAKEN
!       INTO ACCOUNT AFTER THE BEDLOAD TRANSPORT RATE
!       SO PERHAPS RECALL BEDLOAD_BIJKER BEFORE
!        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_BIJKER'
!        CALL SUSPENSION_BIJKER(TAUP,HN,NPOIN,CHARR,QSC,ZREF,
!     &                         ZERO,HMIN,CSTAEQ,XMVE)
!        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_BIJKER'
!
      ELSEIF(ICQ.EQ.3) THEN
!
        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_VANRIJN'
        CALL SUSPENSION_VANRIJN(FDM,TAUP,NPOIN,GRAV,XMVE,XMVS,
     &                          VCE,ZERO,AC,CSTAEQ,ZREF)
        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_VANRIJN'

      ELSEIF(ICQ.EQ.4) THEN
        IF(DEBUG > 0) WRITE(LU,*) 'SUSPENSION_SANDFLOW'
        CALL SUSPENSION_SANDFLOW(FDM,FD90,TAUP,NPOIN,GRAV,XMVE,XMVS,
     &              ZERO,AC,CSTAEQ,ZREF,HN,U2D,V2D,CSRATIO)
        IF(DEBUG > 0) WRITE(LU,*) 'END SUSPENSION_SANDFLOW'
!
      ENDIF
!
!======================================================================!
!======================================================================!
!
      RETURN
      END
!
!#######################################################################
!
