!                    ************************
                     SUBROUTINE INIT_SEDIMENT
!                    ************************
!
     &(NSICLA,ELAY,ZF,ZR,NPOIN,AVAIL,FRACSED_GF,AVA0,
     & LGRAFED,CALWC,XMVS,XMVE,GRAV,VCE,XWC,FDM,
     & CALAC,AC,SEDCO,ES,ES_SABLE, ES_VASE ,NOMBLAY,CONC_VASE,
     & MS_SABLE,MS_VASE,ACLADM,UNLADM,TOCE_SABLE,
     & CONC,NLAYER,DEBU,MIXTE)
!
!***********************************************************************
! SISYPHE   V7P2                                   27/06/2016
!***********************************************************************
!
!brief
!
!history  C. VILLARET (LNHE)
!+        30/12/2008
!+
!+
!
!history  JMH
!+        16/09/2009
!+        V6P0
!+   AVAIL(NPOIN,10,NSICLA)
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  C.VILLARET (EDF-LNHE), P.TASSI (EDF-LNHE)
!+        19/07/2011
!+        V6P1
!+  Name of variables
!+
!
!history  P.TASSI (EDF-LNHE)
!+        30/05/2012
!+        V6P2
!+  Case DSTAR > 150 AC(I) = 0.045D0
!+
!
!history  P.TASSI (EDF-LNHE)
!+        06/07/2012
!+        V6P2
!+  Line MIXTE=.FALSE. added.
!
!history  R.KOPMANN (BAW)
!+        27/06/2016
!+        V7P2
!+  CONTINUOUS APPROSIMATION CRITICAL SHILEDS CURVE DSTAR > 72 AC(I) = 0.045D0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AC             |<->| CRITICAL SHIELDS PARAMETER
!| ACLADM         |-->| MEAN DIAMETER OF SEDIMENT
!| AT0            |<->| TIME IN S
!| AVAIL          |<->| VOLUME PERCENT OF EACH CLASS
!| CALAC          |-->| IF YES, SHIELDS PARAMETER FOUND IN PARAMETER FILE
!| CALWC          |-->| IF YES, SETTLING VELOCITIES FOUND IN PARAMETER FILE
!| CONC_VASE      |<->| MUD CONCENTRATION FOR EACH LAYER
!| ELAY           |<->| THICKNESS OF SURFACE LAYER
!| ES             |<->| LAYER THICKNESSES AS DOUBLE PRECISION
!| ES_SABLE       |<->| LAYER THICKNESSES OF SAND AS DOUBLE PRECISION
!| ES_VASE        |<->| LAYER THICKNESSES OF MUD AS DOUBLE PRECISION
!| FDM            |-->| DIAMETER DM FOR EACH CLASS
!| FRACSED_GF     |-->|(A SUPPRIMER)
!| GRAV           |-->| ACCELERATION OF GRAVITY
!| LGRAFED        |-->|(A SUPPRIMER)
!| MS_SABLE       |<->| MASS OF SAND PER LAYER (KG/M2)
!| MS_VASE        |<->| MASS OF MUD PER LAYER (KG/M2)
!| ES_SABLE       |<->| THICKNESS OF SAND LAYER (M)
!| ES_VASE        |<->| THICKNESS OF MUD LAYER  (M)
!| MIXTE          |<->| SEDIMENT MIXTE  (SABLE /VASE)
!| NOMBLAY        |-->| NUMBER OF BED LAYERS
!| NPOIN          |-->| NUMBER OF POINTS
!| NSICLA         |-->| NUMBER OF SEDIMENT CLASSES
!| SEDCO          |-->| LOGICAL, SEDIMENT COHESIVE OR NOT
!| UNLADM         |-->| MEAN DIAMETER OF ACTIVE STRATUM LAYER
!| VCE            |-->| WATER VISCOSITY
!| XMVE           |-->| FLUID DENSITY
!| XMVS           |-->| WATER DENSITY
!| XWC            |-->| SETTLING VELOCITY
!| ZF             |-->| ELEVATION OF BOTTOM
!| ZR             |-->| NON ERODABLE BED
!| CONC           |<->| CONCENTRATION OF BED LAYER
!| NLAYER         |<->| NUMBER OF BED LAYER
!| DEBU           |-->| FLAG, RESTART ON SEDIMENTOLOGICAL FILE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE INTERFACE_SISYPHE, EX_INIT_SEDIMENT => INIT_SEDIMENT
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER,           INTENT(IN)     :: NSICLA,NPOIN,NOMBLAY
      TYPE(BIEF_OBJ),    INTENT(INOUT)  :: ELAY,ZF,ZR
      TYPE(BIEF_OBJ), INTENT(INOUT)     :: MS_SABLE, MS_VASE
      TYPE(BIEF_OBJ),    INTENT(INOUT)  :: ACLADM, UNLADM
      TYPE(BIEF_OBJ),    INTENT(INOUT)  :: NLAYER
      LOGICAL,           INTENT(IN)     :: LGRAFED,CALWC
      LOGICAL,           INTENT(IN)     :: CALAC
      DOUBLE PRECISION,  INTENT(IN)     :: XMVS,XMVE,GRAV,VCE
      DOUBLE PRECISION,  INTENT(INOUT)  :: AVA0(NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: AVAIL(NPOIN,NOMBLAY,NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: FRACSED_GF(NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: FDM(NSICLA),XWC(NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: AC(NSICLA),TOCE_SABLE
      LOGICAL,           INTENT(IN)     :: SEDCO(NSICLA), DEBU
      LOGICAL,           INTENT(IN)     :: MIXTE
      DOUBLE PRECISION, INTENT(IN)    :: CONC_VASE(NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT) :: ES(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT) :: ES_SABLE(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT) :: ES_VASE(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT) :: CONC(NPOIN,NOMBLAY)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER            :: I,J
      INTEGER            :: ILAYER,IMUD,ISAND
      DOUBLE PRECISION   :: DENS,DSTAR
      DOUBLE PRECISION   :: CHECK_RS,CHECK_RM
      DOUBLE PRECISION   :: RATIO_MUD_VOL !ratio of mud volume to volume of porosity of sand
      DOUBLE PRECISION   :: MASS_SAND_TOT,MASS_MUD_TOT,MASS_TOT
      LOGICAL            :: BED_MIXED_GRADED
!======================================================================!
!======================================================================!
!                               PROGRAM                                !
!======================================================================!
!======================================================================!
      BED_MIXED_GRADED=.FALSE.
!
!  ------ BED COMPOSITION
!
      CALL OS('X=Y-Z   ',X=ELAY,Y=ZF,Z=ZR)
!
!     ONLY ONE CLASS
!
      IF(NSICLA.EQ.1) THEN
        DO I=1,NPOIN
          AVAIL(I,1,1) = 1.D0
          ACLADM%R(I) = FDM(1)
        ENDDO
!       PURE MUD ONLY
        IF(SEDCO(1)) CALL INIT_MIXTE(XMVS,NPOIN,AVAIL,NSICLA,ES,
     &                               ES_SABLE, ES_VASE,
     &                               ELAY%R,NOMBLAY,CONC_VASE,
     &                                MS_SABLE%R,MS_VASE%R,ZF%R,
     &                               ZR%R,AVA0,CONC,DEBU,.FALSE.)
!
      ELSE
        IF(BED_MIXED_GRADED) THEN 
! THIS CONDITION IMPLIES THAT NSAND>1 AND/OR NMUD>1
!
! INITIALISATION OF RATIO_SAND AND RATIO_MUD 
!
          IF(NSAND.EQ.1) THEN
            DO ILAYER=1,NOMBLAYER
              DO I=1,NPOIN
                RATIO_SAND(NSAND,ILAYER,I)=1.D0
              ENDDO
            ENDDO
          ELSEIF(NMUD.EQ.1) THEN
            DO ILAYER=1,NOMBLAYER
              DO I=1,NPOIN
                RATIO_MUD(NMUD,ILAYER,I)=1.D0
              ENDDO
            ENDDO
          ENDIF
! 1) USERS DEFINES WHICH RATIO IN WHICH LAYER
!
!    FOR SAND
!    FIRST LAYER
!          DO I=1,NPOIN
!            RATIO_SAND(1,1,I)=...
!            RATIO_SAND(2,1,I)=...
!            ...
!            RATIO_SAND(NSAND,1,I)=...
!
!    SECOND LAYER
!            RATIO_SAND(1,2,I)=...
!            RATIO_SAND(2,2,I)=...
!            ...
!            RATIO_SAND(NSAND,2,I)=...
!            
!    NOMBLAYER 
!            RATIO_SAND(1,NOMBLAY,I)=...
!            RATIO_SAND(2,NOMBLAY,I)=...
!            ...
!            RATIO_SAND(NSAND,NOMBLAY,I)=...
!    FOR MUD 
!    FIRST LAYER
!            RATIO_MUD(1,1,I)=...
!            RATIO_MUD(2,1,I)=...
!            ...
!            RATIO_MUD(NSAND,1,I)=...
!
!    SECOND LAYER
!            RATIO_MUD(1,2,I)=...
!            RATIO_MUD(2,2,I)=...
!            ...
!            RATIO_MUD(NSAND,2,I)=...
!            
!    NOMBLAYER 
!            RATIO_MUD(1,NOMBLAY,I)=...
!            RATIO_MUD(2,NOMBLAY,I)=...
!            ...
!            RATIO_MUD(NMUD,NOMBLAY,I)=...
!
!          ENDDO  
!         
!  CHECK THE RATIOS
          DO I=1,NPOIN
            DO ILAYER=1,NLAYER
              CHECK_RS=0.D0
              CHECK_RM=0.D0
              DO ISAND=1,NSAND
                CHECK_RS=CHECK_RS+RATIO_SAND(ISAND,ILAYER,I)
              ENDDO
              IF(ABS(CHECK_RS-1.D0).GE.1.D-8) THEN 
                 WRITE(LU,*)'SUM OF SAND RATE COEFF MUST BE EQUAL TO 1!'
                 CALL PLANTE(1)
                 STOP
              ENDIF
              DO IMUD=1,NMUD
                CHECK_RM=CHECK_RM+RATIO_MUD(IMUD,ILAYER,I)
              ENDDO
              IF(ABS(CHECK_RM-1.D0).GE.1.D-8) THEN 
                 WRITE(LU,*)'SUM OF MUD RATE COEFF MUST BE EQUAL TO 1!'
                 CALL PLANTE(1)
                 STOP
              ENDIF
            ENDDO
          ENDDO
!
! 2) USERS DEFINE LAYER THICKNESS
!   a. every layer has the same thickness
!          DO ILAYER=1,NOMBLAY
!            DO IPOIN=1,NPOIN
!               ES(IPOIN,ILAYER)=...
!            ENDDO
!          ENDDO
!   b. different thickness for layers
!          DO IPOIN=1,NPOIN
!            ES(IPOIN,1)=...
!            ES(IPOIN,2)=...
!            ...
!            ES(IPOIN,NOMBLAY)=...
!          ENDDO 
!
! 3) USERS DEFINE RATIO_MUD_SAND
!
!   a. every layer has the same ratio_mud_sand
!          DO ILAYER=1,NOMBLAY
!            DO IPOIN=1,NPOIN
!               RATIO_MUD_SAND(ILAYER,IPOIN)=...
!            ENDDO
!          ENDDO
!   b. different ratio for layers
!          DO IPOIN=1,NPOIN
!            RATIO_MUD_SAND(1,IPOIN)=...
!            RATIO_MUD_SAND(2,IPOIN)=... 
!            ...
!            RATIO_MUD_SAND(NOMBLAY,IPOIN)=...
!          ENDDO 
!
! 4) MASS COMPUTATION
!
        DO I=1,NPOIN
          DO ILAYER = 1,NOMBLAY
            RATIO_MUD_VOL=RATIO_MUD_SAND(ILAYER,I)*XMVS
     &      /(CONC_VASE(ILAYER)*XKV*(1.D0-RATIO_MUD_SAND(ILAYER,I)))
            IF (RATIO_MUD_VOL.GT.1.D0) THEN
              MASS_TOT = ES(I,ILAYER)/
     &        ((1.D0-RATIO_MUD_SAND(ILAYER,I))/XMVS+
     &        RATIO_MUD_SAND(ILAYER,I)/CONC_VASE(ILAYER))
              DO IMUD = 1,NMUD
                MASS_MUD(I,ILAYER,IMUD) = 
     &          RATIO_MUD_SAND(ILAYER,I)*RATIO_MUD(IMUD,ILAYER,I)
     &          *MASS_TOT
              ENDDO
              DO ISAND = 1,NSAND
                MASS_SAND(ISAND,ILAYER,I) =
     &          (1.D0-RATIO_MUD_SAND(ILAYER,I))*MASS_TOT*
     &          RATIO_SAND(ISAND,ILAYER,I)
              ENDDO
            ELSE ! all mud fits inside porosity 
              MASS_SAND_TOT = XMVS*(1.D0-XKV)*ES(I,ILAYER)
              MASS_MUD_TOT = MASS_SAND_TOT/
     &                     (1.D0-RATIO_MUD_SAND(ILAYER,I))
            ENDIF
            DO IMUD = 1,NMUD
              MASS_MUD(I,ILAYER,IMUD) = 
     &        RATIO_MUD(IMUD,ILAYER,I)*MASS_MUD_TOT
            ENDDO
            DO ISAND = 1,NSAND
              MASS_SAND(ISAND,ILAYER,I) = 
     &        RATIO_SAND(ISAND,ILAYER,I)*MASS_SAND_TOT
            ENDDO
          ENDDO
        ENDDO
!
        ENDIF !(BED_MIX_GRADED)
!
!     NON-COHESIVE, MULTI-CLASSES
!
        IF(.NOT.MIXTE) THEN
! 
!
          CALL INIT_AVAI
!         CALL MEAN_GRAIN_SIZE
!         THIS PART CAN BE INTEGRATED INTO INIT_AVAI
          DO J=1,NPOIN
            ACLADM%R(J) = 0.D0
            UNLADM%R(J) = 0.D0
            DO I=1,NSICLA
              IF(AVAIL(J,1,I).GT.0.D0) THEN
                ACLADM%R(J) = ACLADM%R(J) + FDM(I)*AVAIL(J,1,I)
                UNLADM%R(J) = UNLADM%R(J) + FDM(I)*AVAIL(J,2,I)
              ENDIF
            ENDDO
            ACLADM%R(J)=MAX(ACLADM%R(J),0.D0)
            UNLADM%R(J)=MAX(UNLADM%R(J),0.D0)
          ENDDO
        ELSE
!
          CALL INIT_MIXTE(XMVS,NPOIN,AVAIL,NSICLA,ES,
     &               ES_SABLE, ES_VASE, ELAY%R,
     &               NOMBLAY,CONC_VASE,MS_SABLE%R,
     &               MS_VASE%R,ZF%R,ZR%R,AVA0,CONC,DEBU,MIXTE)
          DO I=1,NPOIN
            ACLADM%R(I) = FDM(1)
          ENDDO
        ENDIF
!
      ENDIF !(NSICLA)
!
      IF(LGRAFED) THEN
        DO I=1,NSICLA
          FRACSED_GF(I)=AVA0(I)
        ENDDO
      ENDIF
!
!     SETTLING VELOCITY
!
      IF(.NOT.CALWC) THEN
        DENS = (XMVS - XMVE) / XMVE
        DO I = 1, NSICLA
          CALL VITCHU_SISYPHE(XWC(I),DENS,FDM(I),GRAV,VCE)
          IF(BED_MIXED_GRADED) THEN
!            calcul de la vitesse de chute en 2d ou au fond
!            si T3D: la vitesse de chute est calculée en 3D puis repassée à sisyphe pour le fond		
            IF(TYPE_OF_SEDIMENT.EQ.SED_NCO) THEN
              CALL VITCHU_SISYPHE(XWC(I),DENS,FDM(I),GRAV,VCE)
            ELSE
!           par defaut 1mm/s pour toutes les classes de vase
              XWC(I)= 0.001D0
            ENDIF
          ENDIF
        ENDDO
      ENDIF
!
!     SHIELDS PARAMETER
!
      IF(.NOT.CALAC) THEN
        DENS  = (XMVS - XMVE )/ XMVE
        DO I = 1, NSICLA
          DSTAR = FDM(I)*(GRAV*DENS/VCE**2)**(1.D0/3.D0)
          IF (DSTAR <= 4.D0) THEN
            AC(I) = 0.24D0/DSTAR
          ELSEIF (DSTAR <= 10.D0) THEN
            AC(I) = 0.14D0*DSTAR**(-0.64D0)
          ELSEIF (DSTAR <= 20.D0) THEN
            AC(I) = 0.04D0*DSTAR**(-0.1D0)
!           CORRECTION 27/06/2016
!          ELSEIF (DSTAR <= 150.D0) THEN
          ELSEIF (DSTAR <= 72.D0) THEN 
            AC(I) = 0.013D0*DSTAR**0.29D0
          ELSE
!           CORRECTION 30/05/2012
!           AC(I) = 0.055D0
            AC(I) = 0.045D0
          ENDIF
        ENDDO
      ENDIF
!
!     FOR MIXED SEDIMENTS
!
      IF(MIXTE) TOCE_SABLE=AC(1)*FDM(1)*GRAV*(XMVS - XMVE)
     
      IF(BED_MIXED_GRADED) THEN
        DO ISAND = 1,NSAND
          TOCE_SAND(ISAND)= AC(ISAND)*FDM(ISAND)*GRAV*(XMVS - XMVE)
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
