!                    ************************
                     SUBROUTINE INIT_SEDIMENT
!                    ************************
!
     &(NSICLA,ELAY,ZF,ZR,NPOIN,AVAIL,AVA0,
     & CALWC,XMVS,XMVE,GRAV,VCE,XWC,FDM,
     & CALAC,AC,SEDCO,ES,ES_SABLE, ES_VASE ,NOMBLAY,CONC_MUD,
     & MS_SABLE,MS_VASE,ACLADM,UNLADM,TOCE_SABLE,
     & CONC,NLAYER,DEBU,MIXTE,NEW_BED_MODEL,
     & TOC_MUD,TOCE_MUD,VOLU2D)
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
!| CONC_MUD       |<->| MUD CONCENTRATION FOR EACH LAYER
!| ELAY           |<->| THICKNESS OF SURFACE LAYER
!| ES             |<->| LAYER THICKNESSES AS DOUBLE PRECISION
!| ES_SABLE       |<->| LAYER THICKNESSES OF SAND AS DOUBLE PRECISION
!| ES_VASE        |<->| LAYER THICKNESSES OF MUD AS DOUBLE PRECISION
!| FDM            |-->| DIAMETER DM FOR EACH CLASS
!| GRAV           |-->| ACCELERATION OF GRAVITY
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
!| VOLU2D         |-->| INTEGRAL OF TEST FUNCTIONS (NOT ASSEMBLED IN //)
!| XMVE           |-->| FLUID DENSITY
!| XMVS           |-->| SEDIMENT DENSITY
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
      USE DECLARATIONS_SISYPHE, ONLY: MASS_MUD,MASS_SAND,
     & MASS_MUD_TOT,MASS_SAND_TOT,MASS_MIX_TOT,
     & RATIO_SAND,RATIO_MUD,RATIO_MUD_SAND,NSAND,NMUD,
     & TYPE_SED,XKV,TOCE_SAND,MASSTOT,NUM_IMUD_ICLA,NUM_ISAND_ICLA
      USE INTERFACE_PARALLEL, ONLY: P_DSUM
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER,           INTENT(IN)     :: NSICLA,NPOIN,NOMBLAY
      TYPE(BIEF_OBJ),    INTENT(INOUT)  :: ELAY,ZF,ZR
      TYPE(BIEF_OBJ), INTENT(INOUT)     :: MS_SABLE, MS_VASE
      TYPE(BIEF_OBJ),    INTENT(INOUT)  :: ACLADM, UNLADM
      TYPE(BIEF_OBJ),    INTENT(INOUT)  :: NLAYER,VOLU2D
      LOGICAL,           INTENT(IN)     :: CALWC
      LOGICAL,           INTENT(IN)     :: CALAC
      DOUBLE PRECISION,  INTENT(IN)     :: XMVS,XMVE,GRAV,VCE
      DOUBLE PRECISION,  INTENT(INOUT)  :: AVA0(NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: AVAIL(NPOIN,NOMBLAY,NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: FDM(NSICLA),XWC(NSICLA)
      DOUBLE PRECISION,  INTENT(INOUT)  :: AC(NSICLA),TOCE_SABLE
      LOGICAL,           INTENT(IN)     :: SEDCO(NSICLA), DEBU
      LOGICAL,           INTENT(IN)     :: MIXTE
      DOUBLE PRECISION, INTENT(IN)      :: CONC_MUD(NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT)   :: ES(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT)   :: ES_SABLE(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT)   :: ES_VASE(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT)   :: CONC(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(IN)      :: TOC_MUD(NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT)   :: TOCE_MUD(NOMBLAY,NPOIN)
      LOGICAL,           INTENT(IN)     :: NEW_BED_MODEL
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER            :: I,J,K
      INTEGER            :: ILAYER,IMUD,ISAND
      DOUBLE PRECISION   :: DENS,DSTAR
      DOUBLE PRECISION   :: CHECK_RS,CHECK_RM
      DOUBLE PRECISION   :: TERM,DISCR
      DOUBLE PRECISION   :: MASS_TOT
!======================================================================!
!======================================================================!
!                               PROGRAM                                !
!======================================================================!
!======================================================================!
!
!  ------ BED COMPOSITION
!
      CALL OS('X=Y-Z   ',X=ELAY,Y=ZF,Z=ZR)
!
!     ONLY ONE CLASS : OLD BED MODEL AND NEW BED MODEL
!
      IF(.NOT.NEW_BED_MODEL) THEN
        IF(NSICLA.EQ.1) THEN
          DO I=1,NPOIN
            AVAIL(I,1,1) = 1.D0
            ACLADM%R(I) = FDM(1)
          ENDDO
!         PURE MUD ONLY
          IF(SEDCO(1)) CALL INIT_MIXTE(XMVS,NPOIN,AVAIL,NSICLA,ES,
     &                               ES_SABLE, ES_VASE,
     &                               ELAY%R,NOMBLAY,CONC_MUD,
     &                                MS_SABLE%R,MS_VASE%R,ZF%R,
     &                               ZR%R,AVA0,CONC,DEBU,.FALSE.)
        ENDIF
!
      ELSE
!     NEW BED MODEL!
!
! INITIALISATION OF RATIO_SAND AND RATIO_MUD
!
        IF(NSAND.EQ.1) THEN
          DO I = 1,NPOIN
            DO ILAYER = 1,NOMBLAY
              RATIO_SAND(NSAND,ILAYER,I)=1.D0
            ENDDO
          ENDDO
        ENDIF
        IF(NMUD.EQ.1) THEN
          DO I = 1,NPOIN
            DO ILAYER = 1,NOMBLAY
              RATIO_MUD(NMUD,ILAYER,I)=1.D0
            ENDDO
          ENDDO
        ENDIF
!
! 1) USERS DEFINES WHICH RATIO IN WHICH LAYER (OR WHICH MASS IN WHICH LAYER, to do)
!
!    FOR SAND
!    FIRST LAYER
!        DO I=1,NPOIN
!          RATIO_SAND(1,1,I)=...
!          RATIO_SAND(2,1,I)=...
!          ...
!          RATIO_SAND(NSAND,1,I)=...
!
!    SECOND LAYER
!          RATIO_SAND(1,2,I)=...
!          RATIO_SAND(2,2,I)=...
!          ...
!          RATIO_SAND(NSAND,2,I)=...
!
!    NOMBLAY
!          RATIO_SAND(1,NOMBLAY,I)=...
!          RATIO_SAND(2,NOMBLAY,I)=...
!          ...
!          RATIO_SAND(NSAND,NOMBLAY,I)=...
!
!    FOR MUD
!
!    FIRST LAYER
!          RATIO_MUD(1,1,I)=...
!          RATIO_MUD(2,1,I)=...
!          ...
!          RATIO_MUD(NSAND,1,I)=...
!
!    SECOND LAYER
!          RATIO_MUD(1,2,I)=...
!          RATIO_MUD(2,2,I)=...
!          ...
!          RATIO_MUD(NSAND,2,I)=...
!
!    NOMBLAY
!          RATIO_MUD(1,NOMBLAY,I)=...
!          RATIO_MUD(2,NOMBLAY,I)=...
!          ...
!          RATIO_MUD(NMUD,NOMBLAY,I)=...
!
!        ENDDO
!
!  CHECK THE RATIOS
!  UNCOMMENT THIS PART IF YOU HAVE DEFINED THE RATIO_SAND AND RATIO_MUD
!        DO I=1,NPOIN
!          DO ILAYER=1,NOMBLAY
!            CHECK_RS=0.D0
!            CHECK_RM=0.D0
!            IF(NSAND.GT.0) THEN
!              DO ISAND=1,NSAND
!                CHECK_RS=CHECK_RS+RATIO_SAND(ISAND,ILAYER,I)
!              ENDDO
!              IF(ABS(CHECK_RS-1.D0).GE.1.D-8) THEN
!                 WRITE(LU,*)'SUM OF SAND RATE COEFF MUST BE EQUAL TO 1!'
!                 WRITE(LU,*)'VERIFY YOUR MASS RATE OF SAND'
!                 CALL PLANTE(1)
!                 STOP
!              ENDIF
!            ENDIF
!            IF(NMUD.GT.0) THEN
!              DO IMUD=1,NMUD
!                CHECK_RM=CHECK_RM+RATIO_MUD(IMUD,ILAYER,I)
!              ENDDO
!              IF(ABS(CHECK_RM-1.D0).GE.1.D-8) THEN
!                 WRITE(LU,*)'SUM OF MUD RATE COEFF MUST BE EQUAL TO 1!'
!                 WRITE(LU,*)'VERIFY YOUR MASS RATE OF MUD'
!                 CALL PLANTE(1)
!                 STOP
!              ENDIF
!            ENDIF
!          ENDDO
!        ENDDO
!
! 2) INITIALISATION OF MEAN DIMAETER FOR ACTIVE LAYER AND SUBSTRATUM
!
        IF(NSAND.EQ.1.AND.NMUD.EQ.0) THEN
          DO I = 1,NPOIN
! IF ONLY 1 SAND: MEAN DIAMETER OF THE ACTIVE LAYER IS THE ONE GIVEN IN THE STEERING FILE
! is it the case even if nmud.ne.0? if yes, this case can be included in the
! loop above
            ACLADM%R(I) = FDM(1)
          ENDDO
        ENDIF
!
        IF(NSAND.GT.1.AND.NMUD.EQ.0) THEN
! NOTE THAT IN THIS CASE NSAND=NSICLA
! what happen if nmud.gt.0 for mean diameter?
          DO J=1,NPOIN
            ACLADM%R(J) = 0.D0
            UNLADM%R(J) = 0.D0
            DO I=1,NSAND
              IF(RATIO_SAND(I,1,J).GT.0.D0) THEN
                ACLADM%R(J) = ACLADM%R(J) + FDM(I)*RATIO_SAND(I,1,J)
              ENDIF
            ACLADM%R(J)=MAX(ACLADM%R(J),0.D0)
            ENDDO
          ENDDO
          IF(NOMBLAY.GT.1) THEN
            DO J=1,NPOIN
              DO I=1,NSAND
                IF(RATIO_SAND(I,2,J).GT.0.D0) THEN
                  UNLADM%R(J) = UNLADM%R(J) + FDM(I)*RATIO_SAND(I,2,J)
                ENDIF
                UNLADM%R(J)=MAX(UNLADM%R(J),0.D0)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
!
! 3) USERS DEFINE LAYER THICKNESS
!
!    INITIALISATION:
!    IF USER DOES NOT DEFINE ES, THEN ES WILL BE (ZF-ZR)/NOMBLAY
        DO ILAYER=1,NOMBLAY
          DO I=1,NPOIN
            ES(I,ILAYER)=ELAY%R(I)/NOMBLAY
          ENDDO
        ENDDO
!   a. every layer has the same thickness
!        DO ILAYER=1,NOMBLAY
!          DO I=1,NPOIN
!             ES(I,ILAYER)=...
!          ENDDO
!        ENDDO
!   b. different thickness for layers
!        DO I=1,NPOIN
!          ES(I,1)=...
!          ES(I,2)=...
!          ...
!          ES(I,NOMBLAY)=...
!        ENDDO
!
! to add: if users define es(..,..) then check that sum(es)=zf-zr!!
!
!    FOR SEVERAL SANDS AND LAYERS.GT.2 : ELAY IS THE FIRST UPPER LAYER
        IF(NSAND.GT.1.AND.NOMBLAY.GT.1) THEN
          DO I=1,NPOIN
            ELAY%R(I)=ES(I,1)
! OPTION ACTIVE LAYER NOT CONSTANT: TO DO!
!            ELAY%R(I)=3.D0 * ACLADM%R(J)
          ENDDO
        ENDIF
! 4) USERS DEFINE RATIO_MUD_SAND
!
!     INITIALISATION
        IF(NMUD.EQ.0.AND.NSAND.GE.1) THEN
          DO I=1,NPOIN
            DO ILAYER=1,NOMBLAY
              RATIO_MUD_SAND(ILAYER,I) = 0.D0
            ENDDO
          ENDDO
        ELSEIF(NSAND.EQ.0.AND.NMUD.GE.1) THEN
          DO I=1,NPOIN
            DO ILAYER=1,NOMBLAY
              RATIO_MUD_SAND(ILAYER,I) = 1.D0
            ENDDO
          ENDDO
        ENDIF
!   a. every layer has the same ratio_mud_sand
!        DO ILAYER=1,NOMBLAY
!          DO I=1,NPOIN
!             RATIO_MUD_SAND(ILAYER,I)=...
!          ENDDO
!        ENDDO
!   b. different ratio for layers
!        DO I=1,NPOIN
!          RATIO_MUD_SAND(1,I)=...
!          RATIO_MUD_SAND(2,I)=...
!          ...
!          RATIO_MUD_SAND(NOMBLAY,I)=...
!        ENDDO
!
! 5) MASS COMPUTATION : ALL MASSES IN [kg/m^2]
!     -> MASS_MUD;MASS_SAND;MASS_MUD_TOT;MASS_SAND_TOT
!
!    WARNING : FOR MASS COMPUTATION IS MANDATORY RATIO_MUD_SAND,
!              XMVS,CONC_MUD,XKV,ES(ILAYER)
!
        DO I = 1,NPOIN
          DO ILAYER = 1,NOMBLAY
            TERM=RATIO_MUD_SAND(ILAYER,I)/CONC_MUD(ILAYER)-
     &      (XKV(ILAYER)*(1.D0-RATIO_MUD_SAND(ILAYER,I)))/
     &      (XMVS*(1.D0-XKV(ILAYER)))
!    TERM REPRESENTS THE DIFFERENCE BETWEEN MUD VOLUME AND VOID VOLUME
!    DISCR IS POSITIVE WHEN MUD FILLS ALL THE SAND POROSITY OTHERWISE IS ZERO
            DISCR=MAX(0.D0,TERM)
            MASS_TOT = ES(I,ILAYER)/
     &      ((1.D0-RATIO_MUD_SAND(ILAYER,I))/(XMVS*(1.D0-XKV(ILAYER)))
     &      +DISCR)
!
            MASS_MUD_TOT(ILAYER,I) = RATIO_MUD_SAND(ILAYER,I)
     &      *MASS_TOT
            MASS_SAND_TOT(ILAYER,I) = (1.D0-RATIO_MUD_SAND(ILAYER,I))
     &      *MASS_TOT
!   COMPUTES MASS FOR EVERY MUD
            DO IMUD = 1,NMUD
              MASS_MUD(IMUD,ILAYER,I) = MASS_MUD_TOT(ILAYER,I)
     &        *RATIO_MUD(IMUD,ILAYER,I)
            ENDDO
!   COMPUTES MASS FOR EVERY SAND
            DO ISAND = 1,NSAND
              MASS_SAND(ISAND,ILAYER,I) = MASS_SAND_TOT(ILAYER,I)
     &        *RATIO_SAND(ISAND,ILAYER,I)
            ENDDO
          ENDDO
        ENDDO
!
! 6) REAL MASS COMPUTATION: FROM [kg/m2] to [kg]
! USEFUL FOR FIRST LISTING AND MASS BALANCE
!
        DO I=1,NSICLA
          MASSTOT(I)=0.D0
        ENDDO
!
        IF(NSAND.GT.0) THEN
          DO I=1,NPOIN
            DO ILAYER=1,NOMBLAY
              DO ISAND=1,NSAND
                K=NUM_ISAND_ICLA(ISAND)
                MASSTOT(K)=MASSTOT(K)+MASS_SAND(ISAND,ILAYER,I)
     &                     *VOLU2D%R(I)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NMUD.GT.0) THEN
          DO I=1,NPOIN
            DO ILAYER=1,NOMBLAY
              DO IMUD=1,NMUD
                K=NUM_IMUD_ICLA(IMUD)
                MASSTOT(K)=MASSTOT(K)+MASS_MUD(IMUD,ILAYER,I)
     &                     *VOLU2D%R(I)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
      IF(NCSIZE.GT.1) THEN
        DO I=1,NSICLA
          MASSTOT(I) = P_DSUM(MASSTOT(I))
        ENDDO
      ENDIF
!
      DO I=1,NSICLA
        IF(LNG.EQ.1) THEN
          WRITE(LU,*)'MASSE TOTALE DE LA CLASSE ',I ,' :',MASSTOT(I)
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*)'TOTAL MASS OF CLASS ',I ,' :',MASSTOT(I)
        ENDIF
      ENDDO
!
      ENDIF !(NEW BED MODEL)
!
      IF(.NOT.NEW_BED_MODEL) THEN
!     NON-COHESIVE, MULTI-CLASSES
!
        IF(.NOT.MIXTE) THEN
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
     &               NOMBLAY,CONC_MUD,MS_SABLE%R,
     &               MS_VASE%R,ZF%R,ZR%R,AVA0,CONC,DEBU,MIXTE)
          DO I=1,NPOIN
            ACLADM%R(I) = FDM(1)
          ENDDO
        ENDIF
!
      ELSE
!     NEW_BED_MODEL - Nothing, isn't it?
!
      ENDIF !(NEW_BED_MODEL)
!
!     SETTLING VELOCITY
!
      DENS = (XMVS - XMVE) / XMVE
      DO I = 1, NSICLA
        IF(XWC(I).LT.-1) THEN
! SETTLING VELOCITY IS NOT GIVEN IN THE PARAMETER FILE
          CALL VITCHU_SISYPHE(XWC(I),DENS,FDM(I),GRAV,VCE)
          IF(NEW_BED_MODEL) THEN
!            calcul de la vitesse de chute en 2d ou au fond
!            si T3D: la vitesse de chute est calculée en 3D puis repassée à sisyphe pour le fond
            IF(TYPE_SED(I).EQ.'NCO') THEN
              CALL VITCHU_SISYPHE(XWC(I),DENS,FDM(I),GRAV,VCE)
            ELSE
!           par defaut 1mm/s pour toutes les classes de vase
              XWC(I)= 0.001D0
            ENDIF
          ENDIF
        ENDIF
!
!     SHIELDS PARAMETER
!
        IF(AC(I).LT.-1) THEN
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
        ENDIF
      ENDDO
!
!     FOR MIXED SEDIMENTS
!
      IF(MIXTE) TOCE_SABLE = AC(1)*FDM(1)*GRAV*(XMVS - XMVE)

      IF(NEW_BED_MODEL) THEN
        DO ISAND = 1,NSAND
          TOCE_SAND(ISAND) = AC(ISAND)*FDM(ISAND)*GRAV*(XMVS - XMVE)
        ENDDO
        DO I = 1,NPOIN
          DO ILAYER=1,NOMBLAY
            TOCE_MUD(ILAYER,I) = TOC_MUD(ILAYER)
          ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
