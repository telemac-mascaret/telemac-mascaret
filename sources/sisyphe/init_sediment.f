!                    ************************
                     SUBROUTINE INIT_SEDIMENT
!                    ************************
!
     &(NSICLA,ELAY,ZF,ZR,NPOIN,AVAIL,AVA0,
     & CALWC,XMVS,XMVE,GRAV,VCE,XWC,FDM,
     & CALAC,AC,SEDCO,ES,ES_SABLE, ES_VASE ,NOMBLAY,CONC_MUD,
     & MS_SABLE,MS_VASE,ACLADM,UNLADM,TOCE_SABLE,
     & CONC,NLAYER,DEBU,MIXTE,NEW_BED_MODEL,
     & TOC_MUD,TOCE_MUD,VOLU2D,X,Y)
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
     & TYPE_SED,XKV,TOCE_SAND,MASSTOT,NUM_IMUD_ICLA,NUM_ISAND_ICLA,
     & ELAY0,NOMBSTRAT	 
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
      DOUBLE PRECISION,  INTENT(IN)     :: X(NPOIN),Y(NPOIN)
      !
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER            :: I,J,K,IPOIN
      INTEGER            :: ILAYER,IMUD,ISAND,ICLA,ISTRAT
      DOUBLE PRECISION   :: DENS,DSTAR
      DOUBLE PRECISION   :: CHECK_RS,CHECK_RM
      DOUBLE PRECISION   :: TERM,DISCR
      DOUBLE PRECISION   :: MASS_TOT
      DOUBLE PRECISION   :: tot,tot2,totsand,totmud
      DOUBLE PRECISION   :: RATIO_INIT(NSICLA,NOMBSTRAT,NPOIN)	  
      DOUBLE PRECISION   :: ESTRATUM(NPOIN,NOMBSTRAT)	  	  
!======================================================================!
!======================================================================!
!                               PROGRAM                                !
!======================================================================!
!======================================================================!
!
!  ------ BED COMPOSITION
!
!
!     ONLY ONE CLASS : OLD BED MODEL AND NEW BED MODEL
!
      IF(.NOT.NEW_BED_MODEL) THEN
        CALL OS('X=Y-Z   ',X=ELAY,Y=ZF,Z=ZR)
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
! 1) USERS DEFINE INITIAL LAYER THICKNESS AND INITIAL RATIOS FOR EACH CLASS
      CALL BED_INIT_USER(ESTRATUM,RATIO_INIT)
!
!    INITIALISATION ES:
        DO IPOIN=1,NPOIN
         DO ILAYER=1,NOMBLAY
            ES(IPOIN,ILAYER)=0.D0
         ENDDO
        ENDDO
!
      IF (NSICLA.EQ.1)THEN
        IF (NOMBLAY.EQ.1)THEN
            DO IPOIN=1,NPOIN
                ES(IPOIN,NOMBLAY)= ESTRATUM(IPOIN,NOMBLAY) 
            ENDDO
            ELSE
            WRITE(LU,*) 'NOMBLAY > 1 and only one sediment class'
            CALL PLANTE(1)
            STOP
        ENDIF
      ELSE
!   first layer is active layer
            DO IPOIN=1,NPOIN
                ELAY%R(IPOIN)= ELAY0 ! on pourra mettre une autre option ( 3 * Dm) a gerer dans bed_udate_hirano
                ES(IPOIN,1) = 0.D0 ! la couche active est remplie dans le premeier appel de bed_update_hirano
                DO ILAYER = 2,NOMBLAY
                   ISTRAT=ILAYER-1
                   ES(IPOIN,ILAYER) = ESTRATUM(IPOIN,ISTRAT)
                ENDDO
            ENDDO            
      ENDIF
       DO IPOIN=1,NPOIN 
         ZR%R(IPOIN)=ZF%R(IPOIN)
         DO ILAYER=1,NOMBLAY
           ZR%R(IPOIN)=ZR%R(IPOIN)-ES(IPOIN,ILAYER)
         ENDDO
       ENDDO
! 2) USERS DEFINES: INITIALISATION OF RATIO_SAND 
!************FAIRE LE CALCUL DES RATIO_SAND, RATIO_MUD et RATIO_MUD_SAND à partir des RATIO_INIT(NCLA,NOMBSTRAT,NPOIN)

        DO IPOIN=1,NPOIN
            DO ISTRAT=1,NOMBSTRAT
                tot=0.D0
                tot2=0.D0
                DO ISAND=1,NSAND
                tot=tot+RATIO_INIT(NUM_ISAND_ICLA(ISAND),ISTRAT,IPOIN)
                ENDDO
                DO ISAND=1,NSAND
                    IF(NSICLA.GT.1)THEN
                        ILAYER = ISTRAT+1
!la couche active est remplie dans le premier appel de bed_update_hirano, mais il faut quand même initialiser ses ratios ...
                        IF(ISAND.LT.NSAND)THEN
						  RATIO_SAND(ISAND,1,IPOIN)=0.D0
                        ELSE
						  RATIO_SAND(NSAND,1,IPOIN)=1.D0
                        ENDIF						
                    ELSE
                        ILAYER=ISTRAT
                    ENDIF
                    IF(ISAND.LT.NSAND)THEN
                        IF(tot.GT.0.D0)THEN
                            RATIO_SAND(ISAND,ILAYER,IPOIN)=
     &                RATIO_INIT(NUM_ISAND_ICLA(ISAND),ISTRAT,IPOIN)/tot
                        ELSE
                            RATIO_SAND(NSAND,ILAYER,IPOIN)=0.D0
                        ENDIF
                        tot2=tot2+RATIO_SAND(ISAND,ILAYER,IPOIN)
                    ELSE
                        RATIO_SAND(NSAND,ILAYER,IPOIN)=1.D0-tot2
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
! 3) INITIALISATION OF RATIO_MUD 
        DO IPOIN=1,NPOIN
            DO ISTRAT=1,NOMBSTRAT
                tot=0.D0
                tot2=0.D0
                DO IMUD=1,NMUD
                tot= tot + RATIO_INIT(NUM_IMUD_ICLA(IMUD),ISTRAT,IPOIN)
                ENDDO
                DO IMUD=1,NMUD
                    IF(NSICLA.GT.1)THEN
                        ILAYER = ISTRAT+1
!la couche active est remplie dans le premier appel de bed_update_hirano, mais il faut quand même initialiser ses ratios ...
                        IF(IMUD.LT.NMUD)THEN
						  RATIO_MUD(IMUD,1,IPOIN)=0.D0
                        ELSE
						  RATIO_MUD(NMUD,1,IPOIN)=1.D0
                        ENDIF	
                    ELSE
                        ILAYER=ISTRAT
                    ENDIF
                    IF(IMUD.LT.NMUD)THEN
                        IF(tot.GT.0.D0)THEN
                            RATIO_MUD(IMUD,ILAYER,IPOIN)=
     &                  RATIO_INIT(NUM_IMUD_ICLA(IMUD),ISTRAT,IPOIN)/tot
                        ELSE
                            RATIO_MUD(NMUD,ILAYER,IPOIN)=0.D0
                        ENDIF
                        tot2=tot2+RATIO_MUD(IMUD,ILAYER,IPOIN)
                    ELSE
                        RATIO_MUD(NMUD,ILAYER,IPOIN)=1.D0-tot2
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
! 4) INITIALISATION OF RATIO_MUD_SAND
        DO IPOIN=1,NPOIN
            DO ISTRAT=1,NOMBSTRAT
                totmud=0.D0
                totsand=0.D0
                DO IMUD=1,NMUD
                totmud= totmud + 
     &          RATIO_INIT(NUM_IMUD_ICLA(IMUD),ISTRAT,IPOIN)
                ENDDO
                DO ISAND=1,NSAND
                totsand= totsand + 
     &          RATIO_INIT(NUM_ISAND_ICLA(ISAND),ISTRAT,IPOIN)
                ENDDO               
                IF(NSICLA.GT.1)THEN
                ILAYER = ISTRAT+1
!la couche active est remplie dans le premeier appel de bed_update_hirano
                ELSE
                    ILAYER=ISTRAT
                ENDIF
                IF(totmud+totsand.GT.0.D0)THEN
                   RATIO_MUD_SAND(ILAYER,IPOIN)=totmud/(totmud+totsand)
                ELSE
                   RATIO_MUD_SAND(ILAYER,IPOIN)=1.D0
                ENDIF
            ENDDO
        ENDDO       
!
!  CHECK THE RATIOS
        DO IPOIN=1,NPOIN
          DO ILAYER=1,NOMBLAY
            CHECK_RS=0.D0
            CHECK_RM=0.D0
            IF(NSAND.GT.0) THEN
              DO ISAND=1,NSAND
                CHECK_RS=CHECK_RS+RATIO_SAND(ISAND,ILAYER,IPOIN)
              ENDDO
              IF(ABS(CHECK_RS-1.D0).GE.1.D-8) THEN
                 WRITE(LU,*)'SUM OF SAND RATE COEFF MUST BE EQUAL TO 1!'
                 WRITE(LU,*)'VERIFY YOUR MASS RATE OF SAND'
                 CALL PLANTE(1)
                 STOP
              ENDIF
            ENDIF
            IF(NMUD.GT.0) THEN
              DO IMUD=1,NMUD
                CHECK_RM=CHECK_RM+RATIO_MUD(IMUD,ILAYER,IPOIN)
              ENDDO
              IF(ABS(CHECK_RM-1.D0).GE.1.D-8) THEN
                 WRITE(LU,*)'SUM OF MUD RATE COEFF MUST BE EQUAL TO 1!'
                 WRITE(LU,*)'VERIFY YOUR MASS RATE OF MUD'
                 CALL PLANTE(1)
                 STOP
              ENDIF
            ENDIF
            IF(NSAND.GT.0.AND.NMUD.GT.0) THEN
                 WRITE(LU,*)' you have to define 
     & RATIO_MUD_SAND(ILAY,IPOIN) in sediment_init.f'
                 CALL PLANTE(1)
                 STOP
            ENDIF
          ENDDO
       ENDDO
!
! 5) MASS COMPUTATION : ALL MASSES IN [kg/m^2]
!     -> MASS_MUD;MASS_SAND;MASS_MUD_TOT;MASS_SAND_TOT
!
!    WARNING : FOR MASS COMPUTATION IS MANDATORY RATIO_MUD_SAND,
!              XMVS,CONC_MUD,XKV,ES(ILAYER)
!
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NOMBLAY
            TERM=0.D0
            IF(NMUD.GT.0) THEN
              TERM=RATIO_MUD_SAND(ILAYER,IPOIN)/CONC_MUD(ILAYER)-
     &        (XKV(ILAYER,IPOIN)*(1.D0-RATIO_MUD_SAND(ILAYER,IPOIN)))/
     &        (XMVS*(1.D0-XKV(ILAYER,IPOIN)))
!    TERM REPRESENTS THE DIFFERENCE BETWEEN MUD VOLUME AND VOID VOLUME
            ENDIF
!    DISCR IS POSITIVE WHEN MUD FILLS ALL THE SAND POROSITY OTHERWISE IS ZERO
            DISCR=MAX(0.D0,TERM)
            MASS_TOT = ES(IPOIN,ILAYER)/
     &   ((1.D0-RATIO_MUD_SAND(ILAYER,IPOIN))/
     &   (XMVS*(1.D0-XKV(ILAYER,IPOIN)))
     &      +DISCR)
!
            MASS_MUD_TOT(ILAYER,IPOIN) =
     &      RATIO_MUD_SAND(ILAYER,IPOIN)*MASS_TOT
            MASS_SAND_TOT(ILAYER,IPOIN) =
     &      (1.D0-RATIO_MUD_SAND(ILAYER,IPOIN))*MASS_TOT
!   COMPUTES MASS FOR EVERY MUD
            DO IMUD = 1,NMUD
              MASS_MUD(IMUD,ILAYER,IPOIN) = MASS_MUD_TOT(ILAYER,IPOIN)
     &        *RATIO_MUD(IMUD,ILAYER,IPOIN)
            ENDDO
!   COMPUTES MASS FOR EVERY SAND
            DO ISAND = 1,NSAND
              MASS_SAND(ISAND,ILAYER,IPOIN) =
     &              MASS_SAND_TOT(ILAYER,IPOIN)
     &              *RATIO_SAND(ISAND,ILAYER,IPOIN)
            ENDDO
          ENDDO
        ENDDO
!
! 5) REAL MASS COMPUTATION: FROM [kg/m2] to [kg]
! USEFUL FOR FIRST LISTING AND MASS BALANCE
!
        DO ICLA=1,NSICLA
          MASSTOT(ICLA)=0.D0
        ENDDO
!
        IF(NSAND.GT.0) THEN
          DO IPOIN=1,NPOIN
            DO ILAYER=1,NOMBLAY
              DO ISAND=1,NSAND
                K=NUM_ISAND_ICLA(ISAND)
                MASSTOT(K)=MASSTOT(K)+MASS_SAND(ISAND,ILAYER,IPOIN)
     &                     *VOLU2D%R(IPOIN)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
        IF(NMUD.GT.0) THEN
          DO IPOIN=1,NPOIN
            DO ILAYER=1,NOMBLAY
              DO IMUD=1,NMUD
                K=NUM_IMUD_ICLA(IMUD)
                MASSTOT(K)=MASSTOT(K)+MASS_MUD(IMUD,ILAYER,IPOIN)
     &                     *VOLU2D%R(IPOIN)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
!
      IF(NCSIZE.GT.1) THEN
        DO ICLA=1,NSICLA
          MASSTOT(ICLA) = P_DSUM(MASSTOT(ICLA))
        ENDDO
      ENDIF
!
      DO ICLA=1,NSICLA
        IF(LNG.EQ.1) THEN
       WRITE(LU,*)'MASSE TOTALE DE LA CLASSE ',ICLA ,' :',MASSTOT(ICLA)
        ENDIF
        IF(LNG.EQ.2) THEN
          WRITE(LU,*)'TOTAL MASS OF CLASS ',ICLA ,' :',MASSTOT(ICLA)
        ENDIF
      ENDDO
!
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
      DO ICLA = 1, NSICLA
        IF(XWC(ICLA).LT.-1) THEN
! SETTLING VELOCITY IS NOT GIVEN IN THE PARAMETER FILE
          CALL VITCHU_SISYPHE(XWC(ICLA),DENS,FDM(ICLA),GRAV,VCE)
          IF(NEW_BED_MODEL) THEN
!            calcul de la vitesse de chute en 2d ou au fond
!            si T3D: la vitesse de chute est calculée en 3D puis repassée à sisyphe pour le fond
            IF(TYPE_SED(ICLA).EQ.'NCO') THEN
              CALL VITCHU_SISYPHE(XWC(ICLA),DENS,FDM(ICLA),GRAV,VCE)
            ELSE
!           par defaut 1mm/s pour toutes les classes de vase
              XWC(ICLA)= 0.001D0
            ENDIF
          ENDIF
        ENDIF
!
!     SHIELDS PARAMETER
!
        IF(AC(ICLA).LT.-1) THEN
          DSTAR = FDM(ICLA)*(GRAV*DENS/VCE**2)**(1.D0/3.D0)
          IF (DSTAR <= 4.D0) THEN
            AC(ICLA) = 0.24D0/DSTAR
          ELSEIF (DSTAR <= 10.D0) THEN
            AC(ICLA) = 0.14D0*DSTAR**(-0.64D0)
          ELSEIF (DSTAR <= 20.D0) THEN
            AC(ICLA) = 0.04D0*DSTAR**(-0.1D0)
!           CORRECTION 27/06/2016
!          ELSEIF (DSTAR <= 150.D0) THEN
          ELSEIF (DSTAR <= 72.D0) THEN
            AC(ICLA) = 0.013D0*DSTAR**0.29D0
          ELSE
!           CORRECTION 30/05/2012
!           AC(ICLA) = 0.055D0
            AC(ICLA) = 0.045D0
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
        DO IPOIN = 1,NPOIN
          DO ILAYER=1,NOMBLAY
            TOCE_MUD(ILAYER,IPOIN) = TOC_MUD(ILAYER)
          ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
