!            ***************************
             SUBROUTINE SUSPENSION_ERODE
!            ***************************
!
!
!***********************************************************************
! SISYPHE   V7P3                                             28/03/2017
!***********************************************************************
!
!brief    COMPUTES EROSION FLUX [kg*s-1*m-2]
!+
!
!history  R. WALTHER (ARTELIA), J. FONTAINE (EDF-LNHE)
!+        28/03/2017
!+        V7P3
!+  Creation of the subroutine.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!|                |-->|
!|FLUER           |<->| EROSION FLUX
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
      USE BIEF
      USE DECLARATIONS_SISYPHE
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!

!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IPOIN,JPOIN
      INTEGER ILAYER,ISAND,IMUD
      INTEGER NSAND_VIRTUAL
! IL FAUT DECLARER LES VARIABLES AVEC MAX(NMUD,1)et MAX(NSAND,1)
      DOUBLE PRECISION QER_MUD(NMUD),QER_SAND(NSAND),TIME(NSAND)
      DOUBLE PRECISION TOCE_MIX(NPOIN,NSAND),QE_MOY(NSAND)
      DOUBLE PRECISION FLUER_PUR_MUD(NPOIN),FLUER_PUR_SAND(NSAND,NPOIN)
      DOUBLE PRECISION FLUER_MIX(NSAND,NPOIN)
      DOUBLE PRECISION FLUER_MUD(NMUD,NPOIN),FLUER_SAND(NSAND,NPOIN)
      DOUBLE PRECISION CHECK_POSITIF_MUD_TOT1,CHECK_POSITIF_MUD_TOT2
      DOUBLE PRECISION CHECK_POSITIF_MUD,CHECK_POSITIF_SAND
      DOUBLE PRECISION FLUER_MIX_TOT
!
      IF(NSAND.EQ.0) THEN
        NSAND_VIRTUAL = 1
!       ON FAIT CELA POUR PASSER UNE FOIS DANS LA BOUCLE SABLE MEME S IL N Y EN PAS
!       CELA NE POSE PAS DE PB CAR RATIO_MUD_SAND=1
      ELSE
        NSAND_VIRTUAL = NSAND
      ENDIF
      DO IMUD = 1,NMUD
        QER_MUD(IMUD) = 0.D0
      ENDDO
      DO ISAND = 1,NSAND
        QER_SAND(ISAND) = 0.D0
        TIME(ISAND) = DT
      ENDDO
!
!     COMPUTE EQUILIBRIUM CONCENTRATION
!
      DO ISAND = 1,NSAND
        IF(NSAND.NE.0)THEN      ! TODO: FIX AGUMENTS IN COMPUTE_CAE
          CALL SUSPENSION_COMPUTE_CAE(T4,HN,FDM,FD90,NPOIN,CHARR,
     &         XMVE,XMVS,VCE,GRAV,HMIN,ZERO,XWC,ZERO,
     &         ZREF,AC,CSTAEQ,QS_C,ICQ,U2D,V2D,
     &         CSRATIO,T14,DEBUG)
! TO ADD: MULTIPLICATION XMVS, IN ORDER TO OBTAIN A FLUX IN KG*S-1*M-2
          DO JPOIN = 1,NPOIN
            CSTAEQ%R(JPOIN) = XMVS*CSTAEQ%R(JPOIN)
            FLUER_PUR_SAND(ISAND,JPOIN) =
     &           CSTAEQ%R(JPOIN)*XWC(NUM_ISAND_ICLA(ISAND))
          ENDDO
        ELSE
          DO JPOIN = 1,NPOIN
            CSTAEQ%R(JPOIN) = 0.D0
            FLUER_PUR_SAND(ISAND,JPOIN) = 0.D0
          ENDDO
        ENDIF
      ENDDO
!
      DO IPOIN = 1,NPOIN
!
        DO ILAYER = 1,NOMBLAY
!
!         COMPUTE CRITICAL SHEAR STRESS FOR MIXTURE
          DO ISAND = 1,NSAND_VIRTUAL
!           TOCE_SAND(ISAND) IS COMPUTED IN INIT_SEDIMENT
!           TODO: AJOUTER LE HIDING FACTOR SUR tOC POUR COHERENCE AVEC BEDLOAD!
!           CE N'EST PAS FAIT ACTUELLMENT FAIT DANS SYSIPHE POUR LA SUSPENSION
            IF(RATIO_MUD_SAND(ILAYER,IPOIN).LE.0.3D0)THEN
              TOCE_MIX(ISAND,IPOIN) = TOCE_SAND(ISAND)
            ELSEIF(RATIO_MUD_SAND(ILAYER,IPOIN).GE.0.5D0)THEN
              TOCE_MIX(ISAND,IPOIN) = TOCE_MUD(ILAYER,IPOIN)
            ELSE
              TOCE_MIX(ISAND,IPOIN) =
     &             (RATIO_MUD_SAND(ILAYER,IPOIN)-0.3D0)/(0.5-0.3)
     &             *(TOCE_MUD(ILAYER,IPOIN) - TOCE_SAND(ISAND))
     &             + TOCE_SAND(ISAND)
            ENDIF
          ENDDO
!
!         PUR MUD FLUX
          IF(NMUD.GT.0)THEN
! TOB(IPOIN) : A CALCULER? D'OU IL VIENT?
            IF(TOB(IPOIN).GT.TOCE_MUD(ILAYER,IPOIN))THEN
              FLUER_PUR_MUD(IPOIN) = PARTHENIADES(1) ! to expand partheniades for multilayers
     &             *(TOB(IPOIN)/TOCE_MUD(ILAYER,IPOIN) - 1.D0)
            ELSE
              FLUER_PUR_MUD(IPOIN) = 0.D0
            ENDIF
          ELSE
            FLUER_PUR_MUD(IPOIN) = 0.D0
          ENDIF
!         COMPUTE MIX FLUXES SAND-MUD
          FLUER_MIX_TOT = 0.D0
!
          DO ISAND = 1,NSAND_VIRTUAL
!
            IF(TOB(IPOIN).GT.TOCE_MIX(ILAYER,IPOIN))THEN
              IF(RATIO_MUD_SAND(ILAYER,IPOIN).LE.0.3D0)THEN
!               MUD RATIO < 30%, (PURE SAND FLUX)
                FLUER_MIX(ISAND,IPOIN)= FLUER_PUR_SAND(ISAND,IPOIN)
!               MUD RATIO > 50%, (PURE MUD FLUX)
              ELSEIF(RATIO_MUD_SAND(ILAYER,IPOIN).GE.0.5D0)THEN
                FLUER_MIX(ISAND,IPOIN) = FLUER_PUR_MUD(IPOIN)
!               MUD RATIO >30% AND <50%, (INTERPOLATION)
              ELSE
                FLUER_MIX(ISAND,IPOIN) = (RATIO_MUD_SAND(ILAYER,IPOIN)
     &               - 0.3D0)/(0.5D0-0.3D0)
     &               *(FLUER_PUR_MUD(IPOIN)-FLUER_PUR_SAND(ISAND,IPOIN))
     &               + FLUER_PUR_SAND(ISAND,IPOIN)
              ENDIF
            ELSE
              FLUER_MIX(ISAND,IPOIN) = 0.D0
            ENDIF
!
            FLUER_MIX_TOT = FLUER_MIX_TOT
     &           + FLUER_MIX(ISAND,IPOIN)*RATIO_SAND(ISAND,ILAYER,IPOIN)
!
          ENDDO
!
          IF(FLUER_MIX_TOT.LE.0.D0.AND.MASS_MIX_TOT(ILAYER,IPOIN).GE.
     &         MIN_SED_MASS_COMP) GOTO 10
!         des valeurs infimes de masse ne peuvent pas faire du pavage, on laisse la possibilite d eroder la couche en dessous
!         EXIT LAYER LOOP
          DO ISAND = 1,NSAND_VIRTUAL
!
            QE_MOY(ISAND) = FLUER_MIX(ISAND,IPOIN)*
     &           RATIO_SAND(ISAND,ILAYER,IPOIN)*TIME(ISAND)
!
            IF(QE_MOY(ISAND).GT.0.D0) THEN
!             ELSE NEXT ISAND
              IF(QE_MOY(ISAND).LT.
     &             (MASS_MUD_TOT(ILAYER,IPOIN)
     &             *RATIO_SAND(ISAND,ILAYER,IPOIN)
     &             +MASS_SAND(ISAND,ILAYER,IPOIN))) THEN
                CHECK_POSITIF_MUD_TOT1 =
     &               MIN(MASS_MUD_TOT(ILAYER,IPOIN),
     &               QE_MOY(ISAND)*RATIO_MUD_SAND(ILAYER,IPOIN))
                CHECK_POSITIF_MUD_TOT2 = 0.D0
!
                DO IMUD = 1,NMUD
                  CHECK_POSITIF_MUD =
     &                 MIN(MASS_MUD(IMUD,ILAYER,IPOIN),
     &                 CHECK_POSITIF_MUD_TOT1
     &                *RATIO_MUD(IMUD,ILAYER,IPOIN))
                  MASS_MUD(IMUD,ILAYER,IPOIN) =
     &                 MASS_MUD(IMUD,ILAYER,IPOIN) - CHECK_POSITIF_MUD
                  QER_MUD(IMUD) = QER_MUD(IMUD) + CHECK_POSITIF_MUD
                  CHECK_POSITIF_MUD_TOT2 =
     &                 CHECK_POSITIF_MUD_TOT2 + CHECK_POSITIF_MUD
                ENDDO
!
                CHECK_POSITIF_SAND = MIN(MASS_SAND(ISAND,ILAYER,IPOIN),
     &               (1.D0-RATIO_MUD_SAND(ILAYER,IPOIN))*QE_MOY(ISAND))
                QER_SAND(ISAND) = QER_SAND(ISAND) + CHECK_POSITIF_SAND
                MASS_SAND(ISAND,ILAYER,IPOIN) =
     &               MASS_SAND(ISAND,ILAYER,IPOIN) - CHECK_POSITIF_SAND
                TIME(ISAND) = 0.D0
!
              ELSE
!
                CHECK_POSITIF_MUD_TOT1 = MIN(MASS_MUD_TOT(ILAYER,IPOIN),
     &               MASS_MUD_TOT(ILAYER,IPOIN)
     &               *RATIO_SAND(ISAND,ILAYER,IPOIN))
!
                CHECK_POSITIF_MUD_TOT2 = 0.D0
                DO IMUD = 1,NMUD
                  CHECK_POSITIF_MUD =
     &                 MIN(MASS_MUD(IMUD,ILAYER,IPOIN),
     &                 CHECK_POSITIF_MUD_TOT1
     &                 *RATIO_MUD(IMUD,ILAYER,IPOIN))
!
                  MASS_MUD(IMUD,ILAYER,IPOIN)=
     &                 MASS_MUD(IMUD,ILAYER,IPOIN) - CHECK_POSITIF_MUD
!
                  QER_MUD(IMUD) = QER_MUD(IMUD) + CHECK_POSITIF_MUD
!
                  CHECK_POSITIF_MUD_TOT2 =
     &                 CHECK_POSITIF_MUD_TOT2 + CHECK_POSITIF_MUD
                ENDDO
!
                CHECK_POSITIF_SAND =
     &               MAX(MASS_SAND(ISAND,ILAYER,IPOIN),0.D0)
                MASS_SAND(ISAND,ILAYER,IPOIN) =
     &               MASS_SAND(ISAND,ILAYER,IPOIN) - CHECK_POSITIF_SAND
                QER_SAND(ISAND) = QER_SAND(ISAND) + CHECK_POSITIF_SAND
              ENDIF
!
              TIME(ISAND) = TIME(ISAND) -
     &             ((CHECK_POSITIF_MUD_TOT2 + CHECK_POSITIF_SAND)
     &             /(QE_MOY(ISAND)/DT))
!
            ENDIF
!
          ENDDO                 ! END LOOP ISAND_VIRTUAL
!
          DO ISAND = 1,NSAND_VIRTUAL
            IF(TIME(ISAND).LE.0.D0) GOTO 10
!           EXIT LAYER LOOP
          ENDDO
!
        ENDDO                   ! END LOOP NOMBLAY
10      CONTINUE
!
        DO IMUD = 1,NMUD
          FLUER_MUD(IMUD,IPOIN) = MAX(QER_MUD(IMUD)/DT,0.D0)
          FLUER(NUM_IMUD_ICLA(IMUD),IPOIN) = FLUER_MUD(IMUD,IPOIN)
        ENDDO
        DO ISAND = 1,NSAND_VIRTUAL
          FLUER_SAND(ISAND,IPOIN) = MAX(QER_SABLE(ISAND)/DT,0.D0)
          FLUER(NUM_ISAND_ICLA(ISAND),IPOIN) = FLUER_SAND(ISAND,IPOIN)
!         quand il y aura du bedload ajouter
!!!       FLUER_MUD(IMUD,IPOIN)= FLUER_MUD(IMUD,IPOIN)+FLUER_BEDLOAD_MUD(IMUD,IPOIN)
        ENDDO
!       HERE SAFETY INIT IF NO SAND (COULD BE REMOVED)
        IF(NSAND.EQ.0)THEN
          DO ILAYER = 1,NOMBLAY
            FLUER_SAND(1,IPOIN) = 0.D0
            MASS_SAND(1,ILAYER,IPOIN) = 0.D0
          ENDDO
        ENDIF
!
      ENDDO                     ! END LOOP IPOIN
!
      RETURN
      END SUBROUTINE SUSPENSION_ERODE
