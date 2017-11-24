!                    *********************
                     SUBROUTINE BED_UPDATE
!                    *********************
!
     &(ZR_T3D,ZF_T3D)
!
!***********************************************************************
! SISYPHE   V7P3                                             28/03/2017
!***********************************************************************
!
!brief    COMPUTES BED EVOlUTION;
!+
!+            ACTIVE LAYER IS LAYER 1, IT IS KEPT AT A PRESCRIBED
!
!history  R. WALTHER (ARTELIA), J. FONTAINE (EDF-LNHE)
!+        28/03/2017
!+        V7P3
!+  Creation of the subroutine.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| ZF_T3D             |-->| ELEVATION OF BOTTOM
!| ZR_T3D             |-->| NON ERODABLE BED
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_SISYPHE
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      TYPE (BIEF_OBJ),  INTENT(INOUT)    :: ZR_T3D,ZF_T3D
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IPOIN,ILAYER,ISAND,IMUD
      DOUBLE PRECISION TOT
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! various masses are in [kg/m2]
      DO IPOIN = 1,NPOIN
!
        DO ILAYER = 1,NOMBLAY
          MASS_MIX_TOT(ILAYER,IPOIN) = 0.D0
          MASS_SAND_TOT(ILAYER,IPOIN) = 0.D0
          MASS_MUD_TOT(ILAYER,IPOIN) = 0.D0
!
          IF(NMUD.NE.0)THEN
            DO IMUD = 1,NMUD
              ! THIS CLIPPING IS - A PRIORI - NOT MANDATORY
              IF(MASS_MUD(IMUD,ILAYER,IPOIN).LT.0.D0)THEN
                MASS_MUD(IMUD,ILAYER,IPOIN) = 0.D0
              ENDIF
            ENDDO
          ELSE
! IL FAUT DECLARER MASS_MUD avec  MAX(NMUD,1)
            MASS_MUD(1,ILAYER,IPOIN) = 0.D0
          ENDIF
!
          IF(NSAND.NE.0)THEN
            DO ISAND = 1,NSAND
              ! THIS CLIPPING IS - A PRIORI - NOT MANDATORY
              IF(MASS_SAND(ISAND,ILAYER,IPOIN).LT.0.D0)THEN
                MASS_SAND(ISAND,ILAYER,IPOIN) = 0.D0
              ENDIF
            ENDDO
          ELSE
! IL FAUT DECLARER MASS_SAND avec  MAX(NSAND,1)
            MASS_SAND(1,ILAYER,IPOIN) = 0.D0
          ENDIF
        ENDDO
!
      ENDDO
!
! UPDATES TOT MASS PER LAYER
! it means that after bedload here we receive MASS_SAND(ISAND,ILAYER,IPOIN)/MUD
!
      IF(NSAND.GE.1) THEN
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NOMBLAY
            DO ISAND = 1,NSAND
              MASS_SAND_TOT(ILAYER,IPOIN) = MASS_SAND_TOT(ILAYER,IPOIN)
     &        + MASS_SAND(ISAND,ILAYER,IPOIN)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(NMUD.GE.1) THEN
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NOMBLAY
            DO IMUD = 1,NMUD
              MASS_MUD_TOT(ILAYER,IPOIN) = MASS_MUD_TOT(ILAYER,IPOIN)
     &        + MASS_MUD(IMUD,ILAYER,IPOIN)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      DO IPOIN = 1,NPOIN
        DO ILAYER = 1,NOMBLAY
          MASS_MIX_TOT(ILAYER,IPOIN) = MASS_SAND_TOT(ILAYER,IPOIN)
     &    + MASS_MUD_TOT(ILAYER,IPOIN)
        ENDDO
      ENDDO
!
! COMPUTES MASS RATIOS PER LAYER (SAND CLASSES)
!
      IF(NSAND.GE.1)THEN
!
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NOMBLAY
            TOT = 0.D0
            DO ISAND = 1,NSAND
              IF (ISAND.NE.NSAND)THEN
                IF(MASS_SAND_TOT(ILAYER,IPOIN).GE.MIN_SED_MASS_COMP)THEN
                   RATIO_SAND(ISAND,ILAYER,IPOIN) =
     &             MIN(1.D0,MASS_SAND(ISAND,ILAYER,IPOIN)
     &             / MASS_SAND_TOT(ILAYER,IPOIN))
                   TOT = TOT + RATIO_SAND(ISAND,ILAYER,IPOIN)
                ELSE
                  RATIO_SAND(IPOIN,ILAYER,ISAND) = 0.D0
                ENDIF
              ELSE
                RATIO_SAND(IPOIN,ILAYER,NSAND) = 1.D0-TOT
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!
      ENDIF
!
! COMPUTES MASS RATIOS PER LAYER (MUD CLASSES)
!
      IF(NMUD.GE.1)THEN
!
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NOMBLAY
            TOT = 0.D0
            DO IMUD = 1,NMUD
              IF (IMUD.NE.NMUD)THEN
                IF(MASS_MUD_TOT(ILAYER,IPOIN).GE.MIN_SED_MASS_COMP)THEN
                   RATIO_MUD(IMUD,ILAYER,IPOIN) =
     &             MIN(1.D0,MASS_MUD(IMUD,ILAYER,IPOIN)
     &             / MASS_MUD_TOT(ILAYER,IPOIN))
                   TOT = TOT + RATIO_MUD(IMUD,ILAYER,IPOIN)
                ELSE
                  RATIO_MUD(IPOIN,ILAYER,IMUD) = 0.D0
                ENDIF
              ELSE
                RATIO_MUD(IPOIN,ILAYER,NMUD) = 1.D0-TOT
              ENDIF
            ENDDO
          ENDDO
        ENDDO
!
      ENDIF
!
! COMPUTES SAND-MUD MASS RATIOS PER LAYER
!
      IF(NMUD.EQ.0) THEN
        DO IPOIN=1,NPOIN
          DO ILAYER=1,NOMBLAY
            RATIO_MUD_SAND(ILAYER,IPOIN) = 0.D0
          ENDDO
        ENDDO
      ELSEIF(NSAND.EQ.0) THEN
        DO IPOIN=1,NPOIN
          DO ILAYER=1,NOMBLAY
            RATIO_MUD_SAND(ILAYER,IPOIN) = 1.D0
          ENDDO
        ENDDO
      ELSE ! NSAND AND NMUD GT 0 (IF NSAND AND NMUD EQ 0: WHAT ARE YOU
           ! DOING IS SISYPHE?)
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NOMBLAY
            IF(MASS_MIX_TOT(ILAYER,IPOIN).GE.MIN_SED_MASS_COMP)THEN
              RATIO_MUD_SAND(ILAYER,IPOIN)= MIN(1.D0,
     &          MASS_MUD_TOT(ILAYER,IPOIN)
     &        / MASS_MIX_TOT(ILAYER,IPOIN))
            ELSE
! CHOICE TO DO, BUT NO EROSION OF THIS LAYER IN SUSPENSION_ERODE
              RATIO_MUD_SAND(ILAYER,IPOIN) = 1.D0
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
! COMPUTE THICKNESS [m]
!
      DO ILAYER = 1,NOMBLAY
        DO IPOIN = 1,NPOIN
! ES_PORO_SAND: is the thicnkess of voids in the sand volume.
! ces deux variables sont à garder? sinon on simplifie
! poro_sand est aussi utilise dans update_active_layer_hirano : plutot que le stocker on le recalcul?
          ES_PORO_SAND(IPOIN,ILAYER) = MASS_SAND_TOT(ILAYER,IPOIN)*
     &                         XKV(ILAYER)/((1.D0-XKV(ILAYER))*XMVS)
          ES_MUD_ONLY(IPOIN,ILAYER) = MASS_MUD_TOT(ILAYER,IPOIN)
     &                         /CONC_VASE(ILAYER)
           IF(ES_MUD_ONLY(IPOIN,ILAYER).GE.ES_PORO_SAND(IPOIN,ILAYER))
     &     THEN
             ES(IPOIN,ILAYER) = MASS_SAND_TOT(ILAYER,IPOIN)/XMVS
     &                        + ES_MUD_ONLY(IPOIN,ILAYER)
           ELSE
             ES(IPOIN,ILAYER) = MASS_SAND_TOT(ILAYER,IPOIN)
     &                        /((1.D0-XKV(ILAYER))*XMVS)
           ENDIF
        ENDDO
      ENDDO
!
! COMPUTE ZF (NOTE: THIS COULD BE MOVED)
!
      DO IPOIN = 1,NPOIN
        ZF_T3D%R(IPOIN) = ZR_T3D%R(IPOIN)
        DO ILAYER = 1,NOMBLAY
          ZF_T3D%R(IPOIN) = ZF_T3D%R(IPOIN) + ES(IPOIN,ILAYER)
        ENDDO
      ENDDO
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
! FINALIZE
      DEALLOCATE(MASS_MUD_TOT)
      DEALLOCATE(MASS_SAND_TOT)
      DEALLOCATE(MASS_MIX_TOT)
      DEALLOCATE(RATIO_SAND)
      DEALLOCATE(RATIO_MUD)
      DEALLOCATE(RATIO_MUD_SAND)
      DEALLOCATE(ES_PORO_SAND)
      DEALLOCATE(ES_MUD_ONLY)
!
      RETURN
      END
