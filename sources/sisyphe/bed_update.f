!                    *********************
                     SUBROUTINE BED_UPDATE
!                    *********************
!
     &(ZR,ZF,VOLU2D)
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
!| VOLU2D         |-->| INTEGRAL OF TEST FUNCTIONS
!| ZF             |<->| ELEVATION OF BOTTOM
!| ZR             |<->| NON ERODABLE BED
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_SISYPHE, ONLY: MASS_MIX_TOT,MASS_SAND_TOT,
     & MASS_MUD_TOT,MASS_MUD,MASS_SAND,NSAND,NMUD,RATIO_MUD,RATIO_SAND,
     & RATIO_MUD_SAND,NOMBLAY,CONC_MUD,ES,XKV,XMVS,MIN_SED_MASS_COMP,
     & NPOIN,MASSTOT,NUM_ICLA_IMUD,NUM_ICLA_ISAND,NSICLA,ELAY
      USE DECLARATIONS_SPECIAL
      USE INTERFACE_PARALLEL, ONLY : P_DSUM
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      TYPE (BIEF_OBJ),  INTENT(INOUT)    :: ZR,ZF
      TYPE (BIEF_OBJ),  INTENT(IN)       :: VOLU2D
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IPOIN,ILAYER,ISAND,IMUD,ICLA,K,J,I
      DOUBLE PRECISION TOT,TERM,DISCR
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!--------------------------------------------------------------------------------------
! ONLY MASS_SAND_TOT(1,IPOIN) SET TO ZERO SINCE AFTER BEDLOAD MASS_SAND (BUT ONLY THE
! FIRST LAYER) IS UPDATED IN BEDLOAD_MAIN
! to do for all the layers if masses are updated/computed somewhere before
! VARIOUS MASSES ARE STILL IN [kg/m2]
      DO IPOIN = 1,NPOIN
!!!
!!        DO ILAYER = 1,NOMBLAY
!!          MASS_MIX_TOT(ILAYER,IPOIN) = 0.D0
          MASS_SAND_TOT(1,IPOIN) = 0.D0
!!          MASS_MUD_TOT(ILAYER,IPOIN) = 0.D0
!!!
!!          IF(NMUD.NE.0)THEN
!!            DO IMUD = 1,NMUD
!!              ! THIS CLIPPING IS - A PRIORI - NOT MANDATORY
!!              IF(MASS_MUD(IMUD,ILAYER,IPOIN).LT.0.D0)THEN
!!                MASS_MUD(IMUD,ILAYER,IPOIN) = 0.D0
!!              ENDIF
!!            ENDDO
!!          ELSE
!!! IL FAUT DECLARER MASS_MUD avec MAX(NMUD,1)
!!            MASS_MUD(1,ILAYER,IPOIN) = 0.D0
!!          ENDIF
!!!
!!          IF(NSAND.NE.0)THEN
!!            DO ISAND = 1,NSAND
!!              ! THIS CLIPPING IS - A PRIORI - NOT MANDATORY
!!              IF(MASS_SAND(ISAND,ILAYER,IPOIN).LT.0.D0)THEN
!!                MASS_SAND(ISAND,ILAYER,IPOIN) = 0.D0
!!              ENDIF
!!            ENDDO
!!          ELSE
!!! IL FAUT DECLARER MASS_SAND avec  MAX(NSAND,1)
!!            MASS_SAND(1,ILAYER,IPOIN) = 0.D0
!!          ENDIF
!!        ENDDO
!!!
      ENDDO
! end attention
!-----------------------------------------------------------------------------
! UPDATES TOT MASS PER LAYER
! here we receive MASS_SAND(ISAND,ILAYER,IPOIN) or MASS_MUD(IMUD,ILAYER,IPOIN)
! after bedload/suspension/...
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
              IF(ISAND.NE.NSAND) THEN
                IF(MASS_SAND_TOT(ILAYER,IPOIN).GE.MIN_SED_MASS_COMP)THEN
                  RATIO_SAND(ISAND,ILAYER,IPOIN) =
     &            MIN(1.D0,MASS_SAND(ISAND,ILAYER,IPOIN)
     &            / MASS_SAND_TOT(ILAYER,IPOIN))
                  TOT = TOT + RATIO_SAND(ISAND,ILAYER,IPOIN)
                ELSE
                  RATIO_SAND(ISAND,ILAYER,IPOIN) = 0.D0
                ENDIF
              ELSE
                RATIO_SAND(NSAND,ILAYER,IPOIN) = 1.D0-TOT
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
!       TERM REPRESENTS THE DIFFERENCE BETWEEN MUD VOLUME AND VOID VOLUME
          TERM=RATIO_MUD_SAND(ILAYER,IPOIN)/CONC_MUD(ILAYER)-
     &    (XKV(ILAYER)*(1.D0-RATIO_MUD_SAND(ILAYER,IPOIN)))/
     &    (XMVS*(1.D0-XKV(ILAYER)))
          DISCR=MAX(0.D0,TERM)
!       IF DISCR IS POSITIVE IT MEANS THAT MUD VOLUME IS LARGER THAN VOID VOLUME
!       IF DISCR IS NEGATIVE, THE VOID VOLUME IS NOT COMPLETELY FILLED BY MUD
          ES(IPOIN,ILAYER)=MASS_MIX_TOT(ILAYER,IPOIN)*
     &    ((1.D0-RATIO_MUD_SAND(ILAYER,IPOIN))/(XMVS*(1.D0-XKV(ILAYER)))
     &    + DISCR)
        ENDDO
      ENDDO
! THICKNESS OF THE FIRST LAYER IS THICKNESS OF ACTIVE LAYER
!
      IF(NSAND.GT.1.AND.NOMBLAY.GT.1) THEN
        DO I=1,NPOIN
          ELAY%R(I)=ES(I,1)
! OPTION ACTIVE LAYER NOT CONSTANT: TO DO!
!            ELAY%R(I)=3.D0 * ACLADM%R(J)
        ENDDO
      ENDIF
!
! COMPUTE ZF (NOTE: THIS COULD BE MOVED)
!
      DO IPOIN = 1,NPOIN
        ZF%R(IPOIN) = ZR%R(IPOIN)
        DO ILAYER = 1,NOMBLAY
          ZF%R(IPOIN) = ZF%R(IPOIN) + ES(IPOIN,ILAYER)
        ENDDO
      ENDDO
!
! COMPUTES MASSTOT AT EVERY TIME STEP, FOR EVERY CLASS
! NECESSARY FOR BALANCE
!
      DO ICLA=1,NSICLA
        MASSTOT(ICLA)=0.D0
      ENDDO
!
      DO ICLA=1,NSICLA
        K = NUM_ICLA_IMUD(ICLA)
        J = NUM_ICLA_ISAND(ICLA)
        IF(J.GT.0) THEN
          DO IPOIN=1,NPOIN
            DO I=1,NOMBLAY
              MASSTOT(ICLA)= MASSTOT(ICLA) +
     &                   MASS_SAND(J,I,IPOIN)*VOLU2D%R(IPOIN)
            ENDDO
          ENDDO
        ENDIF
        IF(K.GT.0) THEN
          DO IPOIN=1,NPOIN
            DO I=1,NOMBLAY
              MASSTOT(ICLA)= MASSTOT(ICLA) +
     &                   MASS_MUD(K,I,IPOIN)*VOLU2D%R(IPOIN)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
!
      IF(NCSIZE.GT.1) THEN
        DO I=1,NSICLA
          MASSTOT(I)=P_DSUM(MASSTOT(I))
        ENDDO
      ENDIF
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      RETURN
      END