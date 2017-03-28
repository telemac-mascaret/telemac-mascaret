!                    ******************
                     SUBROUTINE LAYER3D
!                    ******************
!
     &(NLAYER,ZR,ZF,NSICLA,NPOIN)
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
!| ZF             |-->| ELEVATION OF BOTTOM
!| ZR             |-->| NON ERODABLE BED
!| NSICLA         |-->| NUMBER OF SIZE CLASSES FOR BED MATERIALS
!| NLAYER         |<--| NUMBER OF LAYER FOR EACH POINT
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_SISYPHE, ONLY: NOMBLAY
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER,          INTENT(IN)    :: NSICLA,NPOIN
      TYPE (BIEF_OBJ),  INTENT(IN)    :: ZR,ZF
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IPOIN,ILAYER,ISAND,IMUD
      DOUBLE PRECISION TOT,MIN_SED_MASS_COMP
!
      ! ADD TO POINTE
      MASS_MUD(IPOIN,NLAYER,NMUD)
      MASS_SAND(IPOIN,NLAYER,NSAND)
!
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: MASS_MUD_TOT
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: MASS_SAND_TOT
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: MASS_MIX_TOT
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: RATIO_MUD_SAND
      DOUBLE PRECISION,DIMENSION(:,:,:), ALLOCATABLE :: RATIO_SAND
      DOUBLE PRECISION,DIMENSION(:,:,:), ALLOCATABLE :: RATIO_MUD
!      MASS_MUD_TOT(IPOIN,ICOUCH)
!      MASS_SAND_TOT(IPOIN,ICOUCH)
!      MASS_MIX_TOT(IPOIN,ICOUCH)
!      RATIO_SAND(IPOIN,NLAYER,ISAND)
!      RATIO_MUD(IPOIN,NLAYER,ISAND)
!      RATIO_MUD_SAND
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      MIN_SED_MASS_COMP = 1.D-9
      ! INITIALIZATIONS AND ALLOCATIONS
      ALLOCATE(MASS_MUD_TOT(NLAYER,NPOIN))
      ALLOCATE(MASS_SAND_TOT(NLAYER,NPOIN))
      ALLOCATE(MASS_MIX_TOT(NLAYER,NPOIN))
      ALLOCATE(RATIO_SAND(NSAND,NLAYER,NPOIN))
      ALLOCATE(RATIO_MUD(NMUD,NLAYER,NPOIN))
      ALLOCATE(RATIO_MUD_SAND(NLAYER,NPOIN))
!
      DO IPOIN = 1,NPOIN
!
        DO ILAYER = 1,NLAYER
          MASS_MIX_TOT(ILAYER,IPOIN) = 0.D0
          MASS_SAND_TOT(ILAYER,IPOIN) = 0.D0
          MASS_MUD_TOT(ILAYER,IPOIN) = 0.D0
!
          DO IMUD = 1,NMUD
            ! THIS CLIPPING IS - A PRIORI - NOT MANDATORY
            IF(MASS_MUD(IMUD,ILAYER,IPOIN).LT.0.D0)THEN
              MASS_MUD(IMUD,ILAYER,IPOIN) = 0.D0
            ENDIF
          ENDDO
!
          DO ISAND = 1,NSAND
            ! THIS CLIPPING IS - A PRIORI - NOT MANDATORY
            IF(MASS_SAND(ISAND,ILAYER,IPOIN).LT.0.D0)THEN
              MASS_SAND(ISAND,ILAYER,IPOIN) = 0.D0
            ENDIF
          ENDDO
        ENDDO
!
      ENDDO
!
! COMPUTES TOT MASS PER LAYER
!
      IF(NSAND.GE.1) THEN
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NLAYER
            DO ISAND = 1,NSAND
              MASS_SAND_TOT(ILAYER,IPOIN) = MASS_SAND_TOT(ILAYER,IPOIN)
     &        + MASS_SAND(ISAND,ILAYER,IPOIN)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(NMUD.GE.1) THEN
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NLAYER
            DO IMUD = 1,NMUD
              MASS_MUD_TOT(ILAYER,IPOIN) = MASS_MUD_TOT(ILAYER,IPOIN)
     &        + MASS_MUD(IMUD,ILAYER,IPOIN)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      DO IPOIN = 1,NPOIN
        DO ILAYER = 1,NLAYER
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
          DO ILAYER = 1,NLAYER
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
          DO ILAYER = 1,NLAYER
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
          DO ILAYER=1,NLAYER
            RATIO_MUD_SAND(ILAYER,IPOIN) = 0.D0
          ENDDO
        ENDDO
      ELSEIF(NSAND.EQ.0) THEN
        DO IPOIN=1,NPOIN
          DO ILAYER=1,NLAYER
            RATIO_MUD_SAND(ILAYER,IPOIN) = 1.D0
          ENDDO
        ENDDO
      ELSE ! NSAND AND NMUD GT 0 (IF NSAND AND NMUD EQ 0: WHAT ARE YOU
           ! DOING IS SISYPHE?)
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NLAYER
            IF(MASS_SAND_TOT(ILAYER,IPOIN).GE.MIN_SED_MASS_COMP)THEN
              RATIO_MUD_SAND(ILAYER,IPOIN)= MIN(1.D0,
     &          MASS_MUD_TOT(ILAYER,IPOIN)
     &        / MASS_MIX_TOT(ILAYER,IPOIN))
            ELSE
              RATIO_MUD_SAND(ILAYER,IPOIN) = 1.D0
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
! COMPUTE THICKNESS
!
      DO ILAYER = 1,NLAYER
        DO IPOIN = 1,NPOIN
          ES_PORO_SAND(IPOIN,ILAYER) = (MASS_SAND_TOT(ILAYER,IPOIN)*XKV)
     &                               /((1-XKV)*XMVS)
          ES_MUD_ONLY(IPOIN,ILAYER) = MASS_MUD_TOT(ILAYER,IPOIN)
     &                              /CONC_VASE((IPOIN,ILAYER)
           IF(ES_MUD_ONLY(IPOIN,ILAYER).GE.ES_PORO_SAND(IPOIN,ILAYER))
     &     THEN
             ES(IPOIN,ILAYER) = MASS_SAND_TOT(ILAYER,IPOIN)/XMVS
     &                        + ES_MUD_ONLY(IPOIN,ILAYER)
           ELSE
             ES(IPOIN,ILAYER) = MASS_SAND_TOT(ILAYER,IPOIN)
     &                        /((1-XKV)*XMVS)
           ENDIF
        ENDDO
      ENDDO
! COMPUTE ZF (NOTE: THIS COULD BE MOVED)
      DO IPOIN = 1,NPOIN
        ZF(IPOIN) = ZR(IPOIN)
        DO ILAYER = 1,NLAYER
          ZF(IPOIN)= ZF(IPOIN) + ES(IPOIN,ILAYER)
        ENDDO
      ENDDO
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      ! FINALIZE
      DEALLOCATE(MASS_MUD_TOT)
      DEALLOCATE(MASS_SAND_TOT)
      DEALLOCATE(MASS_MIX_TOT)

      RETURN
      END
