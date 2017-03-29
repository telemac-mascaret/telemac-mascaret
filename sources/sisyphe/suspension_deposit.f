!                    *****************************
                     SUBROUTINE SUSPENSION_DEPOSIT
!                    *****************************
!
     &(CA)
!
!***********************************************************************
! SISYPHE   V7P3                                             28/03/2017
!***********************************************************************
!
!brief    COMPUTES FIRST LAYER DEPOSITION;
!+
!
!history  R. WALTHER (ARTELIA), J. FONTAINE (EDF-LNHE)
!+        28/03/2017
!+        V7P3
!+  Creation of the subroutine.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| CA             |-->| BOTTOM CONCENTRATION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_SISYPHE
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      TYPE(BIEF_OBJ),   INTENT(IN) :: CA
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER IPOIN,IMUD,ISAND,ICLA
      INTEGER LAYER_DEPOSIT_SAND,LAYER_DEPOSIT_MUD
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      LAYER_DEPOSIT_SAND = 1
      LAYER_DEPOSIT_MUD = 1
      ISAND = 0
      IMUD = 0

      DO IPOIN = 1,NPOIN
        DO ICLA = 1,NSICLA
          IF(TYPE_OF_SEDIMENT(ICLA).EQ.SED_NCO) THEN
            ISAND = ISAND+1
            MASS_SAND(ISAND,LAYER_DEPOSIT_SAND,IPOIN) =
     &          MASS_SAND(ISAND,LAYER_DEPOSIT_SAND,IPOIN)
     &        + XWC(ICLA)*MAX(CA%ADR(ICLA)%P%R(IPOIN),0.D0)*DT
          ELSE
            IMUD = IMUD+1
            MASS_MUD(IMUD,LAYER_DEPOSIT_MUD,IPOIN) =
     &          MASS_SAND(ISAND,LAYER_DEPOSIT_MUD,IPOIN)
     &        + XWC(ICLA)*MAX(CA%ADR(ICLA)%P%R(IPOIN),0.D0)*DT
          ENDIF
        ENDDO
      ENDDO
!
      RETURN
      END
