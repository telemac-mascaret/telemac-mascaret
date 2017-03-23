!                       *************************
                        SUBROUTINE CONDIM_SISYPHE
!                       *************************
!
     & (U      , V   , QU    , QV  , H   , ZF , Z ,
     &  ESOMT  , THETAWR  , Q     , HWR  , TWR  ,
     &  X      , Y   , NPOIN , AT  , PMAREE)
!
!***********************************************************************
! SISYPHE VERSION 5.3                             E. PELTIER    11/09/95
!                                                 C. LENORMANT
!                                                 J.-M. HERVOUET
!
! COPYRIGHT EDF-DTMPL-SOGREAH-LHF-GRADIENT
!***********************************************************************
!
!     FONCTION  : VALEURS IMPOSEES
!                         - DU DEBIT VECTORIEL    QU, QV
!                         - DE LA HAUTEUR D'EAU   H
!                         - DE LA COTE DU FOND    ZF
!                         - DE LA SURFACE LIBRE   Z
!                         - DE L'EVOLUTION TOTALE ESOMT
!                         - DU DEBIT              Q
!                         - DE LA HAUTEUR DE HOULE HW
!                         - DE LA PERIODE DE HOULE TW
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________
! |      NOM       |MODE|                   ROLE
! |________________|____|______________________________________________
! |   U , V        |<-- | COORDONNEES DES VECTEURS VITESSE
! |   QU , QV      |<-- | DEBIT VECTORIEL SUIVANT X ET SUIVANT Y
! |   H            |<-->| HAUTEUR D'EAU
! |   ZF           |<-->| COTE DU FOND
! |   Z            |<-->| COTE DE SURFACE LIBRE
! |   ESOMT        |<-->| EVOLUTION TOTALE DES FONDS
! |   C            |<-->| CELERITE
! |   Q            |<-->| DEBIT
! |   HW           | -->| HAUTEUR DE HOULE
! |   TW           | -->| PERIODE DE HOULE
! |   X,Y          | -->| COORDONNEES DU MAILLAGE
! |   NPOIN        | -->| NOMBRE DE POINTS DU MAILLAGE
! |   AT           | -->| TEMPS
! |   PMAREE       | -->| PERIODE DE LA MAREE
! |________________|____|______________________________________________
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!-----------------------------------------------------------------------
! PROGRAMME APPELANT : SISYPH
! PROGRAMMES APPELES :
!***********************************************************************
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)::NPOIN
!
      DOUBLE PRECISION, INTENT(IN):: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN):: AT , PMAREE
! SEDIMENT
      DOUBLE PRECISION, INTENT(INOUT) ::  ZF(NPOIN)
      DOUBLE PRECISION, INTENT (INOUT)::  ESOMT(NPOIN)
! HYDRODYNAMICS
      DOUBLE PRECISION, INTENT(INOUT):: Z(NPOIN) , H(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT):: U(NPOIN) , V(NPOIN)
      DOUBLE PRECISION, INTENT (INOUT)::QU(NPOIN), QV(NPOIN), Q(NPOIN)
! WAVES
      DOUBLE PRECISION, INTENT (INOUT):: HWR(NPOIN) , TWR(NPOIN)
      DOUBLE PRECISION, INTENT (INOUT):: THETAWR(NPOIN)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER I
      DOUBLE PRECISION   PI
!
!-----------------------------------------------------------------------
!
!
!  --------------------------------------------------------------
!  INITIALISATION DES TABLEAUX NON LUS DANS LE FICHIER RESULTATS:
!  --------------------------------------------------------------
!
      DO I=1,NPOIN
        QU(I)=.25D0
        Q(I)=.25D0
        QV(I)=0.D0
        Z(I)=.6D0
        ZF(I)=0.D0
      ENDDO
      PI=3.1415926D0
      DO I=1,NPOIN
        IF (X(I) .GE. 2.D0 .AND. X(I) .LE. 10.D0) THEN
          ZF(I)=0.1D0*SIN(PI*(X(I)-2.D0)/8.D0)**2
        ENDIF
        H(I)=Z(I)-ZF(I)
! 25/01/2007
! Definition du U2D, V2D
!
        U(I)=QU(I)/H(I)
        V(I)=QV(I)/H(I)
      ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE CONDIM_SISYPHE

