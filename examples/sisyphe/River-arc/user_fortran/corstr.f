!                       *****************
                        SUBROUTINE CORSTR
!                       *****************
!
!***********************************************************************
!  TELEMAC 2D VERSION 5.2    17/08/94    J-M HERVOUET (LNH) 30 87 80 18
!
!***********************************************************************
!
!      FONCTION: CORRECTION DU COEFFICIENT DE FROTTEMENT SUR LE FOND
!                QUAND IL EST VARIABLE EN TEMPS.
!
!      CE SOUS-PROGRAMME EST SIMPLEMENT UN MODELE
!      IL DOIT ETRE REMPLI PAR L'UTILISATEUR
!
!
!
!
!-----------------------------------------------------------------------
!  EXAMPLE OF POSSIBLE ARGUMENTS
! .________________.____.______________________________________________.
! |      NOM       |MODE|                   ROLE                       |
! |________________|____|______________________________________________|
! |    CHESTR      |<-- |  COEFFICIENT DE FROTTEMENT                   |
! |    X,Y         | -->|  COORDONNEE DU MAILLAGE .                    |
! |    NPOIN       | -->|  NOMBRE DE POINTS DU MAILLAGE                |
! |    PRIVE       | -->|  TABLEAU DE TRAVAIL DEFINI DANS PRINCI       |
! |    ZF          | -->|  COTE DU FOND                                |
! |    KFROT       | -->|  LOI DE FROTTEMENT (LINEAIRE,CHEZY,STRICKLER)|
! |    FFON        | -->|  COEFFICIENT DE FROTTEMENT ASSOCIE A LA LOI  |
! |    H           | -->|  HAUTEUR D'EAU.
! |    AT          | -->|  TIME.
! |________________|____|______________________________________________|
! MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
!
!-----------------------------------------------------------------------
!
!  APPELE PAR : TELMAC
!
!  SOUS-PROGRAMME APPELE :
!
!**********************************************************************
!
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
      INTEGER I
!
      INTEGER NPMAX2,ND,NG
!
      PARAMETER (NPMAX2=200)
      DOUBLE PRECISION XD(NPMAX2),YD(NPMAX2)
      DOUBLE PRECISION XG(NPMAX2),YG(NPMAX2)
!
      DOUBLE PRECISION ZTAMPON,DISTOT,DIST
!
      DOUBLE PRECISION Q,DGRA,DMAX
      EXTERNAL         Q
!
!-----------------------------------------------------------------------
!
      ND=17
!
      XD(01)=462.665
      YD(01)=1186.187
      XD(02)=449.843
      YD(02)=1184.821
      XD(03)=446.900
      YD(03)=1181.878
      XD(04)=463.611
      YD(04)=1148.455
      XD(05)=480.638
      YD(05)=1118.921
      XD(06)=494.617
      YD(06)=1098.005
      XD(07)=504.181
      YD(07)=1082.765
      XD(08)=511.854
      YD(08)=1070.783
      XD(09)=527.199
      YD(09)=1046.504
      XD(10)=544.331
      YD(10)=1017.811
      XD(11)=548.325
      YD(11)=1011.505
      XD(12)=552.529
      YD(12)=1004.358
      XD(13)=558.625
      YD(13)= 993.847
      XD(14)= 568.715
      YD(14)= 972.091
      XD(15)= 620.952
      YD(15)= 861.942
      XD(16)=627.468
      YD(16)= 869.615
      XD(17)= 652.378
      YD(17)= 934.779
!
      NG=13
!
      XG(01)=400.759
      YG(01)=1162.329
      XG(02)=413.372
      YG(02)=1166.112
      XG(03)=431.134
      YG(03)=1129.116
      XG(04)=481.584
      YG(04)=1031.474
      XG(05)=491.043
      YG(05)=1014.868
      XG(06)=494.407
      YG(06)=1008.667
      XG(07)=518.896
      YG(07)= 966.941
      XG(08)=530.142
      YG(08)= 948.232
      XG(09)=557.574
      YG(09)= 899.884
      XG(10)=567.559
      YG(10)= 880.230
      XG(11)=573.129
      YG(11)= 866.146
      XG(12)=581.538
      YG(12)= 844.915
      XG(13)=568.400
      YG(13)= 840.606
!
!
      DO I = 1 , NPOIN
        IF (INPOLY(MESH%X%R(I),MESH%Y%R(I),XD,YD,ND)) THEN
          CHESTR%R(I) = 1.D0
!         CHESTR%R(I) = 0.3D0
        ENDIF

        IF (INPOLY(MESH%X%R(I),MESH%Y%R(I),XG,YG,NG)) THEN
          CHESTR%R(I) = 1.D0
!          CHESTR%R(I) = 0.3D0
        ENDIF

      ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END

