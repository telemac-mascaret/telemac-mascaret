!                    *****************
                     SUBROUTINE CORFON
!                    *****************
!
!
!***********************************************************************
! TELEMAC2D   V6P1                                   21/08/2010
!***********************************************************************
!
!brief    MODIFIES THE BOTTOM TOPOGRAPHY.
!
!warning  USER SUBROUTINE
!
!history  J-M HERVOUET (LNHE)
!+        01/03/1990
!+        V5P2
!+
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      LOGICAL MAS
!
CRWR
      DOUBLE PRECISION DIST
      INTEGER I
      DOUBLE PRECISION AMP_BOSSE,AMP_PENTE_LATERALE
CRWR
!-----------------------------------------------------------------------
!
CRWR
      AMP_BOSSE = 3.5D0
      AMP_PENTE_LATERALE = 10.D0
C
      DO I=1,NPOIN
C***********bosse
        IF(Y(I).LE.0.D0) THEN
         IF(X(I).GE.-300.D0.AND.X(I).LE.-100.D0) THEN
           ZF%R(I)= ZF%R(I)+((X(I)-(-300.D0))/200.D0)*AMP_BOSSE
         ENDIF
         IF(X(I).GE.-100.D0.AND.X(I).LE.100.D0) THEN
           ZF%R(I)= ZF%R(I)+AMP_BOSSE
         ENDIF
         IF(X(I).GE.100.D0.AND.X(I).LE.300.D0) THEN
           ZF%R(I)= (ZF%R(I)+AMP_BOSSE)-
     &             ((X(I)-(100.D0))/200.D0)*AMP_BOSSE
         ENDIF
        ENDIF
C***********etranglement pour banc decouvrant
        IF(Y(I).GT.0.D0.AND.Y(I).LE.700.D0) THEN
            IF(X(I).GE.-300.D0.AND.X(I).LE.-100.D0) THEN
              ZF%R(I)= ((X(I)-(-300.D0))/200.D0)
     &        *ABS(Y(I)-700.D0)/200.D0*AMP_PENTE_LATERALE
            ENDIF
            IF(X(I).GE.-100.D0.AND.X(I).LE.100.D0) THEN
              ZF%R(I)= ABS(Y(I)-700.D0)/200.D0
     &        *AMP_PENTE_LATERALE
            ENDIF
            IF (X(I).GE.100.D0.AND.X(I).LE.300.D0) THEN
              ZF%R(I)= (1.D0-(X(I)-(100.D0))/200.D0)*
     &        ABS(Y(I)-700.D0)/200.D0*AMP_PENTE_LATERALE
            ENDIF
        ELSEIF(Y(I).GT.800.D0) THEN
            IF(X(I).GE.-300.D0.AND.X(I).LE.-100.D0) THEN
              ZF%R(I)=((X(I)-(-300.D0))/200.D0)*
     &        ABS(Y(I)-800.D0)/200.D0*AMP_PENTE_LATERALE
            ENDIF
            IF (X(I).GE.-100.D0.AND.X(I).LE.100.D0) THEN
              ZF%R(I)= ABS(Y(I)-800.D0)/200.D0
     &        *AMP_PENTE_LATERALE
            ENDIF
            IF (X(I).GE.100.D0.AND.X(I).LE.300.D0) THEN
              ZF%R(I)= (1.D0-(X(I)-(100.D0))/200.D0)*
     &        ABS(Y(I)-800.D0)/200.D0*AMP_PENTE_LATERALE
            ENDIF
        ENDIF
      ENDDO
CRWR
!  SMOOTHING(S) OF THE BOTTOM (OPTIONAL)
!
      IF(LISFON.GT.0) THEN
!
        MAS=.TRUE.
        CALL FILTER(ZF,MAS,T1,T2,AM1,'MATMAS          ',
     &              1.D0,T1,T1,T1,T1,T1,T1,MESH,MSK,MASKEL,LISFON)
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(LNG.EQ.1) THEN
        IF(LISFON.EQ.0) THEN
          WRITE(LU,*)
          WRITE(LU,*) 'CORFON (TELEMAC2D) : PAS DE MODIFICATION DU FOND'
          WRITE(LU,*)
        ELSE
          WRITE(LU,*)
          WRITE(LU,*) 'CORFON (TELEMAC2D) : ',LISFON,' LISSAGES DU FOND'
          WRITE(LU,*)
        ENDIF
      ENDIF
      IF(LNG.EQ.2) THEN
        IF(LISFON.EQ.0) THEN
          WRITE(LU,*)
          WRITE(LU,*) 'CORFON (TELEMAC2D): NO MODIFICATION OF BOTTOM'
          WRITE(LU,*)
        ELSE
          WRITE(LU,*)
          WRITE(LU,*) 'CORFON (TELEMAC2D): ',LISFON,' BOTTOM SMOOTHINGS'
          WRITE(LU,*)
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
!                    ****************
                     SUBROUTINE METEO
!                    ****************
!
     &(PATMOS,WINDX,WINDY,FUAIR,FVAIR,X,Y,AT,LT,NPOIN,VENT,ATMOS,
     & HN,TRA01,GRAV,ROEAU,NORD,PRIVE,ATMFILEA,ATMFILEB,FILES,LISTIN,
     & PATMOS_VALUE,AWATER_QUALITY,PLUIE,AOPTWIND,AWIND_SPD)
!
!***********************************************************************
! TELEMAC2D   V7P1
!***********************************************************************
!
!brief    COMPUTES ATMOSPHERIC PRESSURE AND WIND VELOCITY FIELDS
!+               (IN GENERAL FROM INPUT DATA FILES).
!
!warning  CAN BE ADAPTED BY USER
!
!history  J-M HERVOUET (LNHE)
!+        02/01/2004
!+        V5P4
!+
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
!history  J-M HERVOUET (EDF R&D, LNHE)
!+        30/01/2013
!+        V6P3
!+   Now 2 options with an example for reading a file. Extra arguments.
!
!history  C.-T. PHAM (LNHE)
!+        09/07/2014
!+        V7P0
!+   Reading a file of meteo data for exchange with atmosphere
!+   Only the wind is used here
!
!history R.ATA (LNHE)
!+        09/11/2014
!+        V7P0
!+  introducion of water quality option + pluie is introduced as
!+   an optional parameter + remove of my_option which is replaced
!+   by a new keyword + value of patmos managed also with a new keyword
!
!history  J-M HERVOUET (EDF R&D, LNHE)
!+        07/01/2015
!+        V7P0
!+  Adding optional arguments to remove USE DECLARATIONS_TELEMAC2D.
!
!history R.ATA (LNHE)
!+        16/11/2015
!+        V7P0
!+  Adding USE WAQTEL...
!
!history A. LEROY (LNHE)
!+        25/11/2015
!+        V7P1
!+  INTERPMETEO now writes directly in variables of WAQTEL which
!+  can be used by the other modules. This makes it possible to
!+  remove subsequent calls to INTERPMETEO in TELEMAC3D
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AT             |-->| TIME
!| ATMOS          |-->| YES IF PRESSURE TAKEN INTO ACCOUNT
!| FILES          |-->| BIEF_FILES STRUCTURES OF ALL FILES
!| FUAIR          |<->| VELOCITY OF WIND ALONG X, IF CONSTANT
!| FVAIR          |<->| VELOCITY OF WIND ALONG Y, IF CONSTANT
!| GRAV           |-->| GRAVITY ACCELERATION
!| HN             |-->| DEPTH
!| LISTIN         |-->| IF YES, PRINTS INFORMATION
!| LT             |-->| ITERATION NUMBER
!| NORD           |-->| DIRECTION OF NORTH, COUNTER-CLOCK-WISE
!|                |   | STARTING FROM VERTICAL AXIS
!| NPOIN          |-->| NUMBER OF POINTS IN THE MESH
!| PATMOS         |<--| ATMOSPHERIC PRESSURE
!| PATMOS_VALUE   |-->| VALUE OF ATMOSPHERIC PRESSURE IS CONSTANT
!| PRIVE          |-->| USER WORKING ARRAYS (BIEF_OBJ BLOCK)
!| ROEAU          |-->| WATER DENSITY
!| TRA01          |-->| WORKING ARRAY
!| VENT           |-->| YES IF WIND TAKEN INTO ACCOUNT
!| WINDX          |<--| FIRST COMPONENT OF WIND VELOCITY
!| WINDY          |<--| SECOND COMPONENT OF WIND VELOCITY
!| X              |-->| ABSCISSAE OF POINTS
!| Y              |-->| ORDINATES OF POINTS
!| ATMFILEA       |-->| LOGICAL UNIT OF THE ASCII ATMOSPHERIC FILE
!| ATMFILEB       |-->| LOGICAL UNIT OF THE BINARY ATMOSPHERIC FILE
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_WAQTEL,ONLY: PVAP,RAY3,NWIND,NEBU,TAIR,
     &                              TAIR_VALUE,HREL,RAINFALL,
     &                              EVAPORATION,ATMOSEXCH
      USE DECLARATIONS_TELEMAC2D, ONLY : AT1_METEO,AT2_METEO,
     &                                   FUAIR1_METEO,FUAIR2_METEO,
     &                                   FVAIR1_METEO,FVAIR2_METEO
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)             :: LT,NPOIN,ATMFILEA,ATMFILEB
      LOGICAL, INTENT(IN)             :: ATMOS,VENT,LISTIN
      DOUBLE PRECISION, INTENT(IN)    :: X(NPOIN),Y(NPOIN),HN(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: WINDX(NPOIN),WINDY(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: PATMOS(NPOIN),TRA01(NPOIN)
      DOUBLE PRECISION, INTENT(IN)    :: AT,GRAV,ROEAU,NORD,PATMOS_VALUE
      DOUBLE PRECISION, INTENT(INOUT) :: FUAIR,FVAIR
      TYPE(BIEF_OBJ), INTENT(INOUT)   :: PRIVE
      TYPE(BIEF_FILE), INTENT(IN)     :: FILES(*)
!     OPTIONAL
      LOGICAL, INTENT(IN)          ,OPTIONAL :: AWATER_QUALITY
      TYPE(BIEF_OBJ), INTENT(INOUT),OPTIONAL :: PLUIE
      INTEGER, INTENT(IN)          ,OPTIONAL :: AOPTWIND
      DOUBLE PRECISION, INTENT(IN) ,OPTIONAL :: AWIND_SPD(2)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      LOGICAL WATER_QUALITY
      INTEGER UL,OPTWIND
      DOUBLE PRECISION COEF
      DOUBLE PRECISION UAIR,VAIR,WIND_SPD(2)
!     EXCHANGE WITH ATMOSPHERE
      DOUBLE PRECISION PATM,WW,PI,TA
!
CRWR
      INTEGER I
      DOUBLE PRECISION COEF_TEMPS
CRWR
!-----------------------------------------------------------------------
!
!     DEFAULT VALUES OF PARAMETERS WHEN THEY ARE NOT GIVEN
!
      WATER_QUALITY=.FALSE.
      IF(PRESENT(AWATER_QUALITY)) WATER_QUALITY=AWATER_QUALITY
      OPTWIND=1
      IF(PRESENT(AOPTWIND)) OPTWIND=AOPTWIND
      WIND_SPD(1)=0.D0
      WIND_SPD(2)=0.D0
      IF(PRESENT(AWIND_SPD)) THEN
        WIND_SPD(1)=AWIND_SPD(1)
        WIND_SPD(2)=AWIND_SPD(2)
      ENDIF
!
!-----------------------------------------------------------------------
!
!     AT FIRST TIMESTEP
!
      IF(LT.EQ.0) THEN
!
        UL = FILES(ATMFILEA)%LU
!> SEB @ HRW, JR @ RWTH: ALGORITHMIC DIFFERENTIATION
        PI = 4.D0 * ATAN( 1.D0 )
!        PI = ACOS(-1.D0) ! ACOS NOT DIFFERENTIABLE
!< SEB @ HRW, JR @ RWTH
!
!       ATMOSPHERIC PRESSURE AND AIR TEMPERATURE
!
        IF(ATMOS.OR.WATER_QUALITY) THEN
          CALL OV( 'X=C     ' , PATMOS,Y,Y,PATMOS_VALUE,NPOIN )
        ENDIF
        IF(WATER_QUALITY) THEN
          CALL OV( 'X=C     ' , TAIR%R,Y,Y,TAIR_VALUE,NPOIN )
        ENDIF
!
!       WIND :
!
        IF(VENT.OR.WATER_QUALITY) THEN
          IF(OPTWIND.EQ.1)THEN
!           IN THIS CASE THE WIND IS CONSTANT, VALUE GIVEN IN STEERING FILE.
            CALL OV( 'X=C     ' ,WINDX,WINDX,WINDX, FUAIR , NPOIN )
            CALL OV( 'X=C     ' ,WINDY,WINDY,WINDY, FVAIR , NPOIN )
          ELSEIF(OPTWIND.EQ.2) THEN
CRWR
            CALL OV( 'X=C     ' ,WINDX,WINDX,WINDX, FUAIR , NPOIN )
            CALL OV( 'X=C     ' ,WINDY,WINDY,WINDY, FVAIR , NPOIN )
CRWR            IF(FILES(ATMFILEA)%NAME(1:1).NE.' ') THEN
!             JUMPING TWO LINES OF COMMENTS
CRWR              FUAIR = UAIR
CRWR              READ(UL,*)
CRWR              READ(UL,*)
!             READING THE FIRST TWO LINES OF DATA
CRWR              READ(UL,*) AT1_METEO,FUAIR1_METEO,FVAIR1_METEO
CRWR              IF(AT.LT.AT1_METEO) THEN
CRWR                WRITE(LU,*) ' '
CRWR                WRITE(LU,*) 'METEO'
CRWR                IF(LNG.EQ.1) WRITE(LU,*) 'DEBUT TARDIF DU FICHIER DE
CRWR     &                                    VENT'
CRWR                IF(LNG.EQ.2) WRITE(LU,*) 'LATE BEGINNING OF THE WIND
CRWR     &                                    FILE'
CRWR                CALL PLANTE(1)
CRWR                STOP
CRWR              ENDIF
CRWR            ENDIF
          ENDIF
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!
!     FOR THE REMAINING TIME STEPS
!
      IF(VENT.OR.WATER_QUALITY) THEN
!
!       WATER QUALITY
!
        IF(FILES(ATMFILEA)%NAME(1:1).NE.' ')THEN

          IF(WATER_QUALITY) THEN
!         TIME VARYING WATER QUALITY
            IF(ATMOSEXCH.EQ.0)THEN
              CALL INTERPMETEO2(NWIND,UAIR,VAIR,TA,PATM,NEBU,RAINFALL,
     &                          PVAP,RAY3,AT,UL)
!
              CALL OV('X=C     ',WINDX,WINDX,WINDX,UAIR,NPOIN)
              CALL OV('X=C     ',WINDY,WINDY,WINDY,VAIR,NPOIN)
              CALL OV('X=C     ',PATMOS,PATMOS,PATMOS,PATM,NPOIN)
              CALL OV('X=C     ',TAIR%R,TAIR%R,TAIR%R,TA  ,NPOIN)
              IF(PRESENT(PLUIE))THEN
                CALL OS('X=C     ',X = PLUIE, C=RAINFALL) ! MM/S
              ENDIF
!
!         TIME VARYING WATER QUALITY WITH HEAT EXCHANGE WITH ATMOSPHERE
            ELSEIF(ATMOSEXCH.EQ.1.OR.ATMOSEXCH.EQ.2) THEN
              CALL INTERPMETEO(WW,UAIR,VAIR,TA,PATM,
     &                         HREL,NEBU,RAINFALL,EVAPORATION,AT,UL)
              CALL OV('X=C     ',WINDX,Y,Y,UAIR,NPOIN)
              CALL OV('X=C     ',WINDY,Y,Y,VAIR,NPOIN)
!
              CALL OV('X=C     ',PATMOS,Y,Y,PATM,NPOIN)
              CALL OV('X=C     ',TAIR%R,Y,Y,TA  ,NPOIN)
            ENDIF
!
          ELSEIF (VENT) THEN
!
!           WIND VARYING IN TIME CONSTANT IN SPACE

            IF(OPTWIND.EQ.2)THEN
CRWR10            CONTINUE
CRWR              IF(AT.GE.AT1_METEO.AND.AT.LT.AT2_METEO) THEN
CRWR                IF(AT2_METEO-AT1_METEO.GT.1.D-6) THEN
CRWR                  COEF=(AT-AT1_METEO)/(AT2_METEO-AT1_METEO)
CRWR                ELSE
CRWR                  COEF=0.D0
CRWR                ENDIF
CRWR                UAIR=FUAIR1_METEO+COEF*(FUAIR2_METEO-FUAIR1_METEO)
CRWR                VAIR=FVAIR1_METEO+COEF*(FVAIR2_METEO-FVAIR1_METEO)
CRWR                IF(LISTIN) THEN
CRWR                  IF(LNG.EQ.1) WRITE(LU,*) 'VENT A T=',AT,' UAIR=',UAIR,
CRWR     &                                     ' VAIR=',VAIR
CRWR                  IF(LNG.EQ.2) WRITE(LU,*) 'WIND AT T=',AT,' UAIR=',
CRWR     &                                      UAIR,' VAIR=',VAIR
CRWR                ENDIF
CRWR              ELSE
CRWR                AT1_METEO=AT2_METEO
CRWR                FUAIR1_METEO=FUAIR2_METEO
CRWR                FVAIR1_METEO=FVAIR2_METEO
CRWR                READ(UL,*,ERR=100,END=200) AT2_METEO,FUAIR2_METEO,
CRWR     &                                     FVAIR2_METEO
CRWR                GO TO 10
!
!-----------------------------------------------------------------------
!
CRWR100             CONTINUE
CRWR                WRITE(LU,*) ' '
CRWR                WRITE(LU,*) 'METEO'
CRWR                IF(LNG.EQ.1) WRITE(LU,*) 'ERREUR DANS LE FICHIER DE
CRWR     &                                    VENT'
CRWR                IF(LNG.EQ.2) WRITE(LU,*) 'ERROR IN THE WIND FILE'
CRWR                CALL PLANTE(1)
CRWR                STOP
CRWR200             CONTINUE
CRWR                WRITE(LU,*) ' '
CRWR                WRITE(LU,*) 'METEO'
CRWR                IF(LNG.EQ.1)WRITE(LU,*)'FIN PREMATUREE DU FICHIER DE
CRWR     &                                  VENT'
CRWR                IF(LNG.EQ.2)WRITE(LU,*) 'WIND FILE TOO SHORT'
CRWR                CALL PLANTE(1)
CRWR                STOP
CRWR              ENDIF
CRWR              CALL OV('X=C     ',WINDX,Y,Y,UAIR,NPOIN)
CRWR              CALL OV('X=C     ',WINDY,Y,Y,VAIR,NPOIN)
!
CRWR              FUAIR = UAIR
CRWR              FVAIR = VAIR
!
!             WIND VARYING IN TIME AND SPACE
!
            ELSEIF(OPTWIND.EQ.3)THEN
              IF(LNG.EQ.1) THEN
                WRITE(LU,*) 'CETTE OPTION N EST PAS ENCORE PROGRAMMEE'
                WRITE(LU,*) 'VOIR CAS DE VALIDATION WIND_TXY '
                WRITE(LU,*) 'DANS LE DOSSIER EXAMPLES/TELEMAC2D'
              ELSE
                WRITE(LU,*) 'THIS OPTION IS NOT IMPLEMENTED YET'
                WRITE(LU,*) 'SEE VALIDATION CASE WIND_TXY '
                WRITE(LU,*) 'LOCATED AT THE FOLDER EXAMPLES/TELEMAC2D'
              ENDIF
              CALL PLANTE(1)
              STOP
            ENDIF
          ENDIF
        ENDIF
!
!       WIND AND/OR WATER QUALITY VARIABLES
!       VARYING IN SPACE AND TIME, FROM A BINARY FILE
!
        IF(FILES(ATMFILEB)%NAME(1:1).NE.' ') THEN
          IF(FILES(ATMFILEA)%NAME(1:1).NE.' ')THEN
            IF(LNG.EQ.1.AND.LISTIN) THEN
              WRITE(LU,*) 'METEO : LES DONNEES DU FICHIER'
              WRITE(LU,*) 'DE METEO ASCII PRESENTES DANS LE'
              WRITE(LU,*) 'FICHIER METEO BINAIRE VONT ETRE ECRASEES'
            ENDIF
            IF(LNG.EQ.2.AND.LISTIN) THEN
              WRITE(LU,*) 'METEO: THE DATA FROM THE ASCII METEO'
              WRITE(LU,*) 'FILE WILL BE OVERWRITTEN BY THE'
              WRITE(LU,*) 'CORRESPONDING BINARY FILE DATA'
            ENDIF
          ENDIF
          CALL METEO_FROM_BINARY_FILE(PATMOS,WINDX,WINDY,AT,NPOIN,VENT,
     &                                ATMOS,ATMFILEB,FILES,LISTIN,
     &                                WATER_QUALITY,PLUIE,OPTWIND,
     &                                WIND_SPD)
        ENDIF
CRWR
        IF(AT.LE.900.D0)THEN
           COEF_TEMPS=AT/900.D0
C        ELSEIF(AT.LE.20000.D0)THEN
C           COEF_TEMPS=1.D0
C        ELSE
           COEF_TEMPS=1.D0
        ENDIF
CRWR        WRITE(LU,*)'COEF_TEMPS',COEF_TEMPS
        DO I=1,NPOIN
CRWR           WINDX(I)=(-Y(I)/5000.D0*500.D0)*COEF_TEMPS
CRWR           WINDY(I)=(X(I)/5000.D0*500.D0)*COEF_TEMPS
           WINDY(I)=0.D0
           IF(ABS(X(I)).LE.1000.D0.AND.Y(I).LE.0.D0)THEN
               WINDX(I)=30.D0*COEF_TEMPS
           ENDIF
           IF(ABS(X(I)).LE.1000.D0.AND.Y(I).GE.0.D0)THEN
               WINDX(I)=-30.D0*COEF_TEMPS
           ENDIF
        ENDDO
CRWR
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
C***********************************************************************
C***********************************************************************
C*****************************SISYPHE***********************************
C***********************************************************************
C***********************************************************************
C***********************************************************************
C
!                    *****************
                     SUBROUTINE NOEROD
!                    *****************
!
     & (H , ZF , ZR , Z , X , Y , NPOIN , CHOIX , NLISS )
!
!***********************************************************************
! SISYPHE   V6P3                                  21/07/2011
!***********************************************************************
!
!brief    FIXES THE NON-ERODABLE BED ELEVATION ZR.
!
!note     METHODS OF TREATMENT OF NON-ERODABLE BEDS CAN LEAD TO ZF.
!note  CHOOSE TO SMOOTH THE SOLUTION WITH NLISS > 0.
!
!history  C. LENORMANT
!+
!+        V5P1
!+
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
!history  J-M HERVOUET (EDF R&D, LNHE)
!+        21/06/2013
!+        V6P3
!+   Now ZR=ZF-100.D0 by default
!+   previous versions was erronneously ZR=-100.D0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| CHOIX          |-->| SELECTED METHOD FOR THE TREATMENT OF RIGID BEDS
!| H              |-->| WATER DEPTH
!| NLISS          |<->| NUMBER OF SMOOTHINGS
!| NPOIN          |-->| NUMBER OF 2D POINTS
!| X,Y            |-->| 2D COORDINATES
!| Z              |-->| FREE SURFACE
!| ZF             |-->| BED LEVEL
!| ZR             |<--| RIGID BED LEVEL
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN):: NPOIN , CHOIX
      INTEGER, INTENT(INOUT):: NLISS
!
      DOUBLE PRECISION, INTENT(IN)::  Z(NPOIN) , ZF(NPOIN)
      DOUBLE PRECISION , INTENT(IN)::  X(NPOIN) , Y(NPOIN), H(NPOIN)
      DOUBLE PRECISION , INTENT(INOUT)::  ZR(NPOIN)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER I
!
!---------------------
! RIGID BEDS POSITION
!---------------------
!
!     DEFAULT VALUE: ZR=ZF-100.D0
!
      CALL OV('X=Y+C   ',ZR,ZF,ZF,-0.D0,NPOIN)
      DO I=1,NPOIN
C***********bosse 1
         IF(X(I).GE.-500.D0.AND.X(I).LE.500.D0)THEN
            ZR(I)= ZF(I)-0.2D0
         ENDIF
      ENDDO
!
!------------------
! SMOOTHING OPTION
!------------------
!
!     NLISS : NUMBER OF SMOOTHING IF  (ZF - ZR ) NEGATIVE
!             DEFAULT VALUE : NLISS = 0 (NO SMOOTHING)
!
      NLISS = 0
!
!--------------------------------------------------
! CONTROL (CAN BE ACTIVATED IF ZR USER DEFINED...)
!--------------------------------------------------
!
!     DO I=1,NPOIN
!       IF(ZR(I).GT.ZF(I)) THEN
!         WRITE(LU,*) 'POINT ',I,' NON ERODABLE BED HIGHER THAN BED'
!         CALL PLANTE(1)
!         STOP
!       ENDIF
!     ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END
