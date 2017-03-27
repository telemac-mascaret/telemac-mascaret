!                    ***********************
                     MODULE INTERFACE_STBTEL
!                    ***********************
!
!
!***********************************************************************
! TELEMAC2D 7.2
!***********************************************************************
!
!
!-----------------------------------------------------------------------
!
!     DEFINITION OF INTERFACES
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE REFINE_MESH
     &(RLEVELS,MESHINIT,NNELMAX,NPTFRMAX,NTRAC,EXTEND_LIM,
     & CORRESP,LIHBOR,LIUBOR,LIVBOR,LITBOR,HBOR,UBOR,VBOR,
     & CHBORD,TBOR,ATBOR,BTBOR,ZF,H,TEXP,TTILD,TN)
        USE BIEF
        USE DECLARATIONS_TELEMAC
        USE DECLARATIONS_STBTEL
        USE DECLARATIONS_SPECIAL
        IMPLICIT NONE
        TYPE(BIEF_MESH), INTENT(INOUT) :: MESHINIT
        INTEGER        , INTENT(IN)    :: RLEVELS
        INTEGER        , INTENT(IN)    :: NNELMAX
        INTEGER        , INTENT(IN)    :: NPTFRMAX
        INTEGER        , INTENT(IN)    :: NTRAC
        LOGICAL        , INTENT(IN)    :: EXTEND_LIM
        INTEGER,INTENT(INOUT), OPTIONAL :: CORRESP(NNELMAX,RLEVELS)
        INTEGER,INTENT(INOUT), OPTIONAL :: LIHBOR(NPTFRMAX)
        INTEGER,INTENT(INOUT), OPTIONAL :: LIUBOR(NPTFRMAX)
        INTEGER,INTENT(INOUT), OPTIONAL :: LIVBOR(NPTFRMAX)
        TYPE(BIEF_OBJ),INTENT(INOUT), OPTIONAL :: LITBOR
        DOUBLE PRECISION,INTENT(INOUT),OPTIONAL:: UBOR(NPTFRMAX,2)
        DOUBLE PRECISION,INTENT(INOUT),OPTIONAL:: VBOR(NPTFRMAX,2)
        DOUBLE PRECISION,INTENT(INOUT),OPTIONAL:: HBOR(NPTFRMAX)
        DOUBLE PRECISION,INTENT(INOUT),OPTIONAL:: CHBORD(NPTFRMAX)
        TYPE(BIEF_OBJ),INTENT(INOUT),OPTIONAL:: ZF,H
        TYPE(BIEF_OBJ),INTENT(INOUT),OPTIONAL:: TBOR, ATBOR
        TYPE(BIEF_OBJ),INTENT(INOUT),OPTIONAL:: BTBOR
        TYPE(BIEF_OBJ),INTENT(INOUT),OPTIONAL:: TEXP,TTILD,TN
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DIVISE
     &(X,Y,IKLE,NCOLOR,NPOIN,NELEM,NELMAX,NSOM2,SOM2,INDICP,INDICE,
     & SHP,ELT,NPMAX,CORR,LEVEL)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NSOM2,NELMAX,NPMAX
      INTEGER, INTENT(INOUT) :: ELT(NPMAX)
      INTEGER, INTENT(INOUT) :: NPOIN,NELEM
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,*),INDICP(*),INDICE(*)
      INTEGER, INTENT(INOUT) :: NCOLOR(*)
      DOUBLE PRECISION, INTENT(IN) :: SOM2(10,2)
      DOUBLE PRECISION, INTENT(INOUT) :: X(*),Y(*),SHP(NPMAX,3)
      INTEGER, INTENT(INOUT), OPTIONAL :: CORR(NELMAX,*)
      INTEGER, INTENT(IN), OPTIONAL :: LEVEL
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FILL_LIM
     & (NPTFR,NPTFRX,NTRAC,LIHBOR,LIUBOR,LIVBOR,LITBOR,
     &  HBOR,UBOR,VBOR,CHBORD,TBOR,ATBOR,BTBOR)
!
        USE BIEF
        USE DECLARATIONS_SPECIAL
        IMPLICIT NONE
        INTEGER, INTENT(IN)    :: NPTFR,NPTFRX,NTRAC
        INTEGER,INTENT(INOUT) :: LIHBOR(NPTFRX),LIUBOR(NPTFRX)
        INTEGER,INTENT(INOUT) :: LIVBOR(NPTFRX)
        TYPE(BIEF_OBJ)  , INTENT(INOUT) :: LITBOR
        DOUBLE PRECISION, INTENT(INOUT) :: UBOR(NPTFRX,2),VBOR(NPTFRX,2)
        DOUBLE PRECISION, INTENT(INOUT) :: HBOR(NPTFRX)
        DOUBLE PRECISION, INTENT(INOUT) :: CHBORD(NPTFRX)
        TYPE(BIEF_OBJ)  , INTENT(INOUT) :: TBOR, ATBOR, BTBOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE RANBO
     &(NBOR,KP1BOR,IFABOR,IKLE,NCOLOR,TRAV1,NPTFR,X,Y,NCOLFR,
     & NDP,NPOIN,NELEM,NELMAX,MESH)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NDP,NELMAX,MESH,NELEM,NPOIN
      INTEGER, INTENT(INOUT) :: NPTFR
      INTEGER, INTENT(INOUT) :: NBOR(*),KP1BOR(*),NCOLFR(*)
      INTEGER, INTENT(INOUT) :: TRAV1(NPOIN,2)
      INTEGER, INTENT(IN)    :: IFABOR(NELMAX,*),IKLE(NELMAX,4)
      INTEGER, INTENT(IN)    :: NCOLOR(*)
      DOUBLE PRECISION, INTENT(IN) :: X(NPOIN),Y(NPOIN)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CIRCUL
     &(IKLE,ITEST1 ,IELEM,I1,I2,I3,X,Y,NNELMAX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NNELMAX
      INTEGER, INTENT(IN) :: IELEM
      INTEGER, INTENT(INOUT) :: IKLE(NNELMAX,4)
      INTEGER, INTENT(IN) :: I1 , I2 , I3
      INTEGER, INTENT(INOUT) :: ITEST1
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CONVERTER
     &(LOC_INPFILE,LOC_LOGFILE,LOC_BNDFILE,
     & LOC_OUTFILE,LOC_OUTLOGFILE,LOC_OUTBNDFILE)
          USE DECLARATIONS_STBTEL
      IMPLICIT NONE
      CHARACTER(LEN=MAXLENHARD), INTENT(IN) :: LOC_INPFILE
      CHARACTER(LEN=MAXLENHARD), INTENT(IN) :: LOC_LOGFILE
      CHARACTER(LEN=MAXLENHARD), INTENT(IN) :: LOC_BNDFILE
      CHARACTER(LEN=MAXLENHARD), INTENT(IN) :: LOC_OUTFILE
      CHARACTER(LEN=MAXLENHARD), INTENT(IN) :: LOC_OUTLOGFILE
      CHARACTER(LEN=MAXLENHARD), INTENT(IN) :: LOC_OUTBNDFILE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE CORDEP
     &(IKLE,LGVEC)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4)
      INTEGER, INTENT(IN) :: LGVEC
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DECOUP
     &(ISURC,X,Y,IKLE,NCOLOR,IFABOR, NELEM2,NPOIN2,COLOR)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(INOUT) :: X(*) , Y(*)
      INTEGER, INTENT(IN) :: ISURC
      INTEGER, INTENT(INOUT) :: NELEM2 , NPOIN2
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4) , NCOLOR(*)
      INTEGER, INTENT(INOUT) :: IFABOR(NELMAX,*)
      LOGICAL, INTENT(INOUT) :: COLOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DEPARR
     &(IKLE,NDEPAR,LGVEC)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: LGVEC
      INTEGER, INTENT(INOUT) :: NDEPAR
      INTEGER, INTENT(IN) :: IKLE(NELMAX,4)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE DYNAMI
     &(NPTFR,NBOR,LIHBOR,LIUBOR,LIVBOR,LITBOR,NCOLFR,MAILLE,NLIM)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NPTFR, NLIM
      INTEGER, INTENT(IN) :: NBOR(*) , NCOLFR(*)
      INTEGER, INTENT(INOUT) :: LIHBOR(*) , LIUBOR(*)
      INTEGER, INTENT(INOUT) :: LIVBOR(*) ,LITBOR(*)
      CHARACTER(LEN=9), INTENT(IN) :: MAILLE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE ECRRES
     &(VAINIT,IKINIT,NPINIT,NEINIT,SHP,ELT,NPOIN,NPOIN1,NPMAX,W,
     & X,ZF,NSFOND,NCOLOR,COLOR,VAR,NVARIN,NVAROU,STD,NDP,IKLES,
     & STOTOT,TPSFIN,NGEO,NRES)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NPINIT,NEINIT,NPOIN,NPMAX
      INTEGER, INTENT(INOUT) :: NPOIN1
      DOUBLE PRECISION, INTENT(INOUT) :: VAINIT(NPINIT)
      DOUBLE PRECISION, INTENT(IN) :: SHP(NPMAX,3)
      INTEGER, INTENT(IN) :: IKINIT(NEINIT,3),ELT(NPOIN)
      REAL, INTENT(INOUT) :: W(*)
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: ZF(NPOIN)
      INTEGER, INTENT(IN) :: NSFOND
      INTEGER, INTENT(IN) :: NCOLOR(NPOIN)
      LOGICAL, INTENT(IN) :: COLOR
      DOUBLE PRECISION, INTENT(INOUT) :: VAR(NPOIN)
      INTEGER, INTENT(IN) :: NVARIN,NVAROU
      CHARACTER(LEN=3), INTENT(IN) :: STD
      INTEGER, INTENT(IN) :: NDP
      INTEGER, INTENT(INOUT) :: IKLES(NDP,NEINIT)
      LOGICAL, INTENT(IN) :: STOTOT
      DOUBLE PRECISION, INTENT(IN) :: TPSFIN(1)
      INTEGER, INTENT(IN) :: NGEO,NRES
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ECRSEL
     &(VAINIT,IKINIT,NPINIT,NEINIT,SHP,ELT,NPOIN,NPOIN1,NPMAX,W,
     & X,ZF,NSFOND,NCOLOR,COLOR,VAR,NVARIN,NVAROU,NVAR2,STD,FUSION,
     & NRES,NGEO,NFO1,MAILLE)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NPINIT,NEINIT,NPOIN,NPMAX
      INTEGER, INTENT(INOUT) :: NPOIN1,NVAR2
      DOUBLE PRECISION, INTENT(INOUT) :: VAINIT(NPINIT)
      DOUBLE PRECISION, INTENT(IN) :: SHP(NPMAX,3)
      INTEGER, INTENT(IN) :: IKINIT(NEINIT,3),ELT(NPOIN)
      REAL, INTENT(INOUT) :: W(*)
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: ZF(NPOIN)
      INTEGER, INTENT(IN) :: NSFOND
      INTEGER, INTENT(IN) :: NCOLOR(NPOIN)
      LOGICAL, INTENT(IN) :: COLOR,FUSION
      DOUBLE PRECISION, INTENT(INOUT) :: VAR(NPOIN)
      INTEGER, INTENT(IN) :: NVARIN,NVAROU
      CHARACTER(LEN=3), INTENT(IN) :: STD
      INTEGER, INTENT(IN) :: NGEO,NRES,NFO1
      CHARACTER(LEN=9), INTENT(IN) :: MAILLE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE ELMPB
     &(NBPB,NUMPB,X,Y,IKLE,NCOLOR,ISDRY,NEW)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX,NPMAX
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4), ISDRY(NPMAX), NEW(NPMAX)
      INTEGER, INTENT(INOUT) :: NCOLOR(NPMAX)
      INTEGER,INTENT(IN) :: NBPB, NUMPB(100)
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPMAX) , Y(NPMAX)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE ELMSEC
     &( ELPSEC, SEUSEC, TPSFIN, X, Y, IKLE, NCOLOR, ISDRY,
     &  IHAUT, NVAR, H, WORK, NEW, STD, NGEO )
      USE DECLARATIONS_STBTEL, ONLY: NELMAX,NPMAX
      IMPLICIT NONE
!
      LOGICAL, INTENT(IN) :: ELPSEC
      DOUBLE PRECISION, INTENT(IN) :: SEUSEC
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPMAX),Y(NPMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: H(NPMAX),TPSFIN(1)
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4), ISDRY(NPMAX), NEW(NPMAX)
      INTEGER, INTENT(INOUT) :: NCOLOR(NPMAX)
      INTEGER, INTENT(IN) :: IHAUT, NVAR
      REAL, INTENT(INOUT) :: WORK(*)
      INTEGER, INTENT(IN) :: NGEO
      CHARACTER(LEN=3), INTENT(IN) :: STD
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE EXTRAC
     &(X,Y,SOM,IKLE,INDIC,NELEM,NELMAX,NPOIN,NSOM,PROJEC)
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: NELEM
      INTEGER, INTENT(IN) :: NELMAX,NPOIN,NSOM
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,3),INDIC(NPOIN)
      DOUBLE PRECISION, INTENT(INOUT) :: X(NPOIN),Y(NPOIN)
      DOUBLE PRECISION, INTENT(IN) :: SOM(10,2)
      LOGICAL, INTENT(IN) :: PROJEC
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE FM3SEL
     &(X,Y,NPOIN,NBOR,NFIC,STD,NVAR,TEXTE,TEXTLU,VARCLA,NVARCL,
     & TITRE,SORLEO,NSOR,W,IKLE,
     & IKLES,ITRAV,NELEM,NPTFR,NDP,MXPTVS,MXELVS,DATE,TIME,
     & DEBU,SUIT,ECRI,LISTIN,IPARAM,IPOBO)
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(IN) :: X(*),Y(*)
      REAL, INTENT(INOUT) :: W(*)
!                     IKLE(NELEM,NDP) IKLES(NDP,NELEM)
      INTEGER, INTENT(IN) :: NBOR(*)
      INTEGER, INTENT(INOUT) :: IKLE(*),IKLES(*),ITRAV(*)
      INTEGER, INTENT(INOUT) :: NPOIN,NVAR,MXPTVS,MXELVS,TIME(3),DATE(3)
      INTEGER, INTENT(IN) :: NFIC,NVARCL,NSOR
      INTEGER, INTENT(INOUT) :: NELEM,NPTFR,NDP
      INTEGER, INTENT(IN) :: IPARAM(10),IPOBO(*)
      LOGICAL, INTENT(IN) :: DEBU,SUIT,ECRI,LISTIN,SORLEO(*)
      CHARACTER(LEN=3), INTENT(IN) :: STD
      CHARACTER(LEN=72), INTENT(IN) :: TITRE
!                        NSOR      NSOR+NVARCL
      CHARACTER(LEN=32), INTENT(IN) :: TEXTE(*),VARCLA(NVARCL)
      CHARACTER(LEN=32), INTENT(INOUT) :: TEXTLU(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE IMPRIM
     &(NPOIN1,NPOIN,TYPELE,NELEM,TITRE,MAILLE,PRECIS)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NPOIN1, NPOIN, NELEM
!
      CHARACTER(LEN=11), INTENT(IN) :: TYPELE
      CHARACTER(LEN=80), INTENT(IN) :: TITRE
      CHARACTER(LEN=9), INTENT(IN) ::  MAILLE
      CHARACTER(LEN=6), INTENT(IN) ::  PRECIS
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE INIADC
     &(NPOIN1,TYPELE,NSFOND,IHAUT,NGEO,TITRE)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NGEO
      INTEGER, INTENT(INOUT) :: NPOIN1 , NSFOND
      INTEGER, INTENT(INOUT) :: IHAUT
!
      CHARACTER(LEN=80), INTENT(INOUT) :: TITRE
      CHARACTER(LEN=11), INTENT(INOUT) :: TYPELE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE INIFAS
     &(TYPELE,NGEO)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) ::       NGEO
      CHARACTER(LEN=*), INTENT(INOUT) :: TYPELE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
        SUBROUTINE INISEL
     &(NPOIN1,TYPELE,STD,NSFOND,FUSION,IHAUT,NGEO,NFO1)
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: NPOIN1,NSFOND
      CHARACTER(LEN=11), INTENT(INOUT) :: TYPELE
      CHARACTER(LEN=3), INTENT(IN) :: STD
      LOGICAL, INTENT(IN) :: FUSION
      INTEGER, INTENT(INOUT) :: IHAUT
      INTEGER, INTENT(IN) :: NGEO , NFO1
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE INISIM
     &(NPOIN1,TYPELE,INOP5,NGEO)
      IMPLICIT NONE
      !
      INTEGER, INTENT(INOUT) :: NPOIN1,INOP5
      INTEGER, INTENT(IN) :: NGEO
      CHARACTER(LEN=11), INTENT(INOUT) :: TYPELE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE INISTB
     &(NPOIN1,TYPELE,MAILLE,PRECIS,NGEO,NSEC2,NSEC11,NSEC12)
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: NPOIN1
      CHARACTER(LEN=11), INTENT(INOUT) :: TYPELE
      CHARACTER(LEN=9), INTENT(IN) ::  MAILLE
      CHARACTER(LEN=6), INTENT(INOUT) ::  PRECIS
      INTEGER, INTENT(IN) :: NSEC11 , NSEC12 , NGEO, NSEC2
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE INITRI
     &( NPOIN1,TYPELE,NGEO,NFO1)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NGEO, NFO1
      INTEGER, INTENT(INOUT) :: NPOIN1
      CHARACTER(LEN=*), INTENT(INOUT) :: TYPELE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE INTERP
     &(XINIT , YINIT , IKINIT , NPINIT , NEINIT ,
     & X , Y , NPOIN , NPMAX , SHP , ELT)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NPINIT, NEINIT, NPOIN,NPMAX
      DOUBLE PRECISION, INTENT(IN) :: XINIT(NPINIT) , YINIT(NPINIT)
      INTEGER, INTENT(IN) :: IKINIT(NEINIT,3)
      INTEGER, INTENT(INOUT) :: ELT(NPMAX)
      DOUBLE PRECISION, INTENT(IN) :: X(NPMAX) , Y(NPMAX)
      DOUBLE PRECISION, INTENT(INOUT) :: SHP(NPMAX,3)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE LECADC
     &( X , Y , ZF , IKLE , NGEO )
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NGEO
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4)
      DOUBLE PRECISION, INTENT(INOUT) :: X(*) , Y(*),ZF(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE LECFAS
     & (X, Y, IKLE, NCOLOR, TFAST1, TFAST2, ADDFAS,
     &  NGEO , NFO1)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NGEO, NFO1
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4)
      INTEGER, INTENT(INOUT) :: NCOLOR(*)
      INTEGER, INTENT(INOUT) :: TFAST1(*),TFAST2(*)
      LOGICAL, INTENT(IN) :: ADDFAS
      DOUBLE PRECISION, INTENT(INOUT) :: X(*), Y(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE LECFON
     &( XRELV , YRELV , ZRELV , NBAT , NFOND , NBFOND ,  NP ,
     &  NPT , FONTRI , CORTRI , MAILLE, NGEO )
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(INOUT) :: XRELV(*) , YRELV(*) , ZRELV(*)
      INTEGER, INTENT(IN) :: NFOND(*) , NBAT , NBFOND
      INTEGER, INTENT(INOUT) :: NP(5), NPT
      LOGICAL, INTENT(IN) :: FONTRI
      DOUBLE PRECISION, INTENT(IN) :: CORTRI
      CHARACTER(LEN=9), INTENT(IN) ::  MAILLE
      INTEGER, INTENT(IN) :: NGEO
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE LECSEL
     &(XINIT,YINIT,IKINIT,NPINIT,NEINIT,X,Y,IKLE,IKLES,W,TITRE,TEXTE,
     & NVARIN,NVAR2,STD,NCOLOR,FUSION,NGEO,NFO1,IPOBO,IPARAM,DATE,
     & TIME)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX,NELEM,NDP
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(INOUT) :: XINIT(*), YINIT(*), X(*), Y(*)
      REAL, INTENT(INOUT) :: W(*)
      INTEGER, INTENT(IN) :: NGEO , NFO1
      INTEGER, INTENT(INOUT) :: IPARAM(10),DATE(3),TIME(3)
      INTEGER, INTENT(INOUT) :: NEINIT , NPINIT
      INTEGER, INTENT(INOUT) :: NVARIN , NVAR2
      INTEGER, INTENT(INOUT) :: IKINIT(NELEM,NDP)
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,NDP),IKLES(NDP,NELEM)
      INTEGER, INTENT(INOUT) :: NCOLOR(*),IPOBO(*)
      LOGICAL, INTENT(IN) :: FUSION
      CHARACTER(LEN=72), INTENT(INOUT) :: TITRE
      CHARACTER(LEN=32), INTENT(INOUT) :: TEXTE(26)
      CHARACTER(LEN=3), INTENT(IN) ::  STD
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE LECSIM
     &( X , Y , IKLE , NCOLOR , TITRE , NOP5 , NGEO )
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(INOUT) :: X(*) , Y(*)
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4) , NCOLOR(*)
      CHARACTER(LEN=80), INTENT(INOUT) :: TITRE
      INTEGER, INTENT(INOUT) :: NOP5(*)
      INTEGER, INTENT(IN) :: NGEO
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE LECSTB
     &( X , Y ,IKLE , NCOLOR , TITRE , NPOIN1 ,
     &  NGEO , NSEC2,NSEC3,NSEC11,NSEC12)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(INOUT) :: X(*) , Y(*)
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4) , NCOLOR(*)
      CHARACTER(LEN=80), INTENT(INOUT) :: TITRE
      INTEGER, INTENT(IN) :: NPOIN1
      INTEGER, INTENT(IN) :: NGEO
      INTEGER, INTENT(IN) :: NSEC11 , NSEC12 , NSEC2 , NSEC3
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE LECTRI
     & (X, Y, IKLE, NCOLOR,NGEO , NFO1)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NGEO, NFO1
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4)
      INTEGER, INTENT(INOUT) :: NCOLOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: X(*), Y(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE PRESEL
     &(IKLE,TRAV1,NELEM,NELMAX,NDP,TEXTE,NBFOND,SORLEO,COLOR,
     & NSFOND,NVARIN,NVAROU,MAILLE)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NDP,NELEM,NELMAX,NBFOND,NVARIN
      INTEGER, INTENT(INOUT) :: NSFOND,NVAROU
      INTEGER, INTENT(INOUT) :: TRAV1(NELEM,NDP)
      INTEGER, INTENT(IN) :: IKLE(NELMAX,NDP)
      CHARACTER(LEN=32), INTENT(INOUT) :: TEXTE(26)
      CHARACTER(LEN=9), INTENT(IN) :: MAILLE
      LOGICAL, INTENT(INOUT) :: SORLEO(26),COLOR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE PROJEC
     &(X , Y , ZF , XRELV , YRELV , ZRELV , NBAT ,
     & NBOR , NPTFR , NFOND , NBFOND , FOND , DM ,
     & FONTRI , CORTRI , MAILLE,NGEO,KP1BOR)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: NPTFR, NBAT, NBFOND
      INTEGER, INTENT(IN) :: NFOND(*) , NBOR(NPTFR,2)
      INTEGER, INTENT(IN) :: NGEO, KP1BOR(NPTFR)
      DOUBLE PRECISION, INTENT(INOUT) :: XRELV(*) , YRELV(*) , ZRELV(*)
      DOUBLE PRECISION, INTENT(IN) :: X(*) , Y(*) , DM
      DOUBLE PRECISION, INTENT(IN) :: CORTRI
      DOUBLE PRECISION, INTENT(INOUT) :: ZF(*)
      CHARACTER(LEN=72), INTENT(IN) :: FOND(NBFOND)
      CHARACTER(LEN=9), INTENT(IN) ::  MAILLE
      LOGICAL, INTENT(IN) :: FONTRI
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE REMAIL
     &(IKLE,NCOLOR,NEW,X,Y,EPSI,NDP,NPOIN,NELEM,NELMAX)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN)    ::  NDP , NELMAX
      INTEGER, INTENT(INOUT) ::  NPOIN, NELEM
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4) , NEW(*) , NCOLOR(*)
      DOUBLE PRECISION, INTENT(INOUT) :: X(*) , Y(*), EPSI
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE RENUM
     &(X,Y,W,IKLE,NBOR,TRAV1,TRAV2,TAB,NCOLOR,COLOR,NPTFR)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(INOUT) :: X(*) , Y(*) , W(*)
      INTEGER, INTENT(INOUT) :: TRAV1(*) , TRAV2(*)
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,3) , NCOLOR(*) , NBOR(*)
      INTEGER, INTENT(INOUT) :: TAB(*)
      LOGICAL, INTENT(IN) :: COLOR
      INTEGER, INTENT(IN) :: NPTFR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE SHUFLE
     &(IKLE,X)
      USE DECLARATIONS_STBTEL
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4)
      DOUBLE PRECISION, INTENT(IN) :: X(*)
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE STBTEL
     &( NPOIN1 , TYPELE , NFOND , PRECIS , NSFOND , TITRE)
      IMPLICIT NONE
!
      INTEGER, INTENT(INOUT) :: NPOIN1
      CHARACTER(LEN=11), INTENT(INOUT) :: TYPELE
      INTEGER, INTENT(IN) :: NFOND(5)
      CHARACTER(LEN=6), INTENT(INOUT) :: PRECIS
      INTEGER, INTENT(INOUT) :: NSFOND
      CHARACTER(LEN=80), INTENT(INOUT) :: TITRE
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE SURCON
     &(X,Y,IKLE,IPO,NBOR,NPTFR,NCOLOR,IFABOR,COLOR)
      USE DECLARATIONS_STBTEL, ONLY: NELMAX
      IMPLICIT NONE
!
      DOUBLE PRECISION, INTENT(INOUT) :: X(*) , Y(*)
      INTEGER, INTENT(INOUT) :: NBOR(*) , IKLE(NELMAX,4) , NCOLOR(*)
      INTEGER, INTENT(INOUT) :: IFABOR(NELMAX,*) , IPO(*)
      LOGICAL, INTENT(INOUT) :: COLOR
      INTEGER, INTENT(IN) :: NPTFR
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE VERIFI
     &(X,Y,IKLE,NCOLOR,TRAV1,EPSI,MESH,NDP,NPOIN,NELEM,NELMAX)
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) ::  MESH , NDP , NELMAX
      INTEGER, INTENT(INOUT) :: NPOIN, NELEM
      INTEGER, INTENT(INOUT) :: IKLE(NELMAX,4) , NCOLOR(*)
      INTEGER, INTENT(INOUT) :: TRAV1(*)
!
      DOUBLE PRECISION, INTENT(INOUT) :: X(*) , Y(*)
      DOUBLE PRECISION, INTENT(INOUT) :: EPSI
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                        SUBROUTINE VERIFS
     &(IFABOR,IKLE,TRAV1,NPTFR,NUMPB,NBPB)
      USE DECLARATIONS_STBTEL
      IMPLICIT NONE
!
      INTEGER, INTENT(IN) :: IFABOR(NELMAX,*) , IKLE(NELMAX,4)
      INTEGER, INTENT(INOUT) :: TRAV1(NPOIN,2)
      INTEGER, INTENT(INOUT) :: NPTFR
      INTEGER, INTENT(INOUT) :: NUMPB(100), NBPB
        END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                     SUBROUTINE WRITESELLIM
     &(NLIM,LIHBOR,LIUBOR,LIVBOR,HBOR,UBOR,VBOR,
     & CHBORD,NBOR,NPMAX,NPTFR)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NLIM
      INTEGER, INTENT(IN)    :: NPTFR
      INTEGER, INTENT(IN)    :: NPMAX
      INTEGER, INTENT(INOUT) :: LIUBOR(NPMAX),LIVBOR(NPMAX)
      INTEGER, INTENT(INOUT) :: LIHBOR(NPMAX)
      INTEGER, INTENT(INOUT) :: NBOR(NPMAX)
      DOUBLE PRECISION,  INTENT(INOUT) :: UBOR(NPMAX),VBOR(NPMAX)
      DOUBLE PRECISION,  INTENT(INOUT) :: HBOR(NPMAX),CHBORD(NPMAX)
      END SUBROUTINE
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      INTERFACE
                     SUBROUTINE LECSELLIM
     &(NLIM,LIHBOR,LIUBOR,LIVBOR,HBOR,UBOR,VBOR,
     & CHBORD,NBOR,NPTFR,NPTFR2)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: NLIM
      INTEGER, INTENT(IN)    :: NPTFR
      INTEGER, INTENT(INOUT) :: LIUBOR(NPTFR),LIVBOR(NPTFR)
      INTEGER, INTENT(INOUT) :: LIHBOR(NPTFR)
      INTEGER, INTENT(INOUT) :: NBOR(NPTFR)
      INTEGER, INTENT(OUT) :: NPTFR2
      DOUBLE PRECISION,  INTENT(INOUT) :: UBOR(NPTFR,2),VBOR(NPTFR,2)
      DOUBLE PRECISION,  INTENT(INOUT) :: HBOR(NPTFR),CHBORD(NPTFR)
      END SUBROUTINE LECSELLIM
      END INTERFACE
!
!-----------------------------------------------------------------------
!
      END MODULE INTERFACE_STBTEL
