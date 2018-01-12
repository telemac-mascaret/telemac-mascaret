!                    ***************************
                     SUBROUTINE BEDLOAD_SOLVS_VF
!                    ***************************
!
     &(MESH,QSX,QSY,LIMTEC,UNSV2D,EBOR,BREACH,NSEG,NPTFR,NPOIN,
     & KENT,KDIR,KDDL,DT,ZFCL,FLUX,CSF_SABLE,FLBCLA,
     & LIQBOR,QBOR,NUBO,VNOIN,EVCL_M,RATIO_SAND,XMVS)
!
!***********************************************************************
! SISYPHE   V7P0                                      03/06/2014
!***********************************************************************
!
!brief    SOLVES A PART OF EXNER EQUATION (DIVERGENCE TERM)
!         WITH THE FINITE VOLUME METHOD.
!
!history  M. GONZALES DE LINARES
!+        07/05/2002
!+        V5P5
!+   First version.
!
!history  J-M HERVOUET (EDF, LNHE)
!+        30/10/2007
!+        V5P8
!+   UNSV2D +DIRICL DELETED
!
!history  JMH
!+        15/09/2009
!+
!+   KDIR KDDL ADDED (WERE HARD-CODED BEFORE !!!)
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
!history  J-M HERVOUET (EDF-LNHE)
!+        21/02/2012
!+        V6P2
!+  Corrections for a perfect mass balance: FLBTRA built to be used in
!+  bilan_sisyphe, coefficient AVA added in Dirichlet value, QBOR
!+  dealt with.
!
!history  R.ATA (EDF-LNHE)
!+        02/06/2014
!+        V7P0
!+  Corrections of normals and nubo tables
!+  after changes in FV data structure of Telemac2d
!+
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| BREACH         |<->| INDICATOR FOR NON ERODIBLE BED
!| DT             |-->| TIME STEP
!| EBOR           |<->| BOUNDARY CONDITION FOR BED EVOLUTION (DIRICHLET)
!| FLBCLA         |<->| FLUXES AT BOUNDARY FOR THE CLASS
!| FLUX           |<->| SEDIMENT FLUX
!| KDIR           |-->| CONVENTION FOR DIRICHLET VALUE
!| KDDL           |-->| CONVENTION FOR DEGREE OF FREEDOM
!| KENT           |-->| CONVENTION FOR PRESCRIBED VALUE AT ENTRANCE
!| LIMTEC         |-->| TECHNICAL BOUNDARY CONDITIONS FOR BED EVOLUTION
!| LIQBOR         |-->| TYPE OF BOUNDARY CONDITIONS FOR BEDLOAD DISCHARGE
!| EVCL_M         |<->| MASS EVOLUTION DURING A TIME STEP
!| MESH           |<->| MESH STRUCTURE
!| NPOIN          |-->| NUMBER OF POINTS
!| NPTFR          |-->| NUMBER OF BOUNDARY POINTS
!| NSEG           |-->| NUMBER OF SEGMENTS
!| NUBO           |-->| GLOBAL NUMBER OF EDGE EXTREMITIES
!| QBOR           |-->| PRESCRIBED BEDLOAD DISCHARGES
!| QSX            |<->| BEDLOAD TRANSPORT RATE X-DIRECTION
!| QSY            |<->| BEDLOAD TRANSPORT RATE Y-DIRECTION
!| RATIO_SAND     |-->| MASS FRACTION OF SAND, FOR THE CLASS
!| UNSV2D         |-->| INVERSE OF INTEGRALS OF TEST FUNCTIONS
!| VNOIN          |-->| OUTWARD UNIT NORMALS
!| XMVS           |-->| SEDIMENT DENSITY
!| ZFCL           |<->| BEDLOAD EVOLUTION FOR EACH SEDIMENT CLASS
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE INTERFACE_SISYPHE, EX_BEDLOAD_SOLVS_VF => BEDLOAD_SOLVS_VF
      USE BIEF
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      TYPE(BIEF_MESH),  INTENT(INOUT) :: MESH
      TYPE(BIEF_OBJ),   INTENT(IN)    :: QSX,QSY,LIMTEC,UNSV2D,EBOR
      TYPE(BIEF_OBJ),   INTENT(IN)    :: BREACH,LIQBOR,QBOR
      INTEGER,          INTENT(IN)    :: NSEG,NPTFR,NPOIN,KENT,KDIR,KDDL
      DOUBLE PRECISION, INTENT(IN)    :: DT,CSF_SABLE
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: FLBCLA,ZFCL,FLUX
      TYPE(BIEF_OBJ),   INTENT(INOUT) :: EVCL_M
      INTEGER, INTENT(IN)             :: NUBO(2,NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: VNOIN(3,NSEG)
      DOUBLE PRECISION, INTENT(IN)    :: RATIO_SAND(NPOIN),XMVS
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER          :: ISEGIN,K,N,IEL,IEL1,IEL2
      DOUBLE PRECISION :: QSMOY1,QSMOY2,QSP,VNOIN1,VNOIN2,RNORM,XN,YN
      DOUBLE PRECISION :: EVCL_MDIR,PROD_SCAL
!
!======================================================================!
!======================================================================!
!                               PROGRAM                                !
!======================================================================!
!======================================================================!
!
!     INTIALISES THE DIVERGENCE
!
      CALL OS('X=0     ', X=FLUX)
!
!     DETERMINES THE OUTGOING FLUX FOR EACH CELL
!
      DO ISEGIN = 1, NSEG
!
        IEL1 = NUBO(1,ISEGIN)
        IEL2 = NUBO(2,ISEGIN)
!
        ! II.1 - SEGMENT LENGTH (RNORM)
        ! ----------------------------------
        VNOIN1 = VNOIN(1,ISEGIN)
        VNOIN2 = VNOIN(2,ISEGIN)
        RNORM  = VNOIN(3,ISEGIN)
        PROD_SCAL= (MESH%X%R(IEL2)-MESH%X%R(IEL1))*VNOIN1+
     &             (MESH%Y%R(IEL2)-MESH%Y%R(IEL1))*VNOIN2
        IF(PROD_SCAL.LT.0.D0)THEN
          IEL1 = NUBO(2,ISEGIN)
          IEL2 = NUBO(1,ISEGIN)
        ENDIF
        ! II.2 - QS FOR THE SEGMENT, BROKEN UP ACCORDING TO X AND Y
        ! ---------------------------------------------
        QSMOY1 = 0.5D0*(QSX%R(IEL1) + QSX%R(IEL2))
        QSMOY2 = 0.5D0*(QSY%R(IEL1) + QSY%R(IEL2))
        ! II.3 - PROJECTS QS FOR THE SEGMENT ONTO THE SEGMENT NORMAL
        ! ------------------------------------------------------------
        QSP = VNOIN1*QSMOY1 + VNOIN2*QSMOY2
        ! II.4 - UPWIND SCHEME ON NODES WITH A "PROBLEM"
        ! ----------------------------------------------
        IF(BREACH%I(IEL1).EQ.1.AND.QSP.GT.0.D0) THEN
          QSMOY1 = QSX%R(IEL1)
          QSMOY2 = QSY%R(IEL1)
        ENDIF
        IF(BREACH%I(IEL2).EQ.1.AND.QSP.LT.0.D0) THEN
          QSMOY1 = QSX%R(IEL2)
          QSMOY2 = QSY%R(IEL2)
        ENDIF
        QSP = VNOIN1*QSMOY1 + VNOIN2*QSMOY2
        ! II.5 - INTEGRATES BY THE SEGMENT LENGTH
        ! ---------------------------------------------
        FLUX%R(IEL1) = FLUX%R(IEL1) + RNORM*QSP
        FLUX%R(IEL2) = FLUX%R(IEL2) - RNORM*QSP
      ENDDO
!
!     BOUNDARIES
!
      DO K = 1 , NPTFR
        IEL = MESH%NBOR%I(K)
!       VECTOR NORMAL TO A BOUNDARY NODE
!       VERSION WHICH IS NOT NORMED
        XN = MESH%XNEBOR%R(K+NPTFR)
        YN = MESH%YNEBOR%R(K+NPTFR)
!
!       ADDING BOUNDARY FLUX AT OPEN BOUNDARIES
!       QBOR HAS PRIORITY HERE, SO IN CASE OF LIQBOR=KENT
!       LIMTEC=KDIR WILL NOT BE CONSIDERED, SEE BEDLOAD_SOLVS_FE
!       FOR MORE EXPLANATIONS
!
!       QBOR AND FLUX MUST HAVE THE SAME DIMENSION !
!
        IF(LIQBOR%I(K).EQ.KENT) THEN
          FLBCLA%R(K) = QBOR%R(K)
          FLUX%R(IEL) = FLUX%R(IEL) + FLBCLA%R(K)
        ELSEIF(LIMTEC%I(K).EQ.KDDL.OR.LIMTEC%I(K).EQ.KDIR) THEN
!         IF KDIR WILL BE UPDATED LATER
          FLBCLA%R(K)= QSX%R(IEL)*XN + QSY%R(IEL)*YN
!         ADDS THE CONTRIBUTION OF THE FLUX ON THE BOUNDARY SEGMENT
          FLUX%R(IEL) = FLUX%R(IEL) + FLBCLA%R(K)
        ELSE
!         NO SEDIMENT FLUX ACCROSS SOLID BOUNDARIES
          FLBCLA%R(K)=0.D0
        ENDIF
      ENDDO
!
!     ASSEMBLING IN PARALLEL
!
      IF(NCSIZE.GT.1) CALL PARCOM(FLUX, 2, MESH)
!
!     SOLVING, NEGATIVE SIGN BECAUSE OUTGOING FLUX IS POSITIVE
!     NOTE JMH: FLUX MUST BE HERE DIV(QS)
!
      CALL OS('X=CYZ   ', X=EVCL_M, Y=FLUX, Z=UNSV2D, C=-DT)
!
!     THIS PART HAS STILL TO BE CHANGED IN MASS!
!     ebor*xmvs replaced by a unique value given by the user??
!
      DO K=1,NPTFR
!                                  PRIORITY OF LIQBOR, SEE ABOVE
        IF(LIMTEC%I(K).EQ.KDIR.AND.LIQBOR%I(K).NE.KENT) THEN
          N = MESH%NBOR%I(K)
!         ZFCLDIR: DIRICHLET VALUE OF EVOLUTION, ZFCL WILL BE DIVIDED BY
!         CSF_SABLE AFTER, AND THEN IT WILL BE AVA(N)*EBOR...
!         when mass: we do not need to multiply by csf_sable
          EVCL_MDIR = RATIO_SAND(N)*EBOR%R(K)*CSF_SABLE*XMVS
!         CORRECTION OF BOUNDARY FLUX TO GET ZFCLDIR
!         flbcla must be previously multiplied by xmvs!!
          FLBCLA%R(K)=FLBCLA%R(K)-(EVCL_MDIR-EVCL_M%R(N))/
     &                (DT*UNSV2D%R(N))
!         EVCL_MDIR FINALLY IMPOSED
          EVCL_M%R(N) = EVCL_MDIR
        ENDIF
!
      ENDDO
!
!======================================================================!
!======================================================================!
!
      RETURN
      END
