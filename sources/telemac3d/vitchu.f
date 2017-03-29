!                    *****************
                     SUBROUTINE VITCHU
!                    *****************
!
     &(WCHU,WCHU0,TURBWC,U,V,W,H,RUGOF,LISRUF,TURBA,TURBB,
     & TRAV1,TRAV2,TRAV3,TRAV4,S,MESH3D,IELM3,NPOIN2,NPOIN3,NPLAN,
     & NCLASS,MSK,MASKEL,UETCAR,SED_TRA,HN,HMIN,
     & FLOC,FLOC_TYPE,HINDER,HIND_TYPE,CGEL,CINI,TYPE_OF_SEDIMENT,
     & SED_CO,SED_NCO)
!
!***********************************************************************
! TELEMAC3D   V7P2                                  31/08/2016
!***********************************************************************
!
!brief    COMPUTES THE SETTLING VELOCITY AS A FUNCTION
!+                OF TEMPERATURE, SALINITY AND CONCENTRATION OF
!+                SUSPENDED SEDIMENT.
!
!history  C LE NORMANT (LNH)
!+        01/08/91
!+        V5P1
!+
!
!history  JACEK A. JANKOWSKI PINXIT
!+        **/03/99
!+
!+   FORTRAN95 VERSION
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
!history  C. VILLARET & T. BENSON & D. KELLY (HR-WALLINGFORD)
!+        27/02/2014
!+        V7P0
!+   New developments in sediment merged on 25/02/2014.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| IELM3          |-->| DISCRETISATION TYPE FOR 3D
!| H              |-->| WATER DEPTH
!| HN             |-->| WATER DEPTH AT TIME N
!| HMIN           |-->| THRESHOLD FOR EROSION FLUXES ON TIDAL FLATS
!| LISRUF         |-->| TURBULENCE MODEL FOR BOTTOM
!| MASKEL         |-->| MASKING OF ELEMENTS
!|                |   | =1. : NORMAL   =0. : MASKED ELEMENT
!| MESH3D         |<->| 3D MESH
!| MSK            |-->| IF YES, THERE IS MASKED ELEMENTS.
!| NPOIN2         |-->| NUMBER OF POINTS IN 2D
!| NPOIN3         |-->| NUMBER OF POINTS IN 3D
!| NPLAN          |-->| NUMBER OF PLANES IN THE 3D MESH OF PRISMS
!| NCLASS         |-->| NUMBER OF SEDIMENT CLASSES
!| RUGOF          |-->| FRICTION COEFFICIENT
!| S              |-->| VOID STRUCTURE
!| SOULSBYWC      |-->| SWITCH FOR SOULSBY FLOCCULATION FORMULA
!| TA             |-->| TRACER CONCENTRATION (LAST ONE, NTRAC, IS SED)
!| TOB            |-->| BED SHEAR STRESS (INCLUDES TURBULENCE DAMPING)
!| TRAV1          |<->| WORK ARRAY
!| TRAV2          |<->| WORK ARRAY
!| TRAV3          |<->| WORK ARRAY
!| TURBA          |-->| FLOCCULATION COEFFICIENT
!| TURBB          |-->| COEFFICIENT RELATIVE TO FLOC DESTRUCTION
!| TURBWC         |-->| SWITCH FOR TURBULENT FLOC BREAKUP
!| U,V,W          |-->| VELOCITY COMPONENTS
!| UETCAR         |-->| SQUARE OF THE FRICTION VELOCITY
!| WCHU           |<--| SEDIMENT SETTLING VELOCITY (M/S)
!| WCHU0          |-->| CONSTANT SEDIMENT SETTLING VELOCITY (M/S)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      TYPE(BIEF_MESH),INTENT(INOUT) :: MESH3D
      TYPE(BIEF_OBJ), INTENT(INOUT) :: WCHU,TRAV1,TRAV2,TRAV3,TRAV4
      TYPE(BIEF_OBJ), INTENT(IN)    :: MASKEL,S,H,HN,RUGOF,U,V,W
      TYPE(BIEF_OBJ), INTENT(IN)    :: SED_TRA,UETCAR
      LOGICAL, INTENT(IN)           :: TURBWC,MSK,HINDER, FLOC
      INTEGER, INTENT(IN)           :: NPOIN2,NPOIN3,NPLAN,IELM3,LISRUF
      INTEGER, INTENT(IN)           :: NCLASS,HIND_TYPE, FLOC_TYPE
      INTEGER, INTENT(IN)           :: TYPE_OF_SEDIMENT(NCLASS)
      INTEGER, INTENT(IN)           :: SED_CO,SED_NCO 
      DOUBLE PRECISION, INTENT(IN)  :: WCHU0(NCLASS)
      DOUBLE PRECISION, INTENT(IN)  :: TURBA,TURBB,HMIN,CGEL
      DOUBLE PRECISION, INTENT(IN)  :: CINI
      INTEGER                       :: ICLASS
      
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     CONSTANT VALUE GIVEN HERE WORKING FOR COHESIVE AND NON-COHESIVE CLASSES
!
      DO ICLASS=1,NCLASS
        CALL OS( 'X=C     ' , X=WCHU%ADR(ICLASS)%P, C=WCHU0(ICLASS) )
      END DO

! computation the total suspended concentration in TRAV2
! for cohesive and non-cohesive  -- HR !!!
      IF(FLOC_TYPE==2.OR.HINDER) THEN
        CALL OS( 'X=0     ' , X=TRAV4)
        DO ICLASS=1,NCLASS
          CALL OS( 'X=X+Y   ' , X=TRAV4, Y=SED_TRA%ADR(ICLASS)%P )
        END DO
      END IF

!
! 1. FLOCCULATION
!
      IF(FLOC) THEN
        IF(FLOC_TYPE.EQ.1) THEN
!
!         APPLY REDUCTION DUE TO TURBULENT BREAKUP OF FLOCS
!         ONLY FOR COHESIVE SEDIMENTS
!
          DO ICLASS=1,NCLASS
            IF(TYPE_OF_SEDIMENT(ICLASS)==SED_CO) THEN
              CALL WCTURB(WCHU%ADR(ICLASS)%P%R,WCHU0(ICLASS),
     &                U,V,W,H,RUGOF,LISRUF,
     &                TRAV1,TRAV2,TRAV3, S,MESH3D,IELM3,NPOIN2,
     &                NPLAN,TURBA,TURBB,MSK,MASKEL,UETCAR)
            ENDIF
          END DO
!
        ELSEIF(FLOC_TYPE.EQ.2) THEN
!
!         SOULSBY FLOC MODEL
!

          IF(HINDER) THEN
! TRAV4: total concentration including non-cohesive classes
! TRAV1: is limited to the hindering concentration
            CALL OS('X=-(Y,C)',X=TRAV1,Y=TRAV4,C=CINI)
          ELSE
            CALL OS('X=Y     ',X=TRAV1,Y=TRAV4)
          ENDIF
!
          DO ICLASS=1,NCLASS
            IF(TYPE_OF_SEDIMENT(ICLASS)==SED_CO) THEN
              CALL SOULSBYFLOC3D(WCHU%ADR(ICLASS)%P,TRAV1%R,
     &            MESH3D,NPOIN2,NPOIN3,NPLAN,HN,HMIN,UETCAR%R)
!
            ENDIF
          END DO
        ELSE
!
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'FORMULE POUR FLOCULATION INCONNUE : ',FLOC_TYPE
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'UNKNOWN FLOCCULATION FORMULA: ',FLOC_TYPE
          ENDIF
          CALL PLANTE(1)
          STOP
!
        ENDIF
      ENDIF
!
! 2. HINDERED SETTLING
!
! LIMIT THE CONCENTRATION TO CINI (IF HINDERED SETTLING IS ON) (tbe comment: no!
!                                  It only gets limited for the floc model)
!
      IF(HINDER) THEN
!###>TBE - we don't limit the concentration here, otherwise hindering won't happen
        !CALL OS('X=-(Y,C)',X=TRAV1,Y=TA%ADR(NTRAC)%P,C=CINI)
!###<TBE
        ! this copy of concentration is a bit unecessary.
        ! Would be better to pass
        ! a pointer and use double array in WCHIND
        DO ICLASS=1,NCLASS
          IF(TYPE_OF_SEDIMENT(ICLASS)==SED_CO) THEN
! TRAV4: total concentration including non-cohesive classes
           CALL WCHIND(WCHU%ADR(ICLASS)%P,TRAV4,
     &          CINI,CGEL,NPOIN3,HIND_TYPE)
          ENDIF
        END DO
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
