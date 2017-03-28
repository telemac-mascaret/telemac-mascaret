!                    *****************
                     SUBROUTINE SEDI3D_Y
!                    *****************
!
     & (MASBED0,EPAI,CONC,HDEP,T2_01,NPOIN2,NPFMAX,NCOUCH,NPF,
     & TASSE,GIBSON,RHOS,VOLU2D,CFDEP, EPAICO, EPAINCO, MIXTE,
     & MASDEP,ESOMT,DEBUG)
!
!
      USE BIEF
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)             :: NPFMAX,NPOIN2,NCOUCH
      DOUBLE PRECISION, INTENT(IN)    :: CFDEP
      INTEGER, INTENT(IN)             :: NPF(NPOIN2)
      DOUBLE PRECISION, INTENT(INOUT) :: MASBED0
      DOUBLE PRECISION, INTENT(IN)    :: EPAI(NPOIN2,NCOUCH)
      DOUBLE PRECISION, INTENT(IN)    :: EPAICO(NPOIN2), EPAINCO(NPOIN2)
      DOUBLE PRECISION, INTENT(IN)    :: VOLU2D(NPOIN2)
      DOUBLE PRECISION, INTENT(IN)    :: HDEP(NPOIN2)
      DOUBLE PRECISION, INTENT(IN)    :: CONC(NPOIN2,NCOUCH)
      DOUBLE PRECISION, INTENT(IN)    :: RHOS
      DOUBLE PRECISION, INTENT(INOUT) :: MASDEP
      TYPE(BIEF_OBJ), INTENT (INOUT)  :: ESOMT
      TYPE(BIEF_OBJ), INTENT (INOUT)  :: T2_01
      LOGICAL, INTENT(IN)             :: TASSE, GIBSON, MIXTE
      INTEGER, INTENT(IN)             :: DEBUG
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!     
          IF(DEBUG.GT.0) WRITE(LU,*) 'APPEL DE MASSED'
          CALL SEDI3D_MASSED(MASBED0,EPAI,CONC,HDEP,T2_01,NPOIN2,
     &                NPFMAX,
     &                NCOUCH,NPF,TASSE,GIBSON,RHOS,VOLU2D,
     &                CFDEP, EPAICO, EPAINCO, MIXTE)
          IF(DEBUG.GT.0) WRITE(LU,*) 'RETOUR DE MASSED'
          MASDEP = 0.D0
          CALL OS('X=0     ',X=ESOMT)
!         PRINT INITIAL MASS
          IF(LNG.EQ.1) THEN
            WRITE(LU,*) 'MASSE INITIALE DU LIT :',MASBED0
          ENDIF
          IF(LNG.EQ.2) THEN
            WRITE(LU,*) 'INITIAL MASS OF SEDIMENT BED :',MASBED0
          ENDIF
!
      RETURN
      END
