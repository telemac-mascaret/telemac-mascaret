!                    *************************
                     SUBROUTINE INIT_COMPO_COH
!                    *************************
!
     &(ES,CONC_VASE,CONC,NPOIN,NOMBLAY,NSICLA,AVAIL,AVA0)
!
!***********************************************************************
! SISYPHE   V6P2                                   21/07/2011
!***********************************************************************
!
!brief    INITIAL FRACTION DISTRIBUTION, STRATIFICATION,
!+                VARIATION IN SPACE.
!
!warning  USER SUBROUTINE; MUST BE CODED BY THE USER
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE INTERFACE_SISYPHE, EX_INIT_COMPO_COH=> INIT_COMPO_COH
      USE DECLARATIONS_SISYPHE, ONLY : NLAYMAX
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN)              :: NPOIN,NOMBLAY,NSICLA
      DOUBLE PRECISION, INTENT(INOUT)  :: ES(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(IN)     :: CONC_VASE(NOMBLAY)
      DOUBLE PRECISION,  INTENT(INOUT) :: CONC(NPOIN,NOMBLAY)
      DOUBLE PRECISION, INTENT(INOUT)  :: AVAIL(NPOIN,NOMBLAY,NSICLA)
      DOUBLE PRECISION, INTENT(IN)     :: AVA0(NSICLA)
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      DOUBLE PRECISION  EPAI_VASE(NLAYMAX),EPAI_SABLE(NLAYMAX)

      INTEGER I,J
!
!-----------------------------------------------------------------------
!    EXAMPLE FOR NOMBLAY = 10
!
!     EPAI_VASE(1)=0.0525D0
!     EPAI_VASE(2)=0.0385D0
!     EPAI_VASE(3)=0.03995D0
!     EPAI_VASE(4)=0.0437D0
!     EPAI_VASE(5)=0.0517D0
!     EPAI_VASE(6)=0.1259D0
!     EPAI_VASE(7)=0.4889D0
!     EPAI_VASE(8)=1.5071D0
!     EPAI_VASE(9)=0.86410D0
!     EPAI_VASE(9)=0.80D0
!
!
!     HERE A CONSTANT
!
!!      DO J= 1,NOMBLAY
!!        EPAI_VASE(J) = 0.1D0
!        IF(NSICLA.GT.1) THEN
!          EPAI_SABLE(J) = AVA0(1)/AVA0(2)*EPAI_VASE(J)
!        ENDDO
!!      ENDDO
!-----------------------------------------------------------------------
!
!     INITIALISING OF LAYER THICKNESS AND CONC
!

!     BY DEFAULT : UNIFORM BED COMPOSITION (KEY WORDS)
!     V6P3: IT WILL BE POSSIBLE TO HAVE A SPATIAL DISTRIBUTION OF THE BED CONC
!     V6P2: SO FAR THE MUD CONC IS CONSTANT PER LAYER
!     si mixte: calculer aussi les AVAI!
!
      DO I=1,NPOIN
        DO J= 1,NOMBLAY
!
          CONC(I,J) = CONC_VASE(J)
          ES(I,J)   = EPAI_VASE(J)
!          print*, I, J, CONC(I,J), ES(I,J)
!
!          IF(NSICLA.GT.1) THEN
!              ES(I,J)= ES(I,J) + EPAI_SABLE(J)
!              IF(ES(I,J).GE.1.D-6) THEN
! Class 1 is for sand, class 2 is mud
!                AVAIL(I,J,1)= EPAI_SABLE(J)/ES(I,J)
!                AVAIL(I,J,2)= EPAI_VASE(J)/ES(I,J)
!             ELSE
!              AVAIL(I,J,1)= 0.D0
!              AVAIL(I,J,2)= 0.D0
!             ENDIF
!
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      RETURN
      END

