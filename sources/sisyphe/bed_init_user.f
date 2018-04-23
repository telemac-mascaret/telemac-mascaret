!                    ************************
                     SUBROUTINE BED_INIT_USER
!                    ************************
!
     &(ESTRATUM,RATIO_INIT)
!
!***********************************************************************
! 
      USE BIEF
      USE DECLARATIONS_SPECIAL
      USE DECLARATIONS_SISYPHE
      IMPLICIT NONE
!
      DOUBLE PRECISION,INTENT(INOUT):: ESTRATUM(NPOIN,NOMBSTRAT)
      DOUBLE PRECISION,INTENT(INOUT):: 
     &     RATIO_INIT(NSICLA,NOMBSTRAT,NPOIN)
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER            :: IPOIN,ICLA,ISTRAT
!     
!======================================================================!
!======================================================================!
!                               PROGRAM                                !
!======================================================================!
!======================================================================!
!
! DEFAULT CASE : NO STRATIFICATION = ONLY ONE STRATUM

      IF (NOMBSTRAT.EQ.1) THEN
        DO IPOIN=1,NPOIN
            ESTRATUM(IPOIN,1)=100.D0 ! User can change the thickness of sediment here (replaces noerod.f)
            XKV(1,IPOIN)=XKV0(1)            
            DO ICLA=1,NSICLA
              RATIO_INIT(ICLA,1,IPOIN) = AVA0(ICLA)
            ENDDO
        ENDDO
      ELSE
          
!  USERS MUST DEFINE LAYER THICKNESS AND COMPOSITION FOR EACH STRATUM 
        DO IPOIN=1,NPOIN
         DO ISTRAT=1,NOMBSTRAT
            ESTRATUM(IPOIN,ISTRAT)=0.D0 ! par défaut les couches sont vides
            XKV(ISTRAT+1,IPOIN)=XKV0(ISTRAT)
            DO ICLA=1,NSICLA
              RATIO_INIT(ICLA,ISTRAT,IPOIN) = 1.D0/NSICLA 
            ENDDO           
         ENDDO
         XKV(1,IPOIN)=XKV0(1)        
         ! user must specify non-voids layers below 
         !...        
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
