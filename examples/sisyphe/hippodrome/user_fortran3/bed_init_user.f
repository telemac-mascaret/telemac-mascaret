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
		    ESTRATUM(IPOIN,1)=100.D0
            DO ICLA=1,NSICLA
              RATIO_INIT(ICLA,1,IPOIN) = AVA0(ICLA)
			ENDDO
		ENDDO
      ELSE
       	  
!  USERS MUST DEFINE LAYER THICKNESS AND COMPOSITION FOR EACH STRATUM 
		DO IPOIN=1,NPOIN
		 DO ISTRAT=1,NOMBSTRAT
		    ESTRATUM(IPOIN,ISTRAT)=0.D0 ! par défaut les couches sont vides
            DO ICLA=1,NSICLA
              RATIO_INIT(ICLA,ISTRAT,IPOIN) = 1.D0/NSICLA 
			ENDDO			
		 ENDDO		
           IF(MESH%X%R(IPOIN).LT.1000.D0) THEN
            IF(MESH%Y%R(IPOIN).GT.0.D0) THEN ! on met 2 couches, la première étant plus petite que l'active layer
 		    ESTRATUM(IPOIN,1)=0.01D0
             RATIO_INIT(1,1,IPOIN) = 0.25D0
             RATIO_INIT(2,1,IPOIN) = 0.25D0			
             RATIO_INIT(3,1,IPOIN) = 0.25D0
             RATIO_INIT(4,1,IPOIN) = 0.25D0
 		    ESTRATUM(IPOIN,2)=0.2D0
             RATIO_INIT(1,2,IPOIN) = 0.1D0
             RATIO_INIT(2,2,IPOIN) = 0.2D0			
             RATIO_INIT(3,2,IPOIN) = 0.3D0
             RATIO_INIT(4,2,IPOIN) = 0.4D0				 
            ELSE		                     ! on met 3 couches
 		    ESTRATUM(IPOIN,1)=0.1D0
             RATIO_INIT(1,1,IPOIN) = 0.1D0
             RATIO_INIT(2,1,IPOIN) = 0.2D0			
             RATIO_INIT(3,1,IPOIN) = 0.3D0
             RATIO_INIT(4,1,IPOIN) = 0.4D0
 		    ESTRATUM(IPOIN,2)=0.5D0
             RATIO_INIT(1,2,IPOIN) = 0.25D0
             RATIO_INIT(2,2,IPOIN) = 0.25D0			
             RATIO_INIT(3,2,IPOIN) = 0.25D0
             RATIO_INIT(4,2,IPOIN) = 0.25D0	
 		    ESTRATUM(IPOIN,3)=0.1D0
             RATIO_INIT(1,3,IPOIN) = 0.4D0
             RATIO_INIT(2,3,IPOIN) = 0.3D0			
             RATIO_INIT(3,3,IPOIN) = 0.2D0
             RATIO_INIT(4,3,IPOIN) = 0.1D0			 
            ENDIF
           ENDIF		   
		ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
