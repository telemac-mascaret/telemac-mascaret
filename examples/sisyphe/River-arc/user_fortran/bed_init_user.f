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
      INTEGER NPMAX2
      INTEGER ND
      INTEGER NG	
      PARAMETER (NPMAX2=200)
      DOUBLE PRECISION XD(NPMAX2),YD(NPMAX2)
      DOUBLE PRECISION XG(NPMAX2),YG(NPMAX2)	  
!	  
!======================================================================!
!======================================================================!
!                               PROGRAM                                !
!======================================================================!
!======================================================================!


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
! DEFAULT CASE : NO STRATIFICATION = ONLY ONE STRATUM

      IF (NOMBSTRAT.EQ.1) THEN
		DO IPOIN=1,NPOIN
		    ESTRATUM(IPOIN,1)=1000.D0 ! User can change the thickness of sediment here (replaces noerod.f)
        IF (INPOLY(MESH%X%R(IPOIN),MESH%Y%R(IPOIN),XD,YD,ND).OR.
     &      INPOLY(MESH%X%R(IPOIN),MESH%Y%R(IPOIN),XG,YG,NG)) THEN
          ESTRATUM(IPOIN,1)=0.D0
        ENDIF
	
      IF (ZF%R(IPOIN).GT.454.5D0) THEN
!
      RATIO_INIT(1,1,IPOIN) = 0.065
      RATIO_INIT(2,1,IPOIN) = 0.045
      RATIO_INIT(3,1,IPOIN) = 0.145
      RATIO_INIT(4,1,IPOIN) = 0.285
      RATIO_INIT(5,1,IPOIN) = 0.345
      RATIO_INIT(6,1,IPOIN) = 0.115
!
      ELSEIF(ZF%R(IPOIN).LT.453.8D0) THEN
!
      RATIO_INIT(1,1,IPOIN) = 0.065D0
      RATIO_INIT(2,1,IPOIN) = 0.
      RATIO_INIT(3,1,IPOIN) = 0.
      RATIO_INIT(4,1,IPOIN) = 0.
      RATIO_INIT(5,1,IPOIN) = 0.
      RATIO_INIT(6,1,IPOIN) = 0.935D0
!
      ELSEIF(ZF%R(IPOIN).LE.454.5D0.AND.ZF%R(IPOIN).GE.453.8D0) THEN
!      DGRA = DMAX + (zf%R(J)-453.8D0)*(0.022D0-DMAX)/(454.5D0-453.8D0)
!
      RATIO_INIT(1,1,IPOIN) = 0.065
      RATIO_INIT(2,1,IPOIN) = 
     & 0.045*(ZF%R(IPOIN)-453.8D0)/(454.5D0-453.8D0)
      RATIO_INIT(3,1,IPOIN) = 
     & 0.145*(ZF%R(IPOIN)-453.8D0)/(454.5D0-453.8D0)
      RATIO_INIT(4,1,IPOIN) =
     & 0.285*(ZF%R(IPOIN)-453.8D0)/(454.5D0-453.8D0)
      RATIO_INIT(5,1,IPOIN) = 
     & 0.345*(ZF%R(IPOIN)-453.8D0)/(454.5D0-453.8D0)
      RATIO_INIT(6,1,IPOIN) = 1.-RATIO_INIT(1,1,IPOIN)
     &                          -RATIO_INIT(2,1,IPOIN)
     &                          -RATIO_INIT(3,1,IPOIN)
     &                          -RATIO_INIT(4,1,IPOIN)
     &                          -RATIO_INIT(5,1,IPOIN)	 
!
      ENDIF
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
         ! user must specify non-voids layers below	
         !...		 
		ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
