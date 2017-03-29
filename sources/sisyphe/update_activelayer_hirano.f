
!                    *****************
             SUBROUTINE UPDATE_ACTIVELAYER_HIRANO
!                    *****************
!
!***********************************************************************
! SISYPHE   V7P3
!***********************************************************************
!
!     brief    MASS TRANSFER BETWEEN THE ACTIVE LAYER (FIRST LAYER) AND THE
!+             LAYER UNDERNEATH. THIS IS POSSIBLE AND NEEDED ONLY IF THERE
!+             IS NO CONSOLIDATION. 
!+             FOR THE CASE OF CONSOLIDATION, THE MASS TRANSFER BETWEEN LAYERS   
!+             IS ONLY CAUSED BY CONSOLIDATION AND THE ACTIVE LAYER IS VIRTUAL
!+             AND THEREFORE RECOMPUTED WHEN NEEDED FROM THE CONSOLIDATION LAYERS
!+
!+     history  CRWR/MDS/PT
!+        2017
!+        V7P3
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| XXX       |-->| XXX
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE
      USE DECLARATIONS_SPECIAL

!
      IMPLICIT NONE
      
      INTEGER LNG,LU
      COMMON/INFO/LNG,LU

      INTEGER IPOIN,ILAYER,ISAND
      DOUBLE PRECISION FLUX_MASS_MUD(NMUD,NPOIN2)
      DOUBLE PRECISION FLUX_MASS_SAND(NSAND,NPOIN2)
      DOUBLE PRECISION THICK_TRANSFER
!
!     CAREFUL !
!     all three tests below must be removed (and used to decide if this subroutine is called)
!
      IF(NLAYER.LT.2) THEN ! in this case no need to compute transfer of mass between layers, so no need to pass through this routine
         WRITE(LU,*) 'ONLY ONE LAYER : NO ACTIVE LAYER MODEL POSSIBLE'
         CALL PLANTE(1) 
         STOP
      ENDIF

      IF(NSAND+NMUD.LT.2) THEN
         WRITE(LU,*) 'ONLY ONE CLASS: NO ACTIVE LAYER MODEL POSSIBLE'
         CALL PLANTE(1)
         STOP
      ENDIF

      IF(CONSOLIDATION) THEN
        WRITE(LU,*) 'CONSOLIDATION: USE THE OTHER MODEL'
        CALL PLANTE(1)
        STOP
      ENDIF
!
      DO IPOIN = 1,NPOIN
        MASS_SAND_TOT(1,IPOIN) = 0.D0
        MASS_MUD_TOT(1,IPOIN) = 0.D0
        MASS_SAND_TOT(2,IPOIN) = 0.D0
        MASS_MUD_TOT(2,IPOIN) = 0.D0
        IF(NMUD.GE.1) THEN
          DO IMUD = 1,NMUD
            FLUX_MASS_MUD(IMUD,IPOIN) = 0.D0
          ENDDO
        ENDIF
        IF(NSAND.GE.1) THEN
          DO ISAND = 1,NSAND
            FLUX_MASS_SAND(ISAND,IPOIN) = 0.D0
          ENDDO
        ENDIF
      ENDDO

      IF(NSAND.GE.1) THEN
        DO IPOIN = 1,NPOIN
            DO ISAND = 1,NSAND
              MASS_SAND_TOT(1,IPOIN) = MASS_SAND_TOT(1,IPOIN)
     &        + MASS_SAND(ISAND,1,IPOIN)
              MASS_SAND_TOT(2,IPOIN) = MASS_SAND_TOT(2,IPOIN)
     &        + MASS_SAND(ISAND,2,IPOIN)
            ENDDO
        ENDDO
      ENDIF

      IF(NMUD.GE.1) THEN
        DO IPOIN = 1,NPOIN
          DO IMUD = 1,NMUD
            MASS_MUD_TOT(1,IPOIN) = MASS_MUD_TOT(1,IPOIN)
     &      + MASS_MUD(IMUD,1,IPOIN)
            MASS_MUD_TOT(2,IPOIN) = MASS_MUD_TOT(2,IPOIN)
     &      + MASS_MUD(IMUD,2,IPOIN)
          ENDDO
        ENDDO
      ENDIF

      DO IPOIN = 1,NPOIN
        IF(MASS_SAND_TOT(1,IPOIN).GT.1.D-9) THEN
          ES_PORO_SAND(IPOIN,1) = (MASS_SAND_TOT(1,IPOIN)*XKV) ! thickness of first (active) layer
     &                             /((1-XKV)*XMVS)
          ES_MUD_ONLY(IPOIN,1) = MASS_MUD_TOT(1,IPOIN)
     &                             /CONC_MUD((1,IPOIN)
        IF(ES_MUD_ONLY(IPOIN,1).GE.ES_PORO_SAND(IPOIN,1)) THEN
          ES(IPOIN,1) = MASS_SAND_TOT(1,IPOIN)/XMVS
     &                  +ES_MUD_ONLY(IPOIN,1)
        ELSE
	     ES(IPOIN,1) = MASS_SAND_TOT(1,IPOIN)/((1-XKV)*XMVS)
      ENDIF
!      
      THICK_TRANSFER=ESTRAT%R(IPOIN)-ES(IPOIN,1) ! negative sign means from layer 1 to layer 2 
!-----------------------------------------------------------------------
!      DEPOSITION CASE
!-----------------------------------------------------------------------
!          
        IF(THICK_TRANSFER.GT.0.D0) THEN ! active layer too large: transfer of mass needed to layer 2

          IF(NMUD.GE.1) THEN
            DO IMUD = 1,NMUD
              FLUX_MASS_MUD(IMUD,IPOIN)=(THICK_TRANSFER/ES(IPOIN,1))
     &       *MASS_MUD_TOT(1,IPOIN)*RATIO_MUD(IMUD,1,IPOIN)  
            ENDDO
          ENDIF
          IF(NSAND.GE.1) THEN
            DO ISAND=1,NSAND
              FLUX_MASS_SAND(ISAND,IPOIN)=(THICK_TRANSFER/ES(IPOIN,1))
     &       *MASS_SAND_TOT(1,IPOIN) *RATIO_SAND(ISAND,1,IPOIN)  
            ENDDO
          ENDIF
!      
!-----------------------------------------------------------------------
!      EROSION CASE
!-----------------------------------------------------------------------
!    
          ELSEIF(THICK_TRANSFER.LT.0.D0) ! active layer too small: transfer of mass needed from layer 2
             THICK_TRANSFER=MIN(THICK_TRANSFER,ES(IPOIN,2)) ! if not enough sediment in layer 2 active layer thickness will be smaller than ESTRAT
             IF(MASS_SAND_TOT(2,IPOIN).GT.1.D-9) THEN
                ES_PORO_SAND(IPOIN,2) = (MASS_SAND_TOT(2,IPOIN)*XKV) ! thickness of second layer
     &                                  /((1-XKV)*XMVS)
                ES_MUD_ONLY(IPOIN,2) = MASS_MUD_TOT(2,IPOIN)
     &                                  /CONC_MUD((2,IPOIN)
             IF(ES_MUD_ONLY(IPOIN,2).GE.ES_PORO_SAND(IPOIN,2)) THEN
	        ES(IPOIN,2) = MASS_SAND_TOT(2,IPOIN)/XMVS
     &                        +ES_MUD_ONLY(IPOIN,2)
             ELSE
	         ES(IPOIN,2) = MASS_SAND_TOT(2,IPOIN)/((1-XKV)*XMVS)
             ENDIF
             IF(NMUD.GE.1) THEN
               DO IMUD=1,NMUD
                 FLUX_MASS_MUD(IMUD,IPOIN)=(THICK_TRANSFER/ES(IPOIN,2))
     &          *MASS_MUD_TOT(2,IPOIN) *RATIO_MUD(IMUD,2,IPOIN)
               ENDDO
             ENDIF
             IF(NSAND.GE.1) THEN
               DO ISAND=1,NSAND
                  FLUX_MASS_SAND(ISAND,IPOIN)=
     &            (THICK_TRANSFER/ES(IPOIN,2))
     &            *MASS_SAND_TOT(2,IPOIN)*RATIO_SAND(ISAND,2,IPOIN)      
               ENDDO
             ENDIF
          ENDIF ! end test 2nd layer large enough        
        ENDIF ! end test 1st layer large enough 
      ENDDO
!      
!-----------------------------------------------------------------------
!     TRANSFER OF MUD MASSES
!-----------------------------------------------------------------------
!
      IF(NMUD.GE.1) THEN
        DO IPOIN=1,NPOIN
          DO IMUD = 1,NMUD
             MASS_MUD(IMUD,1,IPOIN)= MASS_MUD(IMUD,1,IPOIN)
     &                               +FLUX_MASS_MUD(IMUD,IPOIN)
             MASS_MUD(IMUD,2,IPOIN)= MASS_MUD(IMUD,2,IPOIN)
     &                               -FLUX_MASS_MUD(IMUD,IPOIN)
          ENDDO
        ENDDO
      ENDIF
!      
!-----------------------------------------------------------------------
!     TRANSFER OF SAND MASSES
!-----------------------------------------------------------------------
!
      IF(NSAND.GE.1) THEN
        DO IPOIN=1,NPOIN
          DO ISAND=1,NSAND
            MASS_SAND(ISAND,1,IPOIN)= MASS_SAND(ISAND,1,IPOIN)
     &                                +FLUX_MASS_SAND(ISAND,IPOIN)
            MASS_SAND(ISAND,2,IPOIN)= MASS_SAND(ISAND,2,IPOIN)
     &                                -FLUX_MASS_SAND(ISAND,IPOIN)
          ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE UPDATE_ACTIVELAYER

