!                     ************************************
                      SUBROUTINE BED_UPDATE_ACTIVELAYER_HIRANO
!                     ************************************
!
!***********************************************************************
! SISYPHE   V7P3
!***********************************************************************
!
!     brief    MASS TRANSFER BETWEEN THE ACTIVE LAYER (FIRST LAYER) AND THE
!+             LAYER UNDERNEATH. THIS IS POSSIBLE AND NEEDED ONLY IF THERE
!+             IS NO CONSOLIDATION.
!+             THIS SUBROUTINE UPDATES MASS_SAND AND MASS_MUD USING THEIR MASS
!+             FLUX.
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
!
      USE BIEF
      USE DECLARATIONS_TELEMAC
      USE DECLARATIONS_SISYPHE, ONLY: MIN_SED_MASS_COMP,ES,
     &    MASS_SAND,MASS_MUD,MASS_SAND_TOT,MASS_MUD_TOT,MASS_MIX_TOT,
     &    RATIO_MUD_SAND,RATIO_MUD,RATIO_SAND,NSAND,NMUD,NPOIN,
     &    NOMBLAY,XKV,XMVS,CONC_MUD,ELAY
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!----------------------------------------------------------------------
!
      INTEGER IPOIN,ILAYER,ISAND,IMUD,IERR
      DOUBLE PRECISION THICK_TRANSFER,THICK_TRANSFER_TEMPO
      DOUBLE PRECISION TERM,DISCR,TOT
      DOUBLE PRECISION, ALLOCATABLE :: FLUX_MASS_MUD(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: FLUX_MASS_SAND(:,:,:)
!
!----------------------------------------------------------------------
!
!       ALLOCATE FLUXES
        IF(.NOT.ALLOCATED(FLUX_MASS_MUD)) THEN
          ALLOCATE(FLUX_MASS_MUD(NMUD,NOMBLAY,NPOIN), STAT=IERR)
          IF(IERR.NE.0) THEN
            IF(LNG.EQ.1) THEN
              WRITE(LU,*)'FLUX_MASS_MUD : ERREUR D''ALLOCATION',IERR
            ELSEIF(LNG.EQ.2) THEN
              WRITE(LU,*)'ERROR REALLOCATING FLUX_MASS_MUD:',IERR
            ENDIF
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
        IF(.NOT.ALLOCATED(FLUX_MASS_SAND)) THEN
          ALLOCATE(FLUX_MASS_SAND(NSAND,NOMBLAY,NPOIN), STAT=IERR)
          IF(IERR.NE.0) THEN
            IF(LNG.EQ.1) THEN
              WRITE(LU,*)'FLUX_MASS_SAND : ERREUR D''ALLOCATION',IERR
            ELSEIF(LNG.EQ.2) THEN
              WRITE(LU,*)'ERROR REALLOCATING FLUX_MASS_SAND:',IERR
            ENDIF
            CALL PLANTE(1)
            STOP
          ENDIF
        ENDIF
!
!
         DO IPOIN=1,NPOIN
            DO ILAYER=1,NOMBLAY
                DO IMUD=1,NMUD
                    FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)=0.D0
                ENDDO
                DO ISAND=1,NSAND
                    FLUX_MASS_SAND(ISAND,ILAYER,IPOIN)=0.D0
                ENDDO
            ENDDO
         ENDDO
!     COMPUTATION OF THE THICKNESS FOR TRANSFER
!
      DO IPOIN=1,NPOIN
        THICK_TRANSFER = ES(IPOIN,1)- ELAY%R(IPOIN)		
!       POSITIVE SIGNE MEANS DEPOSITION, NEGATIVE MEANS EROSION
!-----------------------------------------------------------------------
!      DEPOSITION CASE
!-----------------------------------------------------------------------
!
        IF(THICK_TRANSFER.GT.0.D0) THEN
! active layer too large: transfer of mass from layer 1 to layer 2
          IF(NMUD.GE.1) THEN
            DO IMUD = 1,NMUD
             FLUX_MASS_MUD(IMUD,1,IPOIN)= -(THICK_TRANSFER/ES(IPOIN,1))
     &       *MASS_MUD(IMUD,1,IPOIN)
              FLUX_MASS_MUD(IMUD,2,IPOIN)=+(THICK_TRANSFER/ES(IPOIN,1))
     &       *MASS_MUD(IMUD,1,IPOIN)
            ENDDO
          ENDIF
          IF(NSAND.GE.1) THEN
            DO ISAND=1,NSAND
            FLUX_MASS_SAND(ISAND,1,IPOIN)=-(THICK_TRANSFER/ES(IPOIN,1))
     &       *MASS_SAND(ISAND,1,IPOIN)
            FLUX_MASS_SAND(ISAND,2,IPOIN)=+(THICK_TRANSFER/ES(IPOIN,1))
     &       *MASS_SAND(ISAND,1,IPOIN)
            ENDDO
          ENDIF
!
!-----------------------------------------------------------------------
!        EROSION CASE
!-----------------------------------------------------------------------
!
        ELSE   !IF(THICK_TRANSFER.LE.0.D0) THEN
! active layer too small: transfer of mass needed from Sublayers to layer 1
          THICK_TRANSFER_TEMPO = ABS(THICK_TRANSFER)

          DO ILAYER =2,NOMBLAY
           IF (THICK_TRANSFER_TEMPO.GE.ES(IPOIN,ILAYER))THEN

            IF(NMUD.GE.1) THEN
                DO IMUD=1,NMUD
                    FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)=
     &               FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)-
     &               MASS_MUD(IMUD,ILAYER,IPOIN)
    
                    FLUX_MASS_MUD(IMUD,1,IPOIN)=
     &               FLUX_MASS_MUD(IMUD,1,IPOIN)+
     &               MASS_MUD(IMUD,ILAYER,IPOIN)
                ENDDO
            ENDIF
            IF(NSAND.GE.1) THEN
                DO ISAND=1,NSAND
                    FLUX_MASS_SAND(ISAND,ILAYER,IPOIN)=
     &               FLUX_MASS_SAND(ISAND,ILAYER,IPOIN)-
     &               MASS_SAND(ISAND,ILAYER,IPOIN)
    
                    FLUX_MASS_SAND(ISAND,1,IPOIN)=
     &               FLUX_MASS_SAND(ISAND,1,IPOIN)+
     &               MASS_SAND(ISAND,ILAYER,IPOIN)
                ENDDO
            ENDIF
            THICK_TRANSFER_TEMPO = THICK_TRANSFER_TEMPO
     &                              - ES(IPOIN,ILAYER)
          ELSE
            IF(NMUD.GE.1) THEN
                DO IMUD=1,NMUD
              FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)=
     &          FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)-
     &              (THICK_TRANSFER_TEMPO/ES(IPOIN,ILAYER))
     &              *MASS_MUD(IMUD,ILAYER,IPOIN)

              FLUX_MASS_MUD(IMUD,1,IPOIN)=
     &              FLUX_MASS_MUD(IMUD,1,IPOIN)+
     &              (THICK_TRANSFER_TEMPO/ES(IPOIN,ILAYER))
     &              *MASS_MUD(IMUD,ILAYER,IPOIN)     
                ENDDO
            ENDIF
            IF(NSAND.GE.1) THEN
      
                DO ISAND=1,NSAND
               FLUX_MASS_SAND(ISAND,ILAYER,IPOIN)=
     &              FLUX_MASS_SAND(ISAND,ILAYER,IPOIN)-
     &              (THICK_TRANSFER_TEMPO/ES(IPOIN,ILAYER))
     &              *MASS_SAND(ISAND,ILAYER,IPOIN)
     
               FLUX_MASS_SAND(ISAND,1,IPOIN)=
     &              FLUX_MASS_SAND(ISAND,1,IPOIN)+
     &              (THICK_TRANSFER_TEMPO/ES(IPOIN,ILAYER))
     &              *MASS_SAND(ISAND,ILAYER,IPOIN)
                ENDDO
            ENDIF         
          GOTO 100
          ENDIF
          ENDDO 
100   CONTINUE
        ENDIF
      ENDDO ! NPOIN
!
!-----------------------------------------------------------------------
!     TRANSFER OF MUD MASSES
!-----------------------------------------------------------------------
!
      IF(NMUD.GE.1) THEN
        DO IPOIN=1,NPOIN
         DO ILAYER=1,NOMBLAY
          DO IMUD = 1,NMUD
            MASS_MUD(IMUD,ILAYER,IPOIN)= MASS_MUD(IMUD,ILAYER,IPOIN)
     &                              +FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)
           ENDDO
          ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!     TRANSFER OF SAND MASSES
!-----------------------------------------------------------------------
!
!     NEGATIVE FLUX_MASS_SAND MEANS THE FLUX IS FROM UPPER LAYER TO LOWER LAYER
      IF(NSAND.GE.1) THEN
        DO IPOIN=1,NPOIN
            DO ILAYER=1,NOMBLAY
                DO ISAND=1,NSAND
                     MASS_SAND(ISAND,ILAYER,IPOIN)=
     &                      MASS_SAND(ISAND,ILAYER,IPOIN)
     &                      +FLUX_MASS_SAND(ISAND,ILAYER,IPOIN)
                ENDDO
          ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
      DEALLOCATE(FLUX_MASS_MUD)
      DEALLOCATE(FLUX_MASS_SAND)
!-----------------------------------------------------------------------
!
      RETURN
      END

