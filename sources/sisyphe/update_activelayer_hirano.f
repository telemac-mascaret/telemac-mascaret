!                     ************************************
                      SUBROUTINE UPDATE_ACTIVELAYER_HIRANO
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
      USE DECLARATIONS_SISYPHE, ONLY: MIN_SED_MASS_COMP,ES,ESTRAT,
     &    MASS_SAND,MASS_MUD,MASS_SAND_TOT,MASS_MUD_TOT,MASS_MIX_TOT,
     &    RATIO_MUD_SAND,RATIO_MUD,RATIO_SAND,NSAND,NMUD,NPOIN,
     &    XKV,XMVS,CONC_MUD
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!----------------------------------------------------------------------
!
      INTEGER IPOIN,ILAYER,ISAND,IMUD
      DOUBLE PRECISION THICK_TRANSFER,TERM,DISCR
!     TEMPORARY: FLUX_MASS ARE ALLOCATABLE but they are used in
!     bed_consolidation_layer (the same? necessary to put in point_sisyphe?)
      DOUBLE PRECISION, ALLOCATABLE :: FLUX_MASS_MUD(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: FLUX_MASS_SAND(:,:)
!
!----------------------------------------------------------------------
!
      ALLOCATE(FLUX_MASS_MUD(NMUD,NPOIN))
      ALLOCATE(FLUX_MASS_SAND(NSAND,NPOIN))
!
!     INITIALIZATION
      DO IPOIN = 1,NPOIN
! delete this part!
        MASS_SAND_TOT(1,IPOIN) = 0.D0
        MASS_MUD_TOT(1,IPOIN) = 0.D0
        MASS_SAND_TOT(2,IPOIN) = 0.D0
        MASS_MUD_TOT(2,IPOIN) = 0.D0
! end delete
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
! not necessary for the first call:only from the 2nd call (maybe)
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
! end not necessary
! ES for the 1st layer is already computed in init_sediment or bed_update,
! so maybe not necessary here?
!----------------------------------------------------------------------
!     THICKNESS COMPUTATION
!----------------------------------------------------------------------
      DO IPOIN = 1,NPOIN
!     TERM REPRESENTS THE DIFFERENCE BETWEEN MUD VOLUME AND VOID VOLUME
        TERM=RATIO_MUD_SAND(1,IPOIN)/CONC_MUD(1)-
     &  (XKV(1)*(1.D0-RATIO_MUD_SAND(1,IPOIN)))/
     &  (XMVS*(1.D0-XKV(1)))
        DISCR=MAX(0.D0,TERM)
!     IF DISCR IS POSITIVE IT MEANS THAT MUD VOLUME IS LARGER THAN VOID VOLUME
!     IF DISCR IS NEGATIVE, THE VOID VOLUME IS NOT COMPLETELY FILLED BY MUD
        ES(IPOIN,1)=MASS_MIX_TOT(1,IPOIN)*
     &  ((1.D0-RATIO_MUD_SAND(1,IPOIN))/(XMVS*(1.D0-XKV(1)))
     &  + DISCR)
      ENDDO
! end not necessary ES
!
!     COMPUTATION OF THE THICKNESS FOR TRANSFER
!
      DO IPOIN=1,NPOIN
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
     &       *MASS_SAND_TOT(1,IPOIN)*RATIO_SAND(ISAND,1,IPOIN)
            ENDDO
          ENDIF
!
!-----------------------------------------------------------------------
!        EROSION CASE
!-----------------------------------------------------------------------
!
! why ELSEIF,not better ELSE (it includes thick_tr=0.d0
        ELSEIF(THICK_TRANSFER.LT.0.D0) THEN ! active layer too small: transfer of mass needed from layer 2
          THICK_TRANSFER=MIN(THICK_TRANSFER,ES(IPOIN,2)) ! if not enough sediment in layer 2 active layer thickness will be smaller than ESTRAT
! NOT CLEAR!! why do we do the min and then we recompute ES(IPOIN,2)????
!!           IF(MASS_SAND_TOT(2,IPOIN).GT.MIN_SED_MASS_COMP) THEN
!!              ES_PORO_SAND(IPOIN,2) = (MASS_SAND_TOT(2,IPOIN)*XKV(2)) ! thickness of second layer
!!     &                                /((1-XKV(2))*XMVS)
!!              ES_MUD_ONLY(IPOIN,2) = MASS_MUD_TOT(2,IPOIN)
!!     &                                /CONC_MUD((2,IPOIN)
!!           IF(ES_MUD_ONLY(IPOIN,2).GE.ES_PORO_SAND(IPOIN,2)) THEN
!!              ES(IPOIN,2) = MASS_SAND_TOT(2,IPOIN)/XMVS
!!     &                      +ES_MUD_ONLY(IPOIN,2)
!!           ELSE
!!               ES(IPOIN,2) = MASS_SAND_TOT(2,IPOIN)/((1-XKV(2))*XMVS)
!!           ENDIF
! end NOT CLEAR!!
          IF(NMUD.GE.1) THEN
            DO IMUD=1,NMUD
              FLUX_MASS_MUD(IMUD,IPOIN)=(THICK_TRANSFER/ES(IPOIN,2))
     &        *MASS_MUD_TOT(2,IPOIN)*RATIO_MUD(IMUD,2,IPOIN)
            ENDDO
          ENDIF
          IF(NSAND.GE.1) THEN
            DO ISAND=1,NSAND
               FLUX_MASS_SAND(ISAND,IPOIN)=
     &         (THICK_TRANSFER/ES(IPOIN,2))
     &         *MASS_SAND_TOT(2,IPOIN)*RATIO_SAND(ISAND,2,IPOIN)
            ENDDO
          ENDIF
        ENDIF
      ENDDO ! NPOIN
!
!-----------------------------------------------------------------------
!     TRANSFER OF MUD MASSES
!-----------------------------------------------------------------------
!
      IF(NMUD.GE.1) THEN
        DO IPOIN=1,NPOIN
          DO IMUD = 1,NMUD
            MASS_MUD(IMUD,1,IPOIN)= MASS_MUD(IMUD,1,IPOIN)
     &                              +FLUX_MASS_MUD(IMUD,IPOIN)
            MASS_MUD(IMUD,2,IPOIN)= MASS_MUD(IMUD,2,IPOIN)
     &                              -FLUX_MASS_MUD(IMUD,IPOIN)
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
      DEALLOCATE(FLUX_MASS_MUD)
      DEALLOCATE(FLUX_MASS_SAND)
!-----------------------------------------------------------------------
!
      RETURN
      END

