!            **********************************
             SUBROUTINE BED_CONSOLIDATION_LAYER
!            **********************************
!
!            
!***********************************************************************
! SISYPHE   V7P3                                             28/03/2017
!***********************************************************************
!
!brief    COMPUTES BED CONSOLIDATION;
!+
!
!history  R. WALTHER (ARTELIA), J. FONTAINE (EDF-LNHE)
!+        28/03/2017
!+        V7P3
!+  Creation of the subroutine.
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| CA             |-->| BOTTOM CONCENTRATION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
      USE BIEF
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
      INTEGER IPOIN,ILAYER,ISAND
      DOUBLE PRECISION FLUX_MASS_MUD(NPOIN,NLAYER)
      DOUBLE PRECISION FLUX_MASS_SAND(NPOIN,NLAYER,NSAND)
!
!======================================================================
!     le flux de consolidation est dM/dt=a*M
!     je pense que la matrice A(ILAYER) est deja declaree dans sisyphe sous un autre nom
!     
      IF(NMUD.EQ.0)THEN
        WRITE(LU,*)'CONSOLIDATION WITH NO MUD??'
        STOP
      ENDIF
      IF(A(NLAYER).NE.0.D0)THEN
        WRITE(LU,*)'LE TAUX DE TRANFERT POUR LA CONSOLIDATION DOIT', 
     &             ' ETRE NUL POUR LA DERNIERE COUCHE'
        STOP
      ENDIF
!     calcul des flux de transfer de vase
      DO IPOIN = 1,NPOIN
        DO ILAYER = 1,NLAYER
          IF(MASS_MUD_TOT(ILAYER,IPOIN).GE.MIN_SED_MASS_COMP)THEN
!           
            FLUX_MASS_MUD_TOT(ILAYER,IPOIN) =
     &           A(ILAYER)*MASS_MUD_TOT(ILAYER,IPOIN)*DT
!           
            DO IMUD = 1,NMUD
!             
              FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)=
     &             MIN(MASS_MUD(IMUD,ILAYER,IPOIN),
     &             FLUX_MASS_MUD_TOT(ILAYER,IPOIN)
     &             *RATIO_MUD(IMUD,ILAYER,IPOIN))
!             
            ENDDO
          ELSE
            DO IMUD = 1,NMUD
              FLUX_MASS_MUD(IMUD,ILAYER,IPOIN) = 0.D0
            ENDDO
          ENDIF
!         on calcule le transfert total de vase
          FLUX_MASS_MUD_TOT(ILAYER,IPOIN) = 0.D0
          DO IMUD = 1,NMUD
!           
            FLUX_MASS_MUD_TOT(ILAYER,IPOIN) =
     &           FLUX_MASS_MUD_TOT(ILAYER,IPOIN)
     &           +FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)
!           
          ENDDO
        ENDDO
      ENDDO
!     calcul des flux de transfer de sable (le sable accompagne la vase dans les flux)
      IF(NSAND.GE.1)THEN
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NLAYER
            DO ISAND = 1,NSAND
              IF(MASS_MUD(IMUD,ILAYER,IPOIN).GE.MIN_SED_MASS_COMP)THEN
!               
                FLUX_MASS_SAND(ISAND,ILAYER,IPOIN) =
     &               FLUX_MASS_MUD_TOT(ILAYER,IPOIN)/
     &               MASS_MUD_TOT(ILAYER,IPOIN)
     &               *MASS_SAND(ISAND,ILAYER,IPOIN)
!               
                FLUX_MASS_SAND(ISAND,ILAYER,IPOIN) =
     &               MAX(0.D0,FLUX_MASS_SAND(ISAND,ILAYER,IPOIN))
!               
                FLUX_MASS_SAND(ISAND,ILAYER,IPOIN) =
     &               MIN(FLUX_MASS_SAND(ISAND,ILAYER,IPOIN),
     &               MASS_SAND(ISAND,ILAYER,IPOIN))
!               
              ELSE
                FLUX_MASS_SAND(ISAND,ILAYER,IPOIN) = 0.D0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!     transfert des masses de vase
      DO IPOIN = 1,NPOIN2
        DO ILAYER = 1,NLAYER
          DO IMUD = 1,NMUD
            IF(ILAYER.EQ.NLAYER)THEN
!             
              MASS_MUD(IMUD,ILAYER,IPOIN) =
     &             MASS_MUD(IMUD,ILAYER,IPOIN)
     &             +FLUX_MASS_MUD(IMUD,ILAYER-1,IPOIN)
!             
            ELSEIF(ILAYER.EQ.1)THEN
!             
              MASS_MUD(IMUD,ILAYER,IPOIN) =
     &             MASS_MUD(IMUD,ILAYER,IPOIN)
     &             -FLUX_MASS_MUD(IMUD,IPOIN,ILAYER,IPOIN)
!             
            ELSE
!             
              MASS_MUD(IMUD,ILAYER,IPOIN) =
     &             MASS_MUD(IMUD,ILAYER,IPOIN)
     &             +FLUX_MASS_MUD(IMUD,ILAYER-1,IPOIN)
     &             -FLUX_MASS_MUD(IMUD,ILAYER,IPOIN)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
!     transfert des masses de sable
      IF(NSAND.GE.1)THEN
        DO IPOIN = 1,NPOIN
          DO ILAYER = 1,NLAYER
            DO ISAND = 1,NSAND
              IF(ILAYER.EQ.NLAYER)THEN
!               
                MASS_SAND(ISAND,ILAYER,IPOIN) =
     &               MASS_SAND(ISAND,ILAYER,IPOIN)
     &               +FLUX_MASS_SAND(ISAND,ILAYER-1,IPOIN)
!               
              ELSEIF(ILAYER.EQ.1)THEN
!               
                MASS_SAND(ISAND,ILAYER,IPOIN) =
     &               MASS_SAND(ISAND,ILAYER,IPOIN)
     &               -FLUX_MASS_SAND(ISAND,ILAYER,IPOIN)
!               
              ELSE
!               
                MASS_SAND(ISAND,ILAYER,IPOIN) =
     &               MASS_SAND(ISAND,ILAYER,IPOIN)
     &               +FLUX_MASS_SAND(ISAND,ILAYER-1,IPOIN)
     &               -FLUX_MASS_SAND(ISAND,ILAYER,IPOIN)
!               
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!     
      RETURN
      END SUBROUTINE BED_CONSOLIDATION_LAYER
