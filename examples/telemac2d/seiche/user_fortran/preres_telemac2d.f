!                    ***************************
                     SUBROUTINE PRERES_TELEMAC2D
!                    ***************************
     &    (IMP,LEO)
!
!
!***********************************************************************
! TELEMAC2D   V7P1
!***********************************************************************
!
!brief    PREPARES THE VARIABLES WHICH WILL BE WRITTEN TO
!+                THE RESULTS FILE OR TO THE LISTING.
!
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_TELEMAC2D
      USE INTERFACE_TELEMAC2D
!
      USE DECLARATIONS_SPECIAL
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      LOGICAL ,INTENT(INOUT)::IMP,LEO
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      LOGICAL DEJA1,DEJA2
!
      INTEGER LTT,N,IMAX,I,II,JJ
!
      DOUBLE PRECISION HHH,XMAX,ARG1,COEF,HHPLG
      DOUBLE PRECISION, PARAMETER:: EPSS=1.E-10
      DOUBLE PRECISION GPRDTIME,LPRDTIME,RESTE
!
      INTRINSIC MAX,SQRT,CEILING
!
      DOUBLE PRECISION P_DMAX,P_DMIN
      EXTERNAL         P_DMAX,P_DMIN
!
!-----------------------------------------------------------------------
!
      DATA DEJA1/.FALSE./
      DATA DEJA2/.FALSE./
      SAVE DEJA1,DEJA2
!!-----------------------------------------------------------------------
!
!     THE OUTPUT VARIABLES ARE BUILT ONLY IF NECESSARY, HENCE THE
!     FOLLOWING TESTS, WHICH MUST BE THE SAME AS IN BIEF_DESIMP (BIEF LIBRARY)
!
      IMP=.FALSE.
      LEO=.FALSE.
!     THIS WILL TRIGGER THE OUTPUT OF LAST TIMESTEP
!     BUT NOT WITH PARAMETER ESTIMATION (LISPRD WOULD STAY AT 1
!     FOR FURTHER COMPUTATIONS)
!      IF(LT.EQ.NIT.AND.ESTIME(1:1).EQ.' ') THEN
!        IMP=.FALSE.
!        LEO=.FALSE.
!      ENDIF
!     Always write the intial conditions
      IF(LT.EQ.0) THEN
        IMP=.TRUE.
        LEO=.TRUE.
        COMPLEO=0
      ELSE
        IF(EQUA(1:15).NE.'SAINT-VENANT VF') THEN
!         FEM
          LTT=(LT/LISPRD)*LISPRD
          IF(LT.EQ.LTT.AND.LT.GE.PTINIL) IMP=.TRUE.
          LTT=(LT/LEOPRD)*LEOPRD
          IF(LT.EQ.LTT.AND.LT.GE.PTINIG) LEO=.TRUE.
!         FOR GRAPHICAL OUTPUTS
          IF(LEO)COMPLEO=COMPLEO+1
        ELSE
!         FVM
          GPRDTIME=LEOPRD*DTINI
          LPRDTIME=LISPRD*DTINI
          IF(GPRDTIME.LT.EPSS.OR.LPRDTIME.LT.EPSS)THEN
            CALL PLANTE(1)
            STOP
          ENDIF
          IF(LT.GE.PTINIG)THEN
!           GRAPHIC OUTPUT
            LTT=CEILING(AT/GPRDTIME)
            RESTE=(LTT*GPRDTIME-AT)/GPRDTIME
            IF(RESTE.LT.EPSS.OR.ABS(RESTE-1.D0).LT.EPSS)THEN
!                                   CASE WHERE RESTE=1
              LEO=.TRUE.
              COMPLEO=COMPLEO+1
            ENDIF
          ENDIF
          IF(LT.GT.PTINIL)THEN
!           LISTING OUTPUT
            LTT=CEILING(AT/LPRDTIME)
            RESTE=(LTT*LPRDTIME-AT)/LPRDTIME
            IF(RESTE.LT.EPSS.OR.ABS(RESTE-1.D0).LT.EPSS)IMP=.TRUE.
          ENDIF
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!
!
!     PAS D'IMPRESSION, PAS DE SORTIE SUR FICHIER, ON RESSORT
      IF(.NOT.(LEO.OR.IMP)) GO TO 1000
!
!
!=======================================================================
! CALCUL DE LA CELERITE (MISE DANS FU, VOIR LE BLOC VARSOR)
!=======================================================================
!
      IF((LEO.AND.SORLEO(3)).OR.(IMP.AND.SORIMP(3))) THEN
        DO N=1,NPOIN
          FU%R(N) = SQRT ( GRAV * MAX(H%R(N),0.D0) )
        ENDDO
      ENDIF
!
!=======================================================================
! CALCUL DE LA SURFACE LIBRE (= H + ZF, MISE DANS FV)
!=======================================================================
!
      IF((LEO.AND.SORLEO(5)).OR.(IMP.AND.SORIMP(5))) THEN
        CALL OS( 'X=Y+Z   ' , FV , H , ZF , 0.D0 )
      ENDIF
!
!=======================================================================
! CALCUL DU NOMBRE DE FROUDE
!=======================================================================
!
      IF((LEO.AND.SORLEO(7)).OR.(IMP.AND.SORIMP(7))) THEN
        DO N=1,NPOIN
          HHPLG = MAX( H%R(N) , 1.D-8 )
          T2%R(N) = SQRT (( U%R(N)**2 + V%R(N)**2 ) / ( HHPLG*GRAV ))
        ENDDO
      ENDIF
!
!=======================================================================
! CALCUL DU DEBIT SCALAIRE
!=======================================================================
!
      IF((LEO.AND.SORLEO(8)).OR.(IMP.AND.SORIMP(8))) THEN
        DO N=1,NPOIN
          T3%R(N) = SQRT (U%R(N)**2 + V%R(N)**2) * H%R(N)
        ENDDO
      ENDIF
!
!=======================================================================
! CALCUL DU DEBIT VECTORIEL , COMPOSANTE SUIVANT X
!=======================================================================
!
      IF((LEO.AND.SORLEO(13)).OR.(IMP.AND.SORIMP(13))) THEN
        CALL OS( 'X=YZ    ' , T4 , H , U , HHPLG )
      ENDIF
!
!=======================================================================
! CALCUL DU DEBIT VECTORIEL , COMPOSANTE SUIVANT Y
!=======================================================================
!
      IF((LEO.AND.SORLEO(14)).OR.(IMP.AND.SORIMP(14))) THEN
        CALL OS( 'X=YZ    ' , T5 , H , V , HHPLG )
      ENDIF
!
!=======================================================================
! CALCUL DE LA VITESSE SCALAIRE
!=======================================================================
!
      IF((LEO.AND.SORLEO(15)).OR.(IMP.AND.SORIMP(15))) THEN
        CALL OS( 'X=N(Y,Z)' , T6 , U , V , HHPLG )
      ENDIF
!
!=======================================================================
! CALCUL DU NOMBRE DE COURANT
!=======================================================================
!
      IF((LEO.AND.SORLEO(22)).OR.(IMP.AND.SORIMP(22))) THEN
!                             IELM
        CALL CFLPSI(T9,U,V,DT,11,MESH,MSK,MASKEL)
        CALL MAXI(XMAX,IMAX,T9%R,NPOIN)
        IF (LNG.EQ.1) WRITE(LU,78) XMAX
        IF (LNG.EQ.2) WRITE(LU,79) XMAX
78      FORMAT(1X,'PRERES : NOMBRE DE COURANT MAXIMUM :',G16.7)
79      FORMAT(1X,'PRERES: MAXIMUM COURANT NUMBER: ',G16.7)
      ENDIF
!=======================================================================
! CALCUL DE LA HAUTEUR EXACTE
!=======================================================================
!
      COEF = SQRT(1.D0*GRAV)
      IF((LEO.AND.SORLEO(23)).OR.(IMP.AND.SORIMP(23))) THEN
        DO N=1,NPOIN
        ARG1 = COEF-(X(N)-5.98)/(2.D0/MAX(AT,DT))
        PRIVE%ADR(1)%P%R(N) = MAX(0.D0,ARG1)
        PRIVE%ADR(1)%P%R(N) = 4.D0*PRIVE%ADR(1)%P%R(N)**2/9.D0/9.81
        PRIVE%ADR(1)%P%R(N) = MIN(4.D0,PRIVE%ADR(1)%P%R(N))
        ENDDO
      ENDIF
!
!=======================================================================
!
1000  CONTINUE
      RETURN
      END

