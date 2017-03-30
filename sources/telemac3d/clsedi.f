!                    *****************
                     SUBROUTINE CLSEDI
!                    *****************
!
     &( ATABOF , BTABOF  , ATABOS , BTABOS  , SED_TRA     ,
     &  WC     , GRADZFX , GRADZFY, GRADZSX , GRADZSY,
     &  X      , Y       , Z      , HN      , DELTAR ,
     &  TOB    , DENSI   , TRA03  , EPAI    , CFDEP  ,
     &  CONC   , HDEP    , FLUER  , FLUDPT  , LITABF ,
     &  LITABS , KLOG    , NPOIN3 , NPOIN2  , NPLAN  ,
     &  NCOUCH , ITURBV  , DT     , RHO0    ,  RHOS  ,
     &  TOCD   , MPART   , TOCE   , UETCAR  , GRAV   ,
     &  DMOY    , CREF   , ZREF    , CF     ,
     &  AC     , KSPRATIO, ICR    , ICQ     , RUGOF  ,
     &  SETDEP , HMIN    , 
!WCS    , EPAICO  , EPAINCO,
!     &  MIXTE  , SEDNCO  , FLUDPTC, FLUDPTNC, FLUERC ,
!     &  FLUERNC, 
     & NCLASS   , ITRAC, ICLASS, TYPE_OF_SEDIMENT,SED_CO)
!
!***********************************************************************
! TELEMAC3D   V7P0                                   21/08/2010
!***********************************************************************
!
!brief    EXPRESSES THE BOUNDARY CONDITIONS FOR THE SEDIMENT,
!+                AT THE BOTTOM AND SURFACE (FOR COHESIVE SEDIMENT OR NOT).
!
!history  JACEK A. JANKOWSKI PINXIT
!+        **/03/99
!+
!+   FORTRAN95 VERSION
!
!history  CAMILLE LEQUETTE
!+        **/06/03
!+
!+
!
!history  C LE NORMANT (LNH)
!+        12/09/07
!+        V5P0
!+
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        13/07/2010
!+        V6P0
!+   Translation of French comments within the FORTRAN sources into
!+   English comments
!
!history  N.DURAND (HRW), S.E.BOURBAN (HRW)
!+        21/08/2010
!+        V6P0
!+   Creation of DOXYGEN tags for automated documentation and
!+   cross-referencing of the FORTRAN sources
!
!history  C. VILLARET & T. BENSON & D. KELLY (HR-WALLINGFORD)
!+        27/02/2014
!+        V7P0
!+   New developments in sediment merged on 25/02/2014.
!
!history  G. ANTOINE & M. JODEAU & J.M. HERVOUET (EDF - LNHE)
!+        13/10/2014
!+        V7P0
!+   New developments in sediment for mixed sediment transport
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!| AC             |-->| CRITICAL SHIELDS PARAMETER
!| ATABOF         |<->| FOR BOUNDARY CONDITION (BOTTOM)
!| ATABOS         |<->| FOR BOUNDARY CONDITION (SURFACE)
!| BTABOF         |<->| FOR BOUNDARY CONDITION (BOTTOM)
!| BTABOS         |<->| FOR BOUNDARY CONDITION (SURFACE)
!| CF             |-->| QUADRATIC FRICTION COEFFICIENT
!| CFDEP          |-->| MUD DEPOSITION CONCENTRATION (G/L)
!| CONC           |-->| MUD CONCENTRATION FOR EACH LAYER
!| CREF           |<->| EQUILIBRIUM NEAR-BED CONCENTRATION
!| DELTAR         |-->| DELTA RHO / RHO0 = (RHO-RHO0)/RHO0
!| DENSI          |<->| WATER DENSITY
!| DMOY           |-->| DIAMETRE MOYEN DES GRAINS
!| DT             |-->| HYDRODYNAMICS TIME STEP
!| EPAI           |<->| THICKNESS OF BOTTOM LAYERS IN
!|                |   | MATERIAL COORDINATES (EPAI=DZ/(1+IVIDE))
!| SETDEP         |-->| SETTLING SCHEME (0 or 1)
!| FLUDPT         |<->| IMPLICIT DEPOSITION FLUX
!| FLUER          |<->| EROSION FLUX FOR POINTS IN 2D
!| GRADZFX        |-->| GRADIENT-X OF BOTTOM
!| GRADZFY        |-->| GRADIENT-Y OF BOTTOM
!| GRADZSX        |-->| GRADIENT-X OF SURFACE
!| GRADZSY        |-->| GRADIENT-Y OF SURFACE
!| GRAV           |-->| ACCELERATION OF GRAVITY
!| HDEP           |<->| THICKNESS OF FRESH DEPOSIT (FLUID MUD LAYER)
!| HMIN           |-->| THRESHOLD FOR EROSION FLUXES ON TIDAL FLATS
!| HN             |-->| WATER DEPTH AT TIME N
!| ICR            |-->| FLAG FOR THE SKIN FRICTION OPTION
!| ICQ            |-->| FLAG FOR THE REFERENCE CONCENTRATION FORMULA
!| ITRAC          |-->| INDEX OF THE ACTIVE TRACER
!| ITURBV         |-->| VERTICAL TURBULENCE MODEL
!| KLOG           |-->| CONVENTION FOR LOGARITHMIC WALL
!| KSPRATIO       |-->| RELATION BETWEEN SKIN BED ROUGHNESS AND SEDIMENT DIAMETER
!| LITABF         |<->| FOR BOUNDARY CONDITION BOTTOM
!| LITABS         |<->| FOR BOUNDARY CONDITION SURFACE
!| MPART          |-->| EROSION COEFFICIENT (PARTHENIADES'S LAW)
!| NCOUCH         |-->| NUMBER OF LAYERS FOR THE COHESIVE MULTILAYER MODEL
!| NPLAN          |-->| NUMBER OF PLANES IN THE 3D MESH OF PRISMS
!| NPOIN2         |-->| NUMBER OF POINTS IN 2D
!| NPOIN3         |-->| NUMBER OF POINTS IN 3D
!| RHO0           |-->| WATER DENSITY (REFERENCE)
!| RHOS           |-->| MASSE VOLUMIQUE DU SEDIMENT
!| RUGOF          |-->| FRICTION COEFFICIENT
!| TA             |-->| SEDIMENT CONCENTRATION
!| TOB            |-->| BED SHEAR STRESS (TOTAL FRICTION)
!| TOCD           |-->| CRITICAL DEPOSITION SHEAR STRESS
!| TOCE           |-->| CRITICAL EROSION SHEAR STRESS
!| TRA03          |<->| WORK STRUCTURE FOR USER
!| UETCAR         |-->| SQUARE OF THE FRICTION VELOCITY
!| WC             |-->| SETTLING VELOCITY OF MUD
!| WCS            |-->| SETTLING VELOCITY OF SAND
!| X              |-->| FIRST NODE COORDINATE
!| Y              |-->| SECOND NODE COORDINATE
!| Z              |-->| THIRD NODE COORDINATE
!| ZREF           |-->| VERTICAL COORDINATE OF THE HARD BOTTOM
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE INTERFACE_TELEMAC3D, EX_CLSEDI => CLSEDI
      USE DECLARATIONS_SPECIAL
      USE DECLARATIONS_TELEMAC3D, ONLY: KARMAN,PRANDTL,FICT,KFROT
      USE DECLARATIONS_TELEMAC3D, ONLY: T2_01,T2_02,T2_03,T2_04,T2_05
      USE DECLARATIONS_TELEMAC3D, ONLY: DNUVIH
      IMPLICIT NONE
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER, INTENT(IN) :: NPOIN2,NPOIN3,KLOG,ICQ, ITRAC
      INTEGER, INTENT(IN) :: NCOUCH,ITURBV,NPLAN,ICR
!
      DOUBLE PRECISION, INTENT(INOUT) :: ATABOF(NPOIN2), BTABOF(NPOIN2)
      DOUBLE PRECISION, INTENT(INOUT) :: ATABOS(NPOIN2), BTABOS(NPOIN2)
!
      DOUBLE PRECISION, INTENT(IN)  :: X(NPOIN3), Y(NPOIN3), Z(NPOIN3)
      DOUBLE PRECISION, INTENT(IN)  :: SED_TRA(NPOIN3),CFDEP
      DOUBLE PRECISION, INTENT(IN)  :: WC(NPOIN3), DELTAR(NPOIN3)
!
      TYPE(BIEF_OBJ), INTENT(INOUT) :: TOB,CREF,ZREF,RUGOF
      TYPE(BIEF_OBJ), INTENT(IN)    :: DMOY,HN,CF
!
      DOUBLE PRECISION, INTENT(INOUT) :: EPAI(NPOIN2,NCOUCH)
      DOUBLE PRECISION, INTENT(IN)    :: CONC(NPOIN2,NCOUCH)
!
      DOUBLE PRECISION, INTENT(INOUT) :: DENSI(NPOIN2)
      DOUBLE PRECISION, INTENT(INOUT) :: TRA03(NPOIN2),UETCAR(NPOIN2)
      DOUBLE PRECISION, INTENT(INOUT) :: HDEP(NPOIN2),FLUER(NPOIN2)
      DOUBLE PRECISION, INTENT(INOUT) :: FLUDPT(NPOIN2)
!
      DOUBLE PRECISION, INTENT(IN) :: GRADZFX(NPOIN2),GRADZFY(NPOIN2)
      DOUBLE PRECISION, INTENT(IN) :: GRADZSX(NPOIN2),GRADZSY(NPOIN2)
!
      DOUBLE PRECISION, INTENT(IN) :: DT, RHO0, RHOS, HMIN
      DOUBLE PRECISION, INTENT(IN) :: TOCD, GRAV
      DOUBLE PRECISION, INTENT(IN) :: MPART, TOCE(NPOIN2,NCOUCH)
!
      INTEGER, INTENT(INOUT)       :: LITABF(NPOIN2), LITABS(NPOIN2)
      INTEGER, INTENT(IN)          :: SETDEP, TYPE_OF_SEDIMENT(NCLASS)
      INTEGER, INTENT(IN)          :: SED_CO,NCLASS,ICLASS
      DOUBLE PRECISION, INTENT(IN) :: AC, KSPRATIO
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      DOUBLE PRECISION KSP,A,ZERO,HCLIP,MU
      INTEGER IPOIN
!
!-----------------------------------------------------------------------
!
      ZERO = 1.D-6
!
!
! TOBW: wave induced shear stress
! FW, IKS

! calculation of friction from sisyphe routine
! new: no impact of waves (HOULE set to false) OK for you TOM?
! I am not sure if it always works without Nikuradse
! temporaeres 2D Feld
        CALL TOB_SISYPHE
     &    (TOB,T2_01, MU, RUGOF, KSP,T2_05,CF, T2_01,
!     &     CHESTR, UETCAR, CF_TEL,KS_TEL, CODE ,
     &     RUGOF, UETCAR,CF,RUGOF,'TELEMAC3D',
     &     KFROT, ICR, KSPRATIO,
! instead of HOULE COUROU, is that right?? HOULE is the keyword "effect of waves" in sisyphe
     &     .FALSE.,
     &     GRAV,RHO0,  RHOS, DNUVIH, KARMAN,ZERO,
! no roughness predictor: KSPRED=NO, flat bed IKS=1
     &     HMIN,HN,DMOY,T2_02,T2_03,T2_04,NPOIN2,.FALSE.,1)
!     &     HMIN,HN, DMOY, UNORM,UW, TW, NPOIN,KSPRED,IKS)


!      DO IPOIN=1,NPOIN2
!!       COMPUTES THE FLUID DENSITY
!        DENSI(IPOIN) = (DELTAR(IPOIN)+1.D0)*RHO0
!!       COMPUTES THE STRESS AT THE BOTTOM
!        TOB%R(IPOIN) = DENSI(IPOIN)*UETCAR(IPOIN)
!      ENDDO
!
!      IF(ICR.EQ.1) THEN
!
!        DO IPOIN=1,NPOIN2
!!         CORRECTION FOR SKIN FRICTION (SEE TOB_SISYPHE)
!          KSP=KSPRATIO *DMOY%R(IPOIN)
!          IF(CF%R(IPOIN).GT.ZERO.AND.HN%R(IPOIN).GT.KSP) THEN
!            HCLIP=MAX(HN%R(IPOIN),KSP)
!            A = 2.5D0*LOG(12.D0*HCLIP/KSP)
!            MU =2.D0/(A**2*CF%R(IPOIN))
!          ELSE
!            MU=0.D0
!          ENDIF
!          TOB%R(IPOIN) = MU* TOB%R(IPOIN)
!        ENDDO
!!
!      ENDIF
!
!      -----COMPUTES THE EXPLICIT EROSION FLUX-----
! Should be done in Sisyphe but called here, otherwise the FLUER is too late
!       CALL SUSPENSION_EROD()

! cohesive
      IF(TYPE_OF_SEDIMENT(ICLASS)==SED_CO) THEN
        CALL SUSPENSION_EROSION_COH()
! parameter list from sisyphe
!     &(TAUP,NPOIN,XMVS,PARTHENIADES,ZERO,
!     & FLUER,TOCE_VASE,NOMBLAY,DT,MS_VASE)
! MS_VASE: mass per layer MS_VASE(npoin2,nlayer)
     &(TOB,NPOIN2,RHOS,MPART,ZERO,
     & FLUER,TOCE,NCOUCH,DT,MS_VASE)
! new routine from Artelia?? because consolidation needs to be taken into account

      ELSE
! non-cohesive
        CALL SUSPENSION_EROSION()
! parameterlist from sisyphe:
!     &(TAUP,HN,FDM,FD90,AVA,NPOIN,CHARR,XMVE,XMVS,VCE,GRAV,HMIN,XWC,
!     & ZERO,ZREF,AC,FLUER,CSTAEQ,QSC,ICQ,U2D,V2D,CSRATIO,T14,DEBUG)
! dmoy ist die mittlere Korngroesse... es muss aber die aktuelle sein... wo bekomme 
! ich die denn her in T3D??
     &(TOB,HN,FDM,FD90,AVA,NPOIN2,CHARR,RHO0,RHOS,VCE,GRAV,HMIN,XWC,
     & ZERO,ZREF,AC,FLUER,CREF,QSC,ICQ,U2D,V2D,CSRATIO,T2_01,DEBUG)
      ENDIF


!      CALL SEDI3D_EROD !
!     &  (CONC,EPAI,FLUER,TOB,DENSI,
!     &  MPART,DT,NPOIN2,NCOUCH,TOCE,
!     &  HN,HMIN,MIXTE,EPAICO,
!     &  CFDEP,WC,HDEP,
!     &  NPOIN3,KSPRATIO,AC,RHOS,RHO0,
!     &  GRAV,DMOY,CREF,ZREF,CF,ICQ,RUGOF,Z,UETCAR,
!     &  SETDEP,EPAINCO,
!     &  KARMAN,PRANDTL,FICT,FLUERC,FLUERNC,NTRAC,ITRAC,
!     &  SEDCO,SEDNCO)
!
!      -----WRITES THE BOUNDARY CONDITIONS AT THE BOTTOM / SURFACE-----
!      -----                FOR THE SEDIMENT                      -----
!
      CALL FLUSED(ATABOF , BTABOF , ATABOS , BTABOS  ,
     &            LITABF , LITABS , SED_TRA     , WC      ,
     &            X      , Y      , Z      , HN%R    ,
     &            GRADZFX, GRADZFY, GRADZSX, GRADZSY ,
     &            TOB%R  , FLUDPT , FLUER  , TOCD    ,
     &            NPOIN3 , NPOIN2 , NPLAN  , KLOG    ,
     &            HMIN, SETDEP, ICLASS,NCLASS,TYPE_OF_SEDIMENT,SED_CO )
!
!-----------------------------------------------------------------------
!


        DO IPOIN=1,NPOIN2
          ATABOF(IPOIN) = -FLUDPT(IPOIN)
! may be too early?? because FLUER is calculated from sisyphe?!?
! -- asking TOM
          BTABOF(IPOIN) =  FLUER(IPOIN)
        ENDDO

      RETURN
      END
