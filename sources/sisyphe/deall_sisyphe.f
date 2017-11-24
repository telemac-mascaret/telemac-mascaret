!                    **************************
                     SUBROUTINE DEALL_SISYPHE
!                    **************************
!
!
!***********************************************************************
! SISYPHE   V7P1                                   19/05/2016
!***********************************************************************
!
!brief    Memory deallocation of structures, aliases, blocks...
!
!Author  R-S MOURADI (LNHE)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
      USE BIEF
      USE DECLARATIONS_SPECIAL
      USE DECLARATIONS_SISYPHE
!
      IMPLICIT NONE

!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
      INTEGER :: I
!
!+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
!
!
      CALL OUTBIEF(MESH)
!     Deallocate the mesh structure
      CALL DEALMESH(MESH)

      NULLIFY(IKLE)
      NULLIFY(X)
      NULLIFY(Y)
      NULLIFY(NELEM)
      NULLIFY(NELMAX)
      NULLIFY(NPTFR)
      NULLIFY(NPTFRX)
      NULLIFY(TYPELM)
      NULLIFY(NPOIN)
      NULLIFY(NPMAX)
      NULLIFY(MXPTVS)
      NULLIFY(LV)
!


!     Deallocating all the blocks first
      CALL BIEF_DEALLOBJ(E)
      CALL BIEF_DEALLOBJ(ECPL)
      CALL BIEF_DEALLOBJ(Z)
      CALL BIEF_DEALLOBJ(DEL_Z)
      CALL BIEF_DEALLOBJ(ZF_C)
      CALL BIEF_DEALLOBJ(ZF_S)
      CALL BIEF_DEALLOBJ(ESOMT)
      CALL BIEF_DEALLOBJ(EMAX)
      CALL BIEF_DEALLOBJ(QU)
      CALL BIEF_DEALLOBJ(QV)
      CALL BIEF_DEALLOBJ(DEL_QU)
      CALL BIEF_DEALLOBJ(DEL_QV)
      CALL BIEF_DEALLOBJ(DEL_UW)
      CALL BIEF_DEALLOBJ(Q)
      CALL BIEF_DEALLOBJ(QS)
      CALL BIEF_DEALLOBJ(QSX)
      CALL BIEF_DEALLOBJ(QSY)
      CALL BIEF_DEALLOBJ(QS_C)
      CALL BIEF_DEALLOBJ(QSXC)
      CALL BIEF_DEALLOBJ(QSYC)
      CALL BIEF_DEALLOBJ(QS_S)
      CALL BIEF_DEALLOBJ(QSXS)
      CALL BIEF_DEALLOBJ(QSYS)
      CALL BIEF_DEALLOBJ(HN)
      CALL BIEF_DEALLOBJ(HCLIP)
      CALL BIEF_DEALLOBJ(U2D)
      CALL BIEF_DEALLOBJ(V2D)
      CALL BIEF_DEALLOBJ(UNORM)
      CALL BIEF_DEALLOBJ(HCPL)
      CALL BIEF_DEALLOBJ(EBOR)
      CALL BIEF_DEALLOBJ(QBOR)
      CALL BIEF_DEALLOBJ(Q2BOR)
      CALL BIEF_DEALLOBJ(FLBOR)
      CALL BIEF_DEALLOBJ(ZF)
      CALL BIEF_DEALLOBJ(ZR)
      CALL BIEF_DEALLOBJ(ZREF)
      CALL BIEF_DEALLOBJ(VOLU2D)
      CALL BIEF_DEALLOBJ(V2DPAR)
      CALL BIEF_DEALLOBJ(UNSV2D)
      CALL BIEF_DEALLOBJ(CHESTR)
      CALL BIEF_DEALLOBJ(CALFA)
      CALL BIEF_DEALLOBJ(SALFA)
      CALL BIEF_DEALLOBJ(S)
      CALL BIEF_DEALLOBJ(MASKPT)
      CALL BIEF_DEALLOBJ(MASKTR)
      CALL BIEF_DEALLOBJ(MASKB)
      CALL BIEF_DEALLOBJ(MASKEL)
      CALL BIEF_DEALLOBJ(MSKTMP)
      CALL BIEF_DEALLOBJ(W1)
      CALL BIEF_DEALLOBJ(THETAW)
      CALL BIEF_DEALLOBJ(FW)
      CALL BIEF_DEALLOBJ(UW)
      CALL BIEF_DEALLOBJ(HW)
      CALL BIEF_DEALLOBJ(TW)
      CALL BIEF_DEALLOBJ(INDIC)
      CALL BIEF_DEALLOBJ(IFAMAS)
      CALL BIEF_DEALLOBJ(IT1)
      CALL BIEF_DEALLOBJ(IT2)
      CALL BIEF_DEALLOBJ(IT3)
      CALL BIEF_DEALLOBJ(IT4)
      CALL BIEF_DEALLOBJ(LIEBOR)
      CALL BIEF_DEALLOBJ(LIMTEC)
      CALL BIEF_DEALLOBJ(COEFPN)
      CALL BIEF_DEALLOBJ(NUMLIQ)
      CALL BIEF_DEALLOBJ(TOB)
      CALL BIEF_DEALLOBJ(CF)
      CALL BIEF_DEALLOBJ(TOBW)
      CALL BIEF_DEALLOBJ(MU)
      CALL BIEF_DEALLOBJ(KS)
      CALL BIEF_DEALLOBJ(KSP)
      CALL BIEF_DEALLOBJ(KSR)
      CALL BIEF_DEALLOBJ(DZF_GF)
      CALL BIEF_DEALLOBJ(ACLADM)
      CALL BIEF_DEALLOBJ(UNLADM)
      CALL BIEF_DEALLOBJ(NLAYER)
      CALL BIEF_DEALLOBJ(HIDING)
      CALL BIEF_DEALLOBJ(ELAY)
      CALL BIEF_DEALLOBJ(ESTRAT)
      CALL BIEF_DEALLOBJ(FLUDP)
      CALL BIEF_DEALLOBJ(FLUDPT)
      CALL BIEF_DEALLOBJ(FLUER)
      CALL BIEF_DEALLOBJ(FLUERT)
      CALL BIEF_DEALLOBJ(CS)
      CALL BIEF_DEALLOBJ(CST)
      CALL BIEF_DEALLOBJ(CTILD)
      CALL BIEF_DEALLOBJ(CSTAEQ)
      CALL BIEF_DEALLOBJ(CSRATIO)
      CALL BIEF_DEALLOBJ(CBOR)
      CALL BIEF_DEALLOBJ(CSGL)
      CALL BIEF_DEALLOBJ(UCONV)
      CALL BIEF_DEALLOBJ(VCONV)
      CALL BIEF_DEALLOBJ(HPROP)
      CALL BIEF_DEALLOBJ(DISP)
      CALL BIEF_DEALLOBJ(DISP_C)
      CALL BIEF_DEALLOBJ(AFBOR)
      CALL BIEF_DEALLOBJ(BFBOR)
      CALL BIEF_DEALLOBJ(FLBOR_SIS)
      CALL BIEF_DEALLOBJ(FLBORTRA)
      CALL BIEF_DEALLOBJ(LICBOR)
      CALL BIEF_DEALLOBJ(LIHBOR)
      CALL BIEF_DEALLOBJ(LIMPRO)
      CALL BIEF_DEALLOBJ(LIMDIF)
      CALL BIEF_DEALLOBJ(BOUNDARY_COLOUR)
      CALL BIEF_DEALLOBJ(CLT)
      CALL BIEF_DEALLOBJ(CLU)
      CALL BIEF_DEALLOBJ(CLV)
      CALL BIEF_DEALLOBJ(TE1)
      CALL BIEF_DEALLOBJ(TE2)
      CALL BIEF_DEALLOBJ(TE3)
      CALL BIEF_DEALLOBJ(KX)
      CALL BIEF_DEALLOBJ(KY)
      CALL BIEF_DEALLOBJ(KZ)
      CALL BIEF_DEALLOBJ(BREACH)
      CALL BIEF_DEALLOBJ(FLUER_VASE)
      CALL BIEF_DEALLOBJ(TOCE_MIXTE)
      CALL BIEF_DEALLOBJ(MS_SABLE)
      CALL BIEF_DEALLOBJ(MS_VASE)

      CALL BIEF_DEALLOBJ(MBOR)
      CALL BIEF_DEALLOBJ(AM1_S)
      CALL BIEF_DEALLOBJ(AM2_S)

      CALL BIEF_DEALLOBJ(MASK)
      CALL BIEF_DEALLOBJ(TB)
      CALL BIEF_DEALLOBJ(TB2)
      CALL BIEF_DEALLOBJ(PRIVE)
      CALL BIEF_DEALLOBJ(VARCL)
      CALL BIEF_DEALLOBJ(VARHYD)
      !CALL BIEF_DEALLOBJ(VARSOR)

      IF(ALLOCATED(PRO_F)) THEN
        DEALLOCATE(PRO_F)
      ENDIF
      IF(ALLOCATED(PRO_D)) THEN
        DEALLOCATE(PRO_D)
      ENDIF
      DEALLOCATE(AVAIL)
      DEALLOCATE(ES)
      DEALLOCATE(ES_VASE)
      DEALLOCATE(ES_SABLE)

      !CALL BIEF_DEALLOBJ(AVAI)
      !CALL BIEF_DEALLOBJ(LAYTHI)
      !CALL BIEF_DEALLOBJ(LAYCONC)
      CALL BIEF_DEALLOBJ(QSCL)
      CALL BIEF_DEALLOBJ(QSCL_C)
      CALL BIEF_DEALLOBJ(QSCLXC)
      CALL BIEF_DEALLOBJ(QSCLYC)
      CALL BIEF_DEALLOBJ(QSCL_S)
      CALL BIEF_DEALLOBJ(QSCLXS)
      CALL BIEF_DEALLOBJ(QSCLYS)
      CALL BIEF_DEALLOBJ(ZFCL)
      CALL BIEF_DEALLOBJ(ZFCL_C)
      CALL BIEF_DEALLOBJ(ZFCL_S)
      CALL BIEF_DEALLOBJ(ZFCL_MS)
      CALL BIEF_DEALLOBJ(MPM_ARAY)
      CALL BIEF_DEALLOBJ(FLULIM)
      CALL BIEF_DEALLOBJ(FLBCLA)

      DEALLOCATE(IVIDE)
      DEALLOCATE(CONC)

      IF(ALLOCATED(PRO_MAX)) THEN
        DEALLOCATE(PRO_MAX)
      ENDIF

      IF(ALLOCATED(CTRLSC)) THEN
        DEALLOCATE(CTRLSC)
      ENDIF

      DEALLOCATE(OKCGL)

      DEALLOCATE(CBOR_CLASSE)

      DEALLOCATE(SOLDIS)

      CALL BIEF_DEALLOBJ(T1)
      CALL BIEF_DEALLOBJ(T2)
      CALL BIEF_DEALLOBJ(T3)
      CALL BIEF_DEALLOBJ(T4)
      CALL BIEF_DEALLOBJ(T5)
      CALL BIEF_DEALLOBJ(T6)
      CALL BIEF_DEALLOBJ(T7)
      CALL BIEF_DEALLOBJ(T8)
      CALL BIEF_DEALLOBJ(T9)
      CALL BIEF_DEALLOBJ(T10)
      CALL BIEF_DEALLOBJ(T11)
      CALL BIEF_DEALLOBJ(T12)
      CALL BIEF_DEALLOBJ(T13)
      CALL BIEF_DEALLOBJ(T14)
      !CALL BIEF_DEALLOBJ(IKLE)

      IF(ALLOCATED(CHAIN)) THEN
        DEALLOCATE(CHAIN)
      END IF
!
!  to check: all variables of new model deallocated?
!
      DEALLOCATE(RATIO_MUD_SAND)
      DEALLOCATE(MASS_MUD)
      DEALLOCATE(MASS_SAND)
      DEALLOCATE(MASS_MUD_TOT)
      DEALLOCATE(MASS_SAND_TOT)
      DEALLOCATE(MASS_MIX_TOT)
      DEALLOCATE(ES_PORO_SAND)
      DEALLOCATE(ES_MUD_ONLY)
      DEALLOCATE(RATIO_SAND)
      DEALLOCATE(RATIO_MUD)
      DEALLOCATE(TOCE_SAND)
      DEALLOCATE(TOCE_MUD)
!
! Resetting variable
      INIT_FLUXPR = .TRUE.
      DEJA_RFC = .FALSE.
      DEJA_FLUSEC = .FALSE.
      OLD_METHOD_FLUSEC = .FALSE.
      DEJA_FLUSEC2 = .FALSE.

!
!-----------------------------------------------------------------------
!
!=======================================================================
!
! WRITES OUT TO LISTING :
!
!      IF(LISTIN) THEN
!        IF(LNG.EQ.1) WRITE(LU,22)
!        IF(LNG.EQ.2) WRITE(LU,23)
!      ENDIF
!22    FORMAT(1X,///,21X,'****************************************',/,
!     &21X,              '* FIN DE LA DEALLOCATION DE LA MEMOIRE  : *',/,
!     &21X,              '****************************************',/)
!23    FORMAT(1X,///,21X,'*************************************',/,
!     &21X,              '*    END OF MEMORY ORGANIZATION:    *',/,
!     &21X,              '*************************************',/)
!
!-----------------------------------------------------------------------
!
      RETURN
      END SUBROUTINE DEALL_SISYPHE
