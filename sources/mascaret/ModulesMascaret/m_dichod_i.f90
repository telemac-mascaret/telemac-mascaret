!== Copyright (C) 2000-2017 EDF-CEREMA ==
!
!   This file is part of MASCARET.
!
!   MASCARET is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   MASCARET is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with MASCARET.  If not, see <http://www.gnu.org/licenses/>
!

module M_DICHOD_I
! *********************************************************************
! PROGICIEL : MASCARET         N. GOUTAL     F. MAUREL
!
! VERSION : 8.1.3                EDF-CEREMA
! *********************************************************************

   interface

   subroutine DICHOD ( &
        JG           , &
        JD           , &
        SURF         , &
        NOEUD        , &
        SGEOD        , &
        NMLARG       , &
        Erreur         &
                     )

!***********************************************************************
! FONCTION :        RECHERCHE PAR DICHOTOMIE DES BORNES DE L'INTERVALLE
!                   CONTENANT LA SURFACE MOUILLEE SURF SUR LE MAILLAGE
!                   INITIAL
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .___________.____.____.______________________________________________.
! !    NOM    !TYPE!MODE!                   ROLE                       !
! !___________!____!____!______________________________________________!
! !  JG       !  I !  R ! BORNE GAUCHE DE L'INTERVALLE CONTENANT SURF  !
! !  JD       !  I !  R ! BORNE DROITE DE L'INTERVALLE CONTENANT SURF  !
! !  SURF     !  R !  D ! SURFACE MOUILLEE                             !
! !  NOEUD    !  I !  D ! NOEUD CONSIDERE DU MAILLAGE                  !
! !  SGEOD    ! TR !  D ! SURFACE MOUILLEE PLANIMETREE AUX INTERFACES  !
! !  NMLARG   !  I !    ! Nombre de pas des sections                   !
! !___________!____!____!______________________________________________!
!
!                             VARIABLES LOCALES
! .___________.____.____.______________________________________________.
! !  SG       !  R !  A ! SURFACE MOUILLE POUR LA BORNE GAUCHE         !
! !  JMIL     !  I !  A ! MILIEU DE L'INTERVALLE [JG,JD]               !
! !  NMAX     !  I !  A ! NOMBRE MAXIMAL D'ITERATION POUR LA DICHOTOMIE!
! !  NT       !  I !  A ! NOMBRE D'ITERATIONS DEJA EFFECTUEES          !
! !  JGP1     !  I !  A ! JG+1 (<JD TANT QU ON ITERE)                  !
! !___________!____!____!______________________________________________!
!
!     TYPE : I (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : D (DONNEE NON MODIFIEE), R (RESULTAT), M (DONNEE MODIFIEE)
!            A (AUXILIAIRE MODIFIE)
!
!***********************************************************************
!  SGEOD fait partie d'une structure de donnees STRUCTURE_SECTIONS

   !============================= Declarations ===========================
   !.. Modules importes ..
   !----------------------
   use M_PRECISION
   use M_MESSAGE_C           ! Messages d'erreur
   use M_ERREUR_T            ! ERREUR
   use M_TRAITER_ERREUR_I    ! Traitement des erreurs

   !.. Declarations explicites ..
   !-----------------------------
   implicit none

   !.. Arguments ..
   !---------------
   integer     ,                   intent(  out) :: JG,JD
   real(DOUBLE),                   intent(in)    :: SURF
   integer     ,                   intent(in)    :: NOEUD
   ! 1ere dimension IM1, 2nde dimension NMLARG
   real(DOUBLE), dimension(:,:)  , intent(in)    :: SGEOD
   integer     ,                   intent(in)    :: NMLARG
   Type (ERREUR_T)               , intent(inout) :: Erreur

   end subroutine DICHOD

   end interface

end module M_DICHOD_I
