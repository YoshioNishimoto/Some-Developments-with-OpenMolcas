************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      subroutine Cho_X_Dealloc(irc)

      use Definitions, only: iwp
      use ChoArr, only: iSOShl, iBasSh, nBasSh, nBstSh, iSP2F, iAtomShl,
     &                  iShlSO, iRS2F, IntMap, iScr, nDimRS, iL2G,
     &                  iShP2RS, iShP2Q, iQL2G, LQ_Tot, iSimRI

      use ChoSwp, only: iQuAB, iQuAB_L, iQuAB_Hidden, iQuAB_L_Hidden,
     &                  nnBstRSh_Hidden, nnBstRSh,
     &                  nnBstRSh_L_Hidden, nnBstRSh_G,
     &                  iiBstRSh_Hidden, iiBstRSh,
     &                  iiBstRSh_L_Hidden, iiBstRSh_G,
     &                    IndRSh_Hidden,   IndRSh,
     &                    IndRSh_G_Hidden,   IndRSh_G,
     &                    InfRed_Hidden,   InfRed,
     &                    InfRed_G_Hidden,   InfRed_G,
     &                    InfVec_Hidden,   InfVec,
     &                    InfVec_G_Hidden,   InfVec_G,
     &                    IndRed_Hidden,   IndRed,
     &                    IndRed_G_Hidden,   IndRed_G,
     &                    InfVec_Bak
      use stdalloc, only: mma_deallocate
C
C     T.B. Pedersen, July 2004.
C
C     Purpose: deallocate ALL index arrays for the Cholesky utility.
C              On exit, irc=0 signals sucessful completion.
C
      implicit none
      integer(kind=iwp) :: irc
#include "chpari.fh"

      character(len=13), parameter :: SecNam = 'Cho_X_Dealloc'

      irc = 0

C     Deallocate.
C     -----------

      If (Allocated(InfRed_Hidden))
     &    Call mma_deallocate(InfRed_Hidden)
      If (Associated(InfRed)) InfRed=>Null()

      If (Allocated(InfVec_Hidden))
     &    Call mma_deallocate(InfVec_Hidden)
      If (Associated(InfVec)) InfVec=>Null()

      If (Allocated(IndRed_Hidden))
     &    Call mma_deallocate(IndRed_Hidden)
      If (Associated(IndRed)) IndRed=>Null()

      If (Allocated(IndRSh_Hidden))
     &    Call mma_deallocate(IndRSh_Hidden)
      If (Associated(IndRSh)) IndRSh=>Null()

      If (Allocated(iScr)) Call mma_deallocate(iScr)

      If (Allocated(iiBstRSh_Hidden))
     &    Call mma_deallocate(iiBstRSh_Hidden)
      If (Associated(iiBstRSh)) iiBstRSh=>Null()

      If (Allocated(nnBstRSh_Hidden))
     &    Call mma_deallocate(nnBstRSh_Hidden)
      If (Associated(nnBstRSh)) nnBstRSh=>Null()

      If (Allocated(IntMap)) Call mma_deallocate(IntMap)

      If (Allocated(nDimRS)) Call mma_deallocate(nDimRS)

      If (Allocated(iRS2F)) Call mma_deallocate(iRS2F)

      If (Allocated(iSOShl)) Call mma_deallocate(iSOShl)

      If (Allocated(iShlSO)) Call mma_deallocate(iShlSO)

      If (Allocated(iQuAB_Hidden)) Call mma_deallocate(iQuAB_Hidden)
      If (Associated(iQuAB)) iQuAB => Null()

      If (Allocated(iBasSh)) Call mma_deallocate(iBasSh)

      If (Allocated(nBasSh)) Call mma_deallocate(nBasSh)

      If (Allocated(nBstSh)) Call mma_deallocate(nBstSh)

      If (Allocated(iAtomShl)) Call mma_deallocate(iAtomShl)

      If (Allocated(iSP2F)) Call mma_deallocate(iSP2F)

      If (Allocated(iShP2RS)) Call mma_deallocate(iShP2RS)

      If (Allocated(iShP2Q )) Call mma_deallocate(iShP2Q )

      If (Allocated(iQuAB_L_Hidden)) Call mma_deallocate(iQuAB_L_Hidden)
      If (Associated(iQuAB_L)) iQuAB_L => Null()

      If (Allocated(iQL2G )) Call mma_deallocate(iQL2G )

      If (Allocated(LQ_Tot)) Call mma_deallocate(LQ_Tot)

      If (Allocated(InfVec_Bak)) Call mma_deallocate(InfVec_Bak)

      If (Allocated(iSimRI)) Call mma_deallocate(iSimRI)

      If (Allocated(InfVec_G_Hidden))
     &    Call mma_deallocate(InfVec_G_Hidden)
      If (Associated(InfVec_G)) InfVec_G=>Null()

      If (Allocated(IndRed_G_Hidden))
     &    Call mma_deallocate(IndRed_G_Hidden)
      If (Associated(IndRed_G)) IndRed_G=>Null()

      If (Allocated(InfRed_G_Hidden))
     &    Call mma_deallocate(InfRed_G_Hidden)
      If (Associated(InfRed_G)) InfRed_G=>Null()

      If (Allocated(IndRSh_G_Hidden))
     &    Call mma_deallocate(IndRSh_G_Hidden)
      If (Associated(IndRSh_G)) IndRSh_G=>Null()

      If (Allocated(iiBstRSh_L_Hidden))
     &    Call mma_deallocate(iiBstRSh_L_Hidden)
      If (Associated(iiBstRSh_G)) iiBstRSh_G=>Null()

      If (Allocated(nnBstRSh_L_Hidden))
     &    Call mma_deallocate(nnBstRSh_L_Hidden)
      If (Associated(nnBstRSh_G)) nnBstRSh_G=>Null()

      If (Allocated(iL2G)) Call mma_deallocate(iL2G)

      Return
      End
