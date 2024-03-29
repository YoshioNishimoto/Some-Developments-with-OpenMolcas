************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_VecDsk_GetLQ(QVec,l_QVec,LstQSP,nQSP,iV1,nV,mSym)
C
C     Purpose: extract elements corresponding to qualified columns of
C              vectors on disk.
C
      use ChoSwp, only: nnBstRSh
      Implicit Real*8 (a-h,o-z)
      Real*8  QVec(l_QVec)
      Integer LstQSP(nQSP)
      Integer iV1(mSym), nV(mSym)
#include "cholesky.fh"
#include "stdalloc.fh"

      Character(LEN=16), Parameter:: SecNam = 'Cho_VecDsk_GetLQ'

      Integer nVecTot(8)

      Integer, External:: Cho_P_LocalSP

      Real*8, Allocatable:: Scr(:)
      Integer, Allocatable:: iQuAB_2(:)

C     Check input.
C     ------------

      If (nQSP .lt. 1) Return

      If (mSym .lt. nSym) Then
         Call Cho_Quit('mSym<nSym in '//SecNam,104)
      End If

      nVT = nV(1)
      Do iSym = 2,nSym
         nVT = nVT + nV(iSym)
      End Do
      If (nVT .lt. 1) Return

      Do iSym = 1,nSym
         If (nV(iSym) .gt. 0) Then
            If (iV1(iSym) .lt. 1) Then
               Call Cho_Quit('iV1<1 in '//SecNam,104)
            End If
         End If
      End Do

C     Allocate memory for reading.
C     ----------------------------

      l_Scr = 0
      Do iSym = 1,nSym
         If (nV(iSym).gt.0 .and. nQual(iSym).gt.0) Then
            lDim = 0
            Do iQSP = 1,nQSP
               iShlAB = Cho_P_LocalSP(LstQSP(iQSP))
               lDim = lDim + nnBstRSh(iSym,iShlAB,1)
            End Do
            l_Scr = max(l_Scr,lDim)
         End If
      End Do
      Call mma_allocate(Scr,l_Scr,Label='Scr')

C     Allocate extra mapping from qualified to reduced set.
C     -----------------------------------------------------

      l_iQuAB_2 = nQual(1)
      Do iSym = 2,nSym
         l_iQuAB_2 = max(l_iQuAB_2,nQual(iSym))
      End Do
      Call mma_allocate(iQuAB_2,l_iQuAB_2,Label='iQuAB_2')

C     Extract in each symmetry block.
C     -------------------------------

      Call Cho_P_GetGV(nVecTot,nSym)

      jLoc = 2
      iLoc = 3
      iRedC = -1
      kOffQ = 0
      Do iSym = 1,nSym
         iRedQ = -2
         If (nV(iSym).gt.0 .and. nQual(iSym).gt.0) Then
            iV2 = iV1(iSym) + nV(iSym) - 1
            Do jV = iV1(iSym),iV2
               Call Cho_1VecRd_SP(Scr,l_Scr,jV,iSym,LstQSP,nQSP,iRedC,
     &                            iLoc)
               If (iRedQ .ne. iRedC) Then
                  Call Cho_SetQ2(iQuAB_2,LstQSP,nQSP,iSym,jLoc,iLoc)
                  iRedQ = iRedC
               End If
               kQ = kOffQ + nQual(iSym)*(jV-1)
               Do iQ = 1,nQual(iSym)
                  iAB = iQuAB_2(iQ)
                  QVec(kQ+iQ) = Scr(iAB)
               End Do
            End Do
         End If
         kOffQ = kOffQ + nQual(iSym)*nVecTot(iSym)
      End Do

C     Deallocations.
C     --------------

      Call mma_deallocate(iQuAB_2)
      Call mma_deallocate(Scr)

      End
