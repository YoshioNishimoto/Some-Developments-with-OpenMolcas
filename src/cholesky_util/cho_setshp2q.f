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
      SubRoutine Cho_SetShP2Q(irc,iLoc,iShlAB,nAB)
C
C     Purpose: set mapping from shell pair iShlAB to qualified
C              columns within current reduced set (stored at location
C              iLoc = 2 or 3).
C              If a non-zero code (irc) is returned, nothing has been
C              set!!
C
      use ChoArr, only: iSP2F, nBstSh, iShP2Q
      use ChoSwp, only: iQuAB, IndRed
#if defined (_DEBUGPRINT_)
      use ChoSwp, only: IndRSh
#endif
#include "implicit.fh"
      Integer nAB(8)
#include "cholesky.fh"

#if defined (_DEBUGPRINT_)
      Character*12 SecNam
      Parameter (SecNam = 'Cho_SetShP2Q')
#endif

C     Check allocations.
C     ------------------

      Call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.True.)
      If (iShlA .eq. iShlB) Then
         NumAB = nBstSh(iShlA)*(nBstSh(iShlA)+1)/2
      Else
         NumAB = nBstSh(iShlA)*nBstSh(iShlB)
      End If
      lTst = 2*NumAB
      l_iShP2Q=0
      If (Allocated(iShP2Q)) l_iShP2Q=SIZE(iShP2Q)
      If (l_iShP2Q.lt.1 .or. l_iShP2Q.lt.lTst) Then
         irc = 102
         Return
      End If

C     Check iLoc.
C     -----------

      If (iLoc.lt.2 .or. iLoc.gt.3) Then
         irc = 104
         Return
      End If

C     Set mapping array.
C     iShP2Q(1,AB) = index among qualified, symmetry reduced.
C     iShP2Q(2,AB) = symmetry block.
C     Zeros are returned if the element AB is not qualified.
C     -------------------------------------------------------

      iShP2Q(:,1:NumAB)=0

      Do iSym = 1,nSym
         Do lAB = 1,nAB(iSym)
            iAB = iQuAB(iOffQ(iSym)+lAB,iSym) ! addr in current rs
            jAB = IndRed(iAB,iLoc)            ! addr in 1st rs
            kAB = IndRed(jAB,1)               ! addr in full shell pair
#if defined (_DEBUGPRINT_)
            nErr = 0
            If (IndRSh(jAB).ne.iSP2F(iShlAB)) Then
               Write(Lupri,*) SecNam,': inconsistent shell pairs!'
               Write(Lupri,*) SecNam,': from input: ',iSP2F(iShlAB),
     &                        '  from IndRsh: ',IndRSh(jAB)
               nErr = nErr + 1
            End If
            If (kAB.lt.1 .or. kAB.gt.NumAB) Then
               Write(Lupri,*) SecNam,': shell pair address error!'
               Write(Lupri,*) SecNam,': kAB = ',kAB
               Write(Lupri,*) SecNam,': min and max allowed: 1 ',
     &                        NumAB
               nErr = nErr + 1
            End If
            If (nErr .ne. 0) Then
               Write(Lupri,*) SecNam,': Shell A, B, AB: ',
     &                        iShlA,iShlB,iShlAB
               Write(Lupri,*) SecNam,': iLoc: ',iLoc
               Write(Lupri,*) SecNam,': symmetry block: ',iSym
               Write(Lupri,*) SecNam,': red. set address, ',
     &                        'first and current: ',
     &                        jAB,iiBstR(iSym,iLoc)+iAB
               Call Cho_Quit('Error detected in '//SecNam,104)
            End If
#endif
            iShP2Q(1,kAB)   = lAB
            iShP2Q(2,kAB) = iSym
         End Do
      End Do

C     Set return code 0: all ok!
C     --------------------------

      irc = 0

      End
