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
      SubRoutine Cho_P_SetRed_L()
C
C     Purpose: set next local reduced set. The next global reduced set
C              must be available at (global) index array location 2.
C
      Implicit None
#include "cholesky.fh"
#include "choptr.fh"
#include "choptr2.fh"
#include "choglob.fh"
#include "WrkSpc.fh"

      Character*14 SecNam
      Parameter (SecNam = 'Cho_P_SetRed_L')

      Integer irc, iC, nDim
      Integer i, j, k, i0, k0, l, ll
      Integer iSym, iSP, iShlAB
      Integer kOff

      Integer iL2G, mySP
      Integer IndRed, iiBstRSh, nnBstRSh
      Integer IndRed_G, iiBstRSh_G, nnBstRSh_G

      iL2G(i)=iWork(ip_iL2G-1+i)
      mySP(i)=iWork(ip_mySP-1+i)
      IndRed(i,j)=iWork(ip_IndRed-1+mmBstRT*(j-1)+i)
      iiBstRSh(i,j,k)=iWork(ip_iiBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)
      nnBstRSh(i,j,k)=iWork(ip_nnBstRSh-1+nSym*nnShl*(k-1)+nSym*(j-1)+i)
      IndRed_G(i,j)=iWork(ip_IndRed_G-1+mmBstRT_G*(j-1)+i)
      iiBstRSh_G(i,j,k)=
     &            iWork(ip_iiBstRSh_G-1+nSym*nnShl_G*(k-1)+nSym*(j-1)+i)
      nnBstRSh_G(i,j,k)=
     &            iWork(ip_nnBstRSh_G-1+nSym*nnShl_G*(k-1)+nSym*(j-1)+i)


C     Copy current local reduced set (at location 2) to location 3.
C     -------------------------------------------------------------

      Call Cho_X_RSCopy(irc,2,3)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
         Call Cho_Quit('Error in '//SecNam,104)
      End If

C     Re-initialize local reduced set indices at location 2.
C     ------------------------------------------------------

      nDim = nSym*nnShl
      Call Cho_iZero(iWork(ip_IndRed+mmBstRT),mmBstRT)
      Call Cho_iZero(iWork(ip_iiBstRSh+nDim),nDim)
      Call Cho_iZero(iWork(ip_nnBstRSh+nDim),nDim)
      Call Cho_iZero(iiBstR(1,2),nSym)
      Call Cho_iZero(nnBstR(1,2),nSym)
      nnBstRT(2) = 0

C     Set local nnBstRSh counter at location 2.
C     -----------------------------------------

      k0 = ip_nnBstRSh + nDim - 1
      Do iSP = 1,nnShl
         iShlAB = mySP(iSP)
         k = k0 + nSym*(iSP-1)
         Do iSym = 1,nSym
            iWork(k+iSym) = nnBstRSh_G(iSym,iShlAB,2)
         End Do
      End Do

C     Set remaining reduced set indices (excl. IndRed), location 2.
C     -------------------------------------------------------------

      Call Cho_SetRedInd(iWork(ip_iiBstRSh),iWork(ip_nnBstRSh),nSym,
     &                   nnShl,2)

C     Set local IndRed to point to local 1st reduced set.
C     ---------------------------------------------------

      kOff = ip_IndRed + mmBstRT - 1
      iC = 0
      Do iSym = 1,nSym
         Do iSP = 1,nnShl
            iShlAB = mySP(iSP)
            i0 = iiBstR_G(iSym,2) + iiBstRSh_G(iSym,iShlAB,2)
            k0 = iiBstR(iSym,3) + iiBstRSh(iSym,iSP,3)
            Do i = 1,nnBstRSh_G(iSym,iShlAB,2)
               j = IndRed_G(i0+i,2) ! addr in global rs1
               iC = iC + 1
               k = 0
               Do While (k .lt. nnBstRSh(iSym,iSP,3))
                  k = k + 1
                  ll = IndRed(k0+k,3)
                  l = iL2G(ll)
                  If (l .eq. j) Then
                     iWork(kOff+iC) = ll
                     k = nnBstRSh(iSym,iSP,3) ! break while loop
                  End If
               End Do
            End Do
         End Do
      End Do

      End
