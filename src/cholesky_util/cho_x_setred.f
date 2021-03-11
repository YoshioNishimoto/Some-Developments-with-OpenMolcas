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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  Cho_X_SetRed
*
*> @brief
*>   Read and set index arrays for reduced set \p iRed at location \p iLoc
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Reads information for reduced set \p iRed (= ``1``, ``2``, ..., \c MaxRed)
*> and sets up the index arrays
*>
*> - \p nnBstRT(iLoc)      &rarr; stored in cholesky.fh
*> - \p nnBstR(:,iLoc)     &rarr; stored in cholesky.fh
*> - \p iiBstR(:,iLoc)     &rarr; stored in cholesky.fh
*> - \p nnBstRSh(:,:,iLoc) &rarr; accesible via \c ip_nnBstRSh in choptr.fh
*> - \p iiBstRSh(:,:,iLoc) &rarr; accesible via \c ip_iiBstRSh in choptr.fh
*> - \p IndRed(:,iLoc)     &rarr; accesible via \c ip_IndRed in choptr.fh
*>
*> On succesful completion, \p irc = ``0`` is returned.
*> Note that the only allowed \p iLoc values are ``2`` and ``3``; any other
*> value gives rise to error code \p irc = ``1`` and nothing is done!
*> If \p iRed is out of bounds, \p irc = ``2`` is returned and nothing is done!
*>
*> @note
*> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
*>
*> @param[out] irc  return code
*> @param[in]  iLoc location in index arrays
*> @param[in]  iRed reduced set on disk
************************************************************************
      Subroutine Cho_X_SetRed(irc,iLoc,iRed)
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      If (iLoc.eq.2 .or. iLoc.eq.3) Then
         If (iRed.lt.1 .or. iRed.gt.MaxRed) Then
            irc = 2
         Else
            kOff1 = ip_nnBstRSh + nSym*nnShl*(iLoc - 1)
            kOff2 = ip_IndRed   + nnBstRT(1)*(iLoc - 1)
            Call Cho_GetRed(iWork(ip_InfRed),iWork(kOff1),
     &                      iWork(kOff2),iWork(ip_IndRSh),
     &                      iWork(ip_iSP2F),
     &                      MaxRed,nSym,nnShl,nnBstRT(1),iRed,.false.)
            Call Cho_SetRedInd(iWork(ip_iiBstRSh),iWork(ip_nnBstRSh),
     &                         nSym,nnShl,iLoc)
            irc = 0
            If (iRed .eq. 1) Then ! set correct IndRed array
               kOff2 = kOff2 - 1
               Do iab = 1,nnBstRT(1)
                  kOff = kOff2 + iab
                  iWork(kOff) = iab
               End Do
            End If
         End If
      Else
         irc = 1
      End If

      End
