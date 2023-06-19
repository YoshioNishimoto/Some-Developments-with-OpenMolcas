!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SubRoutine Cho_VecBuf_Copy(Vec,nVec,iSym)
!
!     Purpose: copy as many vectors as fit directly into buffer using
!              current reduced set dimension for symmetry iSym. nVec
!              is the number of vectors in array Vec.
!         NB!  It is important that the vector counter NumCho does NOT
!              include the nVec vectors in array Vec.
!
      use ChoVecBuf, only: CHVBUF, ip_CHVBUF_SYM, l_CHVBUF_SYM,
     &                     nVec_in_Buf
      Implicit None
      Real*8  Vec(*)
      Integer nVec, iSym
#include "cholesky.fh"

      Integer mUsed, Left, nCopy, lCopy, kOff

!     Check if there is anything to do at all.
!     1) any buffer allocated?
!     2) any vectors to copy?
!     3) any elements?
!     ----------------------------------------

      If (l_ChVBuf_Sym(iSym) .lt. 1) Return
      If (nVec .lt. 1) Return
      If (nnBstR(iSym,2) .lt. 1) Return

!     Copy vectors
!     1) if all previous vectors are in the buffer, and
!     2) if there is sufficient free space in the buffer.
!     ---------------------------------------------------

      If (nVec_in_Buf(iSym) .eq. NumCho(iSym)) Then
         mUsed = nnBstR(iSym,2)*nVec_in_Buf(iSym)
         Left  = l_ChVBuf_Sym(iSym) - mUsed
         nCopy = min(Left/nnBstR(iSym,2),nVec)
         If (nCopy .gt. 0) Then
            lCopy = nnBstR(iSym,2)*nCopy
            kOff  = ip_ChVBuf_Sym(iSym) + mUsed
            Call dCopy_(lCopy,Vec,1,CHVBUF(kOff),1)
            nVec_in_Buf(iSym) = nVec_in_Buf(iSym) + nCopy
         End If
      End If

      End
