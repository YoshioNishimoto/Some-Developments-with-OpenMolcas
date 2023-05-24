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
      SubRoutine Cho_GnVc_GetInt(xInt,lInt,nVecRS,iVecRS,ListSp,
     &                           mSym,mPass,mmShl,iPass1,NumPass,NumSP)
!
!     Purpose: compute integrals for NumPass integral passes starting at
!              pass iPass1.
!
      use ChoSwp, only: IndRSh, InfVec
      Implicit Real*8 (a-h,o-z)
      Real*8  xInt(lInt)
      Integer nVecRS(mSym,mPass), iVecRS(mSym,mPass), ListSP(mmShl)
#include "cholesky.fh"
#include "stdalloc.fh"

      Character(LEN=15), Parameter:: SecNam = 'Cho_GnVc_GetInt'

      Integer, Allocatable:: SPTmp(:)

      Integer, External:: Cho_F2SP

!     Initialization and input check.
!     -------------------------------

      If (NumPass .lt. 1) Then
         NumSP = 0
         return
      End If

      If (mSym .ne. nSym) Then
         Call Cho_Quit('Input error [1] in '//SecNam,103)
      End If

      If (iPass1 .lt. 1) Then
         Call Cho_Quit('Input error [2] in '//SecNam,103)
      End If

      iPass2 = iPass1 + NumPass - 1
      If (iPass2 .gt. mPass) Then
         Call Cho_Quit('Input error [3] in '//SecNam,103)
      End If

      If (mmShl .lt. nnShl) Then
         Call Cho_Quit('Input error [4] in '//SecNam,103)
      End If

!     Set up list of shell pairs to compute.
!     --------------------------------------

      Call mma_allocate(SPTmp,nnShl,Label='SPTmp')
      SPTmp(:)=0
      NumSP = 0
      Do iPass = iPass1,iPass2
         Do iSym = 1,nSym
            iV1 = iVecRS(iSym,iPass)
            iV2 = iV1 + nVecRS(iSym,iPass) - 1
            Do iV = iV1,iV2
               iAB   = InfVec(iV,1,iSym) ! addr in 1st reduced set
               jShAB = IndRsh(iAB) ! shell pair (full)
               iShAB = Cho_F2SP(jShAB) ! reduced shell pair
               If (iShAB .gt. 0) Then
                  If (SPTmp(iShAB) .eq. 0) Then ! register SP
                     SPTmp(iShAB) = 1
                     NumSP = NumSP + 1
                     ListSP(NumSP) = iShAB
                  End If
               Else
                  Call Cho_Quit('SP not found in reduced list!',103)
               End If
            End Do
         End Do
      End Do
      Call mma_deallocate(SPTmp)

!     Set memory used by Seward.
!     --------------------------

      Call mma_maxDBLE(lSewInt)
      Call xSetMem_Ints(lSewInt)

!     Loop through shell pair list.
!     -----------------------------

      Do iSP = 1,NumSP

!        Get shell pair index.
!        ---------------------

         iShAB = ListSP(iSP)

!        Compute integral distribution (**|A B).
!        ---------------------------------------

         Call Cho_MCA_CalcInt_3(xInt,lInt,iShAB)

      End Do

!     Deallocation.
!     -------------

      Call xRlsMem_Ints
      End
