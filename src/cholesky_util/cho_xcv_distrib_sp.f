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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Cho_XCV_Distrib_SP(mySP,l_mySP,N_mySP)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Determine distribution of Shell Pairs according to their
C     dimension.
C
      Use Para_Info, Only: MyRank, nProcs
      use ChoSwp, only: nnBstRSh
      Implicit None
      Integer l_mySP
      Integer mySP(l_mySP)
      Integer N_mySP
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "stdalloc.fh"

      Integer  Cho_iFindSmallest
      External Cho_iFindSmallest

      Integer iSP, iNode, n, iSym
      Integer, Allocatable:: ProcDim(:)

#if defined (_DEBUGPRINT_)
      If (l_mySP.lt.nnShl) Then
         Call Cho_Quit('Dimension error in Cho_XCV_Distrib_SP',103)
      End If
#endif

      If (Cho_Real_Par) Then
         Call mma_allocate(ProcDim,nProcs,Label='ProcDim')
         ProcDim(:)=0

         N_mySP=0
         Do iSP=1,nnShl
            n=nnBstRSh(1,iSP,1)
            Do iSym=2,nSym
               n=n+nnBstRSh(iSym,iSP,1)
            End Do
            If (n.gt.0) Then
               iNode=Cho_iFindSmallest(ProcDim,nProcs)
               ProcDim(iNode)=ProcDim(iNode)+n
               If (iNode-1.eq.myRank) Then
                  N_mySP=N_mySP+1
                  mySP(N_mySP)=iSP
               End If
            End If
         End Do

         Call mma_deallocate(ProcDim)
      Else
         N_mySP=0
         Do iSP=1,nnShl
            n=nnBstRSh(1,iSP,1)
            Do iSym=2,nSym
               n=n+nnBstRSh(iSym,iSP,1)
            End Do
            If (n.gt.0) Then
               N_mySP=N_mySP+1
               mySP(N_mySP)=iSP
            End If
         End Do
      End If

      End
