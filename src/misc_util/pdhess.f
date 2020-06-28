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
* (c) Copyright, Roland Lindh 2020                                     *
************************************************************************
      Subroutine pdHess(H,nH)
      Implicit None
#include "stdalloc.fh"
      Integer nH
      Real*8 H(nH,nH)
      Integer i, j, k
      Real*8, Allocatable:: EVal(:), EVec(:,:)
      Real*8  Hii
*
*#define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('pdHess: H(Start)',' ',H,nH,nH)
#endif
*
      Call mma_allocate(EVal,nH*(nH+1)/2,Label='Eval')
      Call mma_allocate(EVec,nH,nH,Label='Evec')
*
*---- Copy elements for H
*
      Do i = 1, nH
         Do j = 1, i
            EVal(i*(i-1)/2+j)=H(i,j)
         End Do
      End Do
*
      EVec(:,:)=0.0D0
      ForAll (i=1:nH) EVec(i,i)=1.0D0
#ifdef _DEBUG_
      Call TriPrt('pdHess: EVal(init)',' ',EVal,nH)
      Call RecPrt('pdHess: EVec(init)',' ',EVec,nH,nH)
#endif
*
*---- Compute eigenvalues and eigenvectors
*
      Call NIdiag(EVal,EVec,nH,nH,0)
#ifdef _DEBUG_
      Call TriPrt('pdHess: EVal',' ',EVal,nH)
      Call RecPrt('pdHess: EVec(init)',' ',EVec,nH,nH)
#endif
*
*---- Apply corrections if any ...
*
      Do i = 1, nH
         Hii=EVal(i*(i+1)/2)
         If (Hii.gt.0.0D0) Cycle
         Write (6,*) 'i, Hii=',i,Hii
         Hii=2.0D0*Abs(Hii)
         Do j=1,nH
            Do k=1,nH
            H(j,k)=H(j,k) + Hii*EVec(j,i)*EVec(k,i)
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(EVal)
      Call mma_deallocate(EVec)
#ifdef _DEBUG_
      Call RecPrt('pdHess: H(Final)',' ',H,nH,nH)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
