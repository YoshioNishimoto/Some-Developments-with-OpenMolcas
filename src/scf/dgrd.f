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
*
*     Compute the difference with the previous gradient
*
      Subroutine dGrd()
      use LnkLst, only: SCF_V, LLGrad, LLdGrd
      use InfSCF
      Implicit None
#include "real.fh"
#include "stdalloc.fh"
#include "file.fh"
      Integer jpgrd,inode, i
      Real*8, Dimension(:), Allocatable:: Scr
      Integer, External :: LstPtr
      Call mma_allocate(Scr,mOV,Label='Scr')

#define _NEW_CODE_
#ifdef _NEW_CODE_
*     Do i = 2, iter
      Do i = iter_ref+1, iter
         jpgrd=LstPtr(i,LLGrad)
         Call GetNod(i-1,LLGrad,inode)
         If (inode.eq.0) Then
            Write (6,*) 'inode.eq.0'
            Call Abend()
         End If
        Call iVPtr(Scr,mOV,inode)
        Call DaXpY_(mOV,-One,SCF_V(jpgrd)%A,1,Scr,1)
        Call DScal_(mOV,-One,Scr,1)
        Call PutVec(Scr,mOV,i-1,'OVWR',LLdGrd)
      End Do
#else
      jpgrd=LstPtr(iter,LLGrad)
      Call GetNod(iter-1,LLGrad,inode)
      If (inode.eq.0) Then
         Write (6,*) 'inode.eq.0'
         Call Abend()
      End If
      Call iVPtr(Scr,mOV,inode)
      Call DaXpY_(mOV,-One,SCF_V(jpgrd)%A,1,Scr,1)
      Call DScal_(mOV,-One,Scr,1)
      Call PutVec(Scr,mOV,iter-1,'OVWR',LLdGrd)
#endif
      Call mma_deallocate(Scr)
      Return
      End Subroutine dGrd
