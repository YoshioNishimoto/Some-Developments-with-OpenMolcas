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
      SubRoutine Cho_RstD_ChkSP2F(iSP2F,l_iSP2F,nErr)
      Implicit None
      Integer l_iSP2F, nErr
      Integer iSP2F(l_iSP2F)
#include "stdalloc.fh"

      Integer i

      Integer, Allocatable:: iChk(:)

      Call mma_allocate(iChk,l_iSP2F,Label='iChk')

      Call Cho_RstD_GetInd3(iChk,l_iSP2F)

      nErr = 0
      Do i = 1, l_iSP2F
         If (iChk(i) .ne. iSP2F(i)) nErr = nErr + 1
      End Do

      Call mma_deallocate(iChk)

      End
