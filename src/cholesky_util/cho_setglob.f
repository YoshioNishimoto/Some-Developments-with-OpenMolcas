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
      SubRoutine Cho_SetGlob()
C
C     Purpose: define entries in choglob.fh
C
      Implicit None
#include "choglob.fh"

      Integer N, iSym

      Integer, Parameter:: iLarge=999999

      nnShl_G = 0
      mmBstRT_G = 0

      N = 8*nLoc_G
      Call iZero(iiBstR_G,N)
      Call iZero(nnBstR_G,N)
      Call iZero(nnBstRT_G,nLoc_G)
      Call iZero(NumCho_G,8)
      NumChT_G = 0

      Do iSym = 1,8
         LuCho_G(iSym) = -iLarge
      End Do
      LuRed_G = -iLarge
      LuRst_G = -iLarge

      End
