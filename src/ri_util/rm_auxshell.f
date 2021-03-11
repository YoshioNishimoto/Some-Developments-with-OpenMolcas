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
      Subroutine rm_AuxShell(Info,nInfo,iCnttp)
************************************************************************
*                                                                      *
*     Remove an auxiliary basis set by making it empty.                *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
#include "real.fh"
#include "WrkSpc.fh"
*                                                                      *
************************************************************************
*                                                                      *
      Do k = 0, nTot_Shells(iCnttp)-1
         iShll = ipVal(iCnttp) + k
*
         nExp(iShll) = 0
         nBasis(iShll) = 0
         nBasis_Cntrct(iShll) = 0
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(Info)
         Call Unused_integer(nInfo)
      End If
      End
