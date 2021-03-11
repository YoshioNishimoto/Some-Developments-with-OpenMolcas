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
      Subroutine Get_D2AV(ipD2AV,nDens)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character*24 Label
      Logical      Found

      Label='D2av'
      Call qpg_dArray(Label,Found,nDens)
      If(.not.Found .or. nDens.eq.0) Then
         Call SysAbendMsg('get_d2av','Did not find',Label)
      End If
      Call GetMem('D2av','Allo','Real',ipD2AV,nDens)
      Call get_dArray(Label,Work(ipD2AV),nDens)

      Return
      End
