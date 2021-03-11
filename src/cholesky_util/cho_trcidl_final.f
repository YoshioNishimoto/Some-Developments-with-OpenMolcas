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
      SubRoutine Cho_TrcIdl_Final()
C
C     Thomas Bondo Pedersen, May 2010.
C
C     Deallocate array for tracing idle processors
C
      Implicit None
#include "choptr2.fh"

      If (l_Idle.gt.0) Then
         Call GetMem('TrcIdl','Free','Inte',ip_Idle,l_Idle)
         l_Idle=0
      End If

      End
