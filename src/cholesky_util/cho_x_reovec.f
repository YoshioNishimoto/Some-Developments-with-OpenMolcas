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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  Cho_X_ReoVec
*
*> @brief
*>   Reorder Cholesky vectors on disk
*> @author Thomas Bondo Pedersen
*>
*> @details
*> @note
*> You should consider using on-the-fly read/reorder routine ::Cho_X_getVFull instead!!
*>
*> This routine makes it possible to reorder Cholesky vectors
*> from reduced storage to full storage in exactly the same
*> manner as when specifying the keywork ``REORder`` in the
*> Cholesky input section (in Seward). If the vectors have
*> already been reordered (checked through runfile), the routine
*> returns immediately.
*> The reordered vectors
*> are stored in the files ``CHFVnm`` where ``n`` is the symmetry of
*> the first AO index, ``m`` that of the second. The resulting files
*> are split (if needed).
*>
*> @note
*> The Cholesky procedure must have been successfully initialized (by ::Cho_X_Init).
*>
*> @param[out] irc return code
************************************************************************
      SubRoutine Cho_X_ReoVec(irc)
      Implicit None
      Integer irc
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Integer ip_Temp, l_Temp, ip_Wrk, l_Wrk, iReo

      irc = 0

      Call Get_iScalar('Cholesky Reorder',iReo)
      If (iReo .eq. 0) Then
         l_Temp = 3*nnBstRT(1)
         Call GetMem('Temp','Allo','Inte',ip_Temp,l_Temp)
         Call GetMem('Maxi','Max ','Real',ip_Wrk,l_Wrk)
         Call GetMem('Work','Allo','Real',ip_Wrk,l_Wrk)
         Call Cho_ReoVec(iWork(ip_Temp),3,nnBstRT(1),Work(ip_Wrk),l_Wrk)
         Call GetMem('Work','Free','Real',ip_Wrk,l_Wrk)
         Call GetMem('Temp','Free','Inte',ip_Temp,l_Temp)
         iReo = 1
         Call Put_iScalar('Cholesky Reorder',iReo)
      End If

      End
