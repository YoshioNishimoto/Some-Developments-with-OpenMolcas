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
      SubRoutine Cho_TrcIdl_Update(IAmIdle)
C
C     Thomas Bondo Pedersen, May 2010.
C
C     Update array for tracing idle processors
C
      Use Para_Info, Only: MyRank
      use ChoArr, only: Idle
      Implicit None
      Logical IAmIdle
#include "cho_para_info.fh"
#if defined (_DEBUGPRINT_)
#include "cholesky.fh"
#endif

#if defined (_DEBUGPRINT_)
      If (.NOT.Allocated(Idle) .or. .not.Trace_Idle) Then
         Write(LuPri,'(A)')
     &   'Cho_TrcIdl_Update should not be called in this run!'
         Write(LuPri,*) 'Trace_Idle=',Trace_Idle
         Call Cho_Quit('Illegal call to Cho_TrcIdl_Update',103)
      End If
#endif

      If (IAmIdle) Then
         If (Cho_Real_Par) Then
            Idle(1+myRank) = Idle(1+myRank)+1
         Else
            Idle(1)=Idle(1)+1
         End If
      End If

      End
