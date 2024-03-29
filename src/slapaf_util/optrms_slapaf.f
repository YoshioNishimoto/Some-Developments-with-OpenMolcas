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
      Subroutine OptRMS_Slapaf(x,y,nAt,RMS,RMSMax)
      Use Symmetry_Info, only: VarR, VarT
      Implicit Real*8 (a-h,o-z)
#include "print.fh"
#include "real.fh"
      Real*8 x(3,nAt), y(3,nAt)
*
*---- Only align if energy is rotational and translational invariant.
*     (no weighting)
*
      If (.not.(VarR.or.VarT)) Then
         Call Superpose(x,y,nAt,RMS,RMSMax)
      Else
*
*---- Otherwise, just compute RMS
*
         RMS = Zero
         RMSMax = Zero
         Do i = 1, nAt
            disp = Zero
            Do ixyz = 1, 3
               diff = x(ixyz,i)-y(ixyz,i)
               RMS = RMS + diff**2
               disp = disp + diff**2
            End Do
            If (Sqrt(disp).gt.RMSMax) RMSMax=Sqrt(disp)
         End Do
         RMS = Sqrt( RMS/DBLE(nAt) )
      End If
*
      Return
      End
