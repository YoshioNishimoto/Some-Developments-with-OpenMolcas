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
      SubRoutine Integral_WrOut_Cho(
#define _FIXED_FORMAT_
#define _CALLING_
#include "int_wrout_interface.fh"
     &                             )
*     calls the proper routines IndSft/PLF
*     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
*                      else  integral order within symblk: ijkl
      Implicit Real*8 (A-H,O-Z)
*
#include "cholesky.fh"

      Character*18 SecNam
      Parameter (SecNam = 'Integral_WrOut_Cho')
*
#include "int_wrout_interface.fh"
*
* call sorting routine
*
      If (IfcSew .eq. 1) Then
         If (nSym.eq.1) Then
           Call PLF_Cho(TInt,nTInt,
     &              AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &              iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &              iBas,jBas,kBas,lBas,kOp)
         Else
           Call IndSft_Cho(TInt,nTInt,
     &                  iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,
     &                  iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
         End If
      Else If (IfcSew .eq. 2) Then
         If (nSym.eq.1) Then
           Call PLF_Cho_2(TInt,nTInt,
     &              AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &              iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &              iBas,jBas,kBas,lBas,kOp)
         Else
           Call IndSft_Cho_2(TInt,nTInt,
     &                  iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,
     &                  iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
         End If
      Else If (IfcSew .eq. 3) Then
         If (nSym.eq.1) Then
           Call PLF_Cho_3(TInt,nTInt,
     &              AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &              iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &              iBas,jBas,kBas,lBas,kOp)
         Else
           Call IndSft_Cho_3(TInt,nTInt,
     &                  iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,
     &                  iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
         End If
      Else
         Write(6,*)
         Write(6,*)
         Write(6,*) '!!!!!!!!!! IfcSew=',IfcSew,' !!!!!!!!!!'
         Call Cho_Quit('IfcSew out of bounds in '//SecNam,105)
      End If
*
      Return
* Avoid unused argument warnings
      IF (.False.) Then
         Call Unused_integer_array(MapOrg)
         Call Unused_integer(nSkal)
         Call Unused_integer_array(itOffs)
      End If
      End
