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
      SubRoutine Int_LDF_JK_2P(
#define _FIXED_FORMAT_
#define _CALLING_
#include "int_wrout_interface.fh"
     &                        )
*     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
*                      else  integral order within symblk: ijkl
      Implicit Real*8 (a-h,o-z)
*
#include "localdf_int2.fh"
*
#include "int_wrout_interface.fh"
*
      External LDF_nShell, LDF_nAuxShell
*
* call sorting routine
*
      If (mSym==1) Then
         nS_Val=LDF_nShell()
         nS_Aux=LDF_nAuxShell()
         iS_Dum=nS_Val+nS_Aux+1
         If (SHA.eq.iS_Dum .and.
     &       SHB.gt.nS_Val .and. SHB.lt.iS_Dum .and.
     &       SHC.eq.iS_Dum .and.
     &       SHD.gt.nS_Val .and. SHD.lt.iS_Dum) Then
            ! type (J|L)
            Call PLF_LDF_JK_2P_1(TInt,nTInt,MapOrg,
     &                       AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                       iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
         Else If (SHA.eq.iS_Dum .and.
     &            SHB.gt.nS_Val .and. SHB.lt.iS_Dum .and.
     &            SHC.le.nS_Val .and.
     &            SHD.le.nS_Val) Then
            ! type (J|uv)
            Call PLF_LDF_JK_2P_2(TInt,nTInt,MapOrg,
     &                       AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                       iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
         Else If (SHA.le.nS_Val .and.
     &            SHB.le.nS_Val .and.
     &            SHC.eq.iS_Dum .and.
     &            SHD.gt.nS_Val .and. SHD.lt.iS_Dum) Then
            ! type (uv|J)
            Call PLF_LDF_JK_2P_3(TInt,nTInt,MapOrg,
     &                       AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                       iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
         Else If (SHA.le.nS_Val .and.
     &            SHB.le.nS_Val .and.
     &            SHC.le.nS_Val .and.
     &            SHD.le.nS_Val) Then
            ! type (uv|kl)
            Call PLF_LDF_JK_2P_4(TInt,nTInt,MapOrg,
     &                       AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                       iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
         Else
            Call WarningMessage(2,
     &             'Shell combination not implemented in Int_LDF_JK_2P')
            Write(6,'(A,4I9)')
     &      'SHA,SHB,SHC,SHD........',SHA,SHB,SHC,SHD
            Write(6,'(A,3I9)')
     &      'nS_Val,nS_Aux,iS_Dum...',nS_Val,nS_Aux,iS_Dum
            Call LDF_Quit(1)
         End If
      Else
         Call WarningMessage(2,
     &                      'Symmetry not implemented in Int_LDF_JK_2P')
         Call LDF_Quit(1)
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer_array(iShell)
         Call Unused_logical(Shijij)
         Call Unused_logical(IJeqKL)
         Call Unused_real_array(SOInt)
         Call Unused_integer(nSOint)
         Call Unused_integer_array(iSOSym)
         Call Unused_integer(nSkal)
         Call Unused_integer_array(itOffs)
      End If
      End
