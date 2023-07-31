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
      SubRoutine Integral_WrOut_Cho_diag(
#define _FIXED_FORMAT_
#define _CALLING_
#include "int_wrout_interface.fh"
     &                                  )
*     calls the proper routines IndSft/PLF
*     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
*                      else  integral order within symblk: ijkl
      Implicit Real*8 (A-H,O-Z)
*
#include "int_wrout_interface.fh"
*
* call sorting routine
*
      If (mSym.eq.1) Then
        Call PLF_Cho_Diag(TInt,nTInt,
     &           AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &           iShell,iAO,iAOst,Shijij.and.IJeqKL,
     &           iBas,jBas,kBas,lBas,kOp)
      Else
        Call IndSft_Cho_Diag(TInt,nTInt,
     &               iCmp,iShell,iBas,jBas,kBas,lBas,Shijij,
     &               iAO,iAOst,ijkl,SOInt,nSOint,iSOSym,nSOs)
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
