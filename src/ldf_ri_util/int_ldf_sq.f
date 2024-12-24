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
      SubRoutine Int_LDF_SQ(
#define _FIXED_FORMAT_
#define _CALLING_
#include "int_wrout_interface.fh"
     &                     )
*     calls the proper routines IndSft/PLF
*     if IntOrd_jikl==.TRUE. integral order within symblk: jikl
*                      else  integral order within symblk: ijkl
      Implicit Real*8 (a-h,o-z)
*
#include "localdf_int.fh"
*
#include "int_wrout_interface.fh"
*
      External LDF_nShell, LDF_nAuxShell
*
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
*
* check that shells have not been re-ordered
*
      If (MapOrg(1).ne.1 .or. MapOrg(2).ne.2 .or.
     &    MapOrg(3).ne.3 .or. MapOrg(4).ne.4) Then
         Call WarningMessage(2,
     &   'Int_LDF_SQ: Shell reordering not implemented!')
         Write(6,'(A,4I9)')
     &   'MapOrg.................',(MapOrg(i),i=1,4)
         Write(6,'(A,4I9)')
     &   'SHA,SHB,SHC,SHD........',SHA,SHB,SHC,SHD
         Write(6,'(A,2(9X,I9))')
     &   'SHAB,SHCD..............',iTri(SHA,SHB),iTri(SHC,SHD)
         nS_Val=LDF_nShell()
         Write(6,'(A,I9)')
     &   'nS_Val.................',nS_Val
         Call LDF_Quit(1)
      End If
*
* call sorting routine
*
      If (mSym==1) Then
         nS_Val=LDF_nShell()
         If (SHA.le.nS_Val .and.
     &       SHB.le.nS_Val .and.
     &       SHC.le.nS_Val .and.
     &       SHD.le.nS_Val) Then
            ! type (ij|kl)
            Call PLF_LDF_SQ(TInt,nTInt,
     &                      AOInt,ijkl,iCmp(1),iCmp(2),iCmp(3),iCmp(4),
     &                      iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
         Else
            Call WarningMessage(2,
     &      'Shell combination not implemented in Int_LDF_SQ')
            Write(6,'(A,4I9)')
     &      'SHA,SHB,SHC,SHD........',SHA,SHB,SHC,SHD
            Write(6,'(A,I9)')
     &      'nS_Val.................',nS_Val
            Call LDF_Quit(1)
         End If
      Else
         Call WarningMessage(2,
     &               'Symmetry not implemented in Int_LDF_SQ')
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
