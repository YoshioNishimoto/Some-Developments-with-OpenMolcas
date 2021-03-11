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
* Copyright (C) 1996, Anders Bernhardsson                              *
************************************************************************
      SubRoutine TrGrd_Alaska(GradIn,nGrad)
********************************************************************
*                                                                  *
*      Transforms a symmetry adapted gradient to unsymmetric  form *
*                                                                  *
*       Written by Anders Bernhardsson                             *
*       960427                                                     *
********************************************************************
      Implicit Real*8(a-h,o-z)
      parameter (tol=1d-8)
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      Real*8 CGrad(3,MxAtom)
      dimension GradIn(nGrad),A(3)
      Logical TF,TstFnc
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                       nIrrep/nStab(mdc),iChTbl,iIrrep,iComp,
     &                       nStab(mdc))
      mdc=0
      iIrrep=0
*
      call dcopy_(3*MxAtom,[Zero],0,CGrad,1)
      iCen=0
      nCnttp_Valence=0
      Do iCnttp = 1, nCnttp
         If (AuxCnttp(iCnttp)) Go To 999
         nCnttp_Valence = nCnttp_Valence+1
      End Do
 999  Continue
*
      Do iCnttp=1,nCnttp_Valence
         ixyz = ipCntr(iCnttp)
         Do iCnt=1,nCntr(iCnttp)
            mdc=mdc+1
            call dcopy_(3,Work(ixyz),1,A,1)
            Do iCo=0,nIrrep/nStab(mdc)-1
               kop=iCoSet(iCo,0,mdc)
               nDispS = IndDsp(mdc,iIrrep)
               A1=DBLE(iPrmt(NrOpr(kop,iOper,nIrrep),1))*A(1)
               A2=DBLE(iPrmt(NrOpr(kop,iOper,nIrrep),2))*A(2)
               A3=DBLE(iPrmt(NrOpr(kop,iOper,nIrrep),4))*A(3)
               iCen=iCen+1
               Do iCar=0,2
                  iComp = 2**iCar
                  If ( TF(mdc,iIrrep,iComp)) Then
                     nDispS = nDispS + 1
                     XR=DBLE(iPrmt(NrOpr(kop,iOper,nIrrep),icomp))
                     CGrad(iCar+1,iCen)=XR*GradIn(nDispS)
                  End If
               End Do
            End Do
         End Do
         ixyz=ixyz+3
      End Do
*
*     Call RecPrt('CGrad',' ',CGrad,3,iCen)
      Return
      End
