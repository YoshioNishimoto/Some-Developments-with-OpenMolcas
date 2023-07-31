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
* Copyright (C) 1990, IBM                                              *
*               1991, Roland Lindh                                     *
************************************************************************
      Logical Function TstFnc(iCoSet,iIrrep,iBsFnc,nStab)
************************************************************************
*                                                                      *
* Object: to establish if a function is a basis function of a          *
*         irreducible representation.                                  *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             September 1991                                           *
************************************************************************
      use Symmetry_Info, only: nIrrep, iOper, iChTbl
      Implicit Real*8 (A-H,O-Z)
      Integer iCoSet(0:7,0:7), iAcc(0:7)
      Integer iBsFnc
*
      TstFnc = .True.
      nCoSet=nIrrep/nStab
      iAcc(0:nCoSet-1)=0
*
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) 'TstFnc'
      Write (6,*)
      Write (6,*) 'Coset:'
      Do i = 0, nCoSet-1
         Write (6,'(8I4)') (iCoSet(i,j),j=0,nStab-1)
      End Do

      Write (6,*)
      Write (6,*) 'iOper:'
      Write (6,'(8I4)') (iOper(i),i=0,nIrrep-1)
      Write (6,*)
      Write (6,*) 'iBsFnc=',iBsFnc
      Write (6,*)
      Write (6,*) 'iChTbl:'
      Write (6,'(8I4)') (iChTbl(iIrrep,i),i=0,nIrrep-1)
#endif
*
*     Loop over operators
*
      Do i = 0, nIrrep-1
*
*        Find index of the generated center
*
         n = -1
         Do j = 0, nCoSet-1
            If (n.ge.0) Cycle
            Do k = 0, nStab-1
               If (iOper(i).eq.iCoSet(j,k)) n = j
            End Do
         End Do
*
         If (n.lt.0 .or. n.gt.nCoSet-1) Then
            Call WarningMessage(2,'TstFnc: n.lt.0 .or. n.gt.nCoSet-1')
            Write (6,*) ' Coset index',n,' is wrong!'
            Call Abend()
         End If
*
         iCom=iAnd(iOper(i),iBsFnc)
         iAcc(n) = iAcc(n) + iChTbl(iIrrep,i)*iPrmt_(iCom)
*
      End Do
      Do i = 0, nCoSet-1
         If (iAcc(i).eq.0) TstFnc = .False.
      End Do
*
      Return
      End Function TstFnc
*
      Logical Function TF(mdc,iIrrep,iComp)
      Use Center_Info, Only : dc
      Implicit Real*8 (a-h,o-z)
      Logical, External :: TstFnc
      TF = TstFnc(dc(mdc)%iCoSet,iIrrep,iComp,dc(mdc)%nStab)
      End Function TF
*
      Integer Function iPrmt_(iCom)
************************************************************************
*     Returns the phase factor of a basis function under a symmetry    *
*     operation, jOper. iChct contains the information about the       *
*     character of the basis function.                                 *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer i, iCom
      iPrmt_= 1
      Do i = 1, 3
         If (iAnd(iCom,2**(i-1)).ne.0) iPrmt_= iPrmt_*(-1)
      End Do
      Return
      End Function iPrmt_
