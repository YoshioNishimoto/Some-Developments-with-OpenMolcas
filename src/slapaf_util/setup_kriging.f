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
* Copyright (C) 2020, Roland Lindh                                     *
************************************************************************
      Subroutine SetUp_Kriging(nRaw,iFirst)
      Use kriging_mod, only: blavAI, set_l, layer_U, nSet
      use Slapaf_Info, only: qInt, dqInt, Energy, dqInt_Aux,
     &                       Energy0
      Implicit None
      Integer nRaw, iFirst

#include "stdalloc.fh"
#include "real.fh"
      Integer i,iInter,jInter,ij, iS, iE, nInter
      Real*8 Value_l
      Real*8, Allocatable:: Array_l(:), HTri(:), Hessian(:,:),
     &                      qInt_s(:,:), dqInt_s(:,:,:),
     &                      Hessian_HMF(:,:), Energy_s(:,:)
*                                                                      *
************************************************************************
*                                                                      *
      nInter=Size(qInt,1)
      iS = iFirst
      iE = iFirst + nRaw - 1
*                                                                      *
************************************************************************
*                                                                      *
*     Pick up the HMF Hessian, used to set the characteristic lengths
*
      Call mma_Allocate(Hessian_HMF,nInter,nInter,Label='Hessian_HMF')
      Call Mk_Hss_Q()
      Call Get_dArray('Hss_Q',Hessian_HMF,nInter**2)
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) 'setup_Kriging, nRaw=',nRaw
      Call RecPrt('Setup_kriging: Energy',' ',Energy(iS:iE),1,nRaw)
      Call RecPrt('Setup_kriging: qInt',' ',qInt(:,iS:iE),nInter,nRaw)
      Call RecPrt('Setup_kriging: dqInt',' ',dqInt(:,iS:iE),nInter,nRaw)
      Call RecPrt('HMF Hessian',' ',Hessian_HMF,nInter,nInter)
      If (nSet>=2) Then
         Call RecPrt('Setup_kriging: Energy0',' ',Energy0(iS:iE),1,nRaw)
         Call RecPrt('Setup_kriging: dqInt_Aux(:,:,1)',' ',
     &               dqInt_Aux(:,iS:iE,1),nInter,nRaw)
      End If
      If (nSet==3) Call RecPrt('Setup_kriging: dqInt_Aux(:,:,2)',' ',
     &               dqInt_Aux(:,iS:iE,2),nInter,nRaw)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(layer_U,nInter,nInter,Label='layer_U')
*
      layer_U(:,:)=Zero
      Do i=1,nInter
         layer_U(i,i)=One
      End Do
*
      Call mma_allocate(Hessian,nInter,nInter,Label='Hessian')
      Call mma_allocate(HTri,nInter*(nInter+1)/2,Label='HTri')
*
      Do iInter = 1, nInter
         Do jInter = 1, iInter
            ij = iInter*(iInter-1)/2 + jInter
            HTri(ij)=Hessian_HMF(iInter,jInter)
         End Do
      End Do

      Call mma_deallocate(Hessian_HMF)

*     Call TriPrt('HTri(raw)',' ',HTri,nInter)
      Call NIDiag_new(HTri,layer_U,nInter,nInter)
*     layer_U(:,:)=Zero
*     Do i=1,nInter
*        layer_U(i,i)=One
*     End Do
*     Call TriPrt('HTri',' ',HTri,nInter)
      Hessian(:,:) = Zero
      Do i=1,nInter
         Hessian(i,i)=HTri(i*(i+1)/2)
      End Do
*
      Call mma_deallocate(HTri)
*                                                                      *
************************************************************************
*                                                                      *
*     Select between setting all ls to a single value or go in
*     multiple l-value mode in which the l-value is set such that
*     the kriging hessian reproduce the diagonal value of the HMF
*     Hessian of the current structure.
*
      Call mma_Allocate(Array_l,nInter,Label='Array_l')
      If (Set_l) Then
         Call Get_dScalar('Value_l',Value_l)
         Array_l(:)=Value_l
      Else
         Call Set_l_Array(Array_l,nInter,blavAI,Hessian)
      End If
      Call mma_deallocate(Hessian)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_Allocate(qInt_s,nInter,nRaw,Label="qInt_s")
      Call mma_Allocate(dqInt_s,nInter,nRaw,nSet,Label="dqInt_s")
      Call mma_Allocate(Energy_s,nRaw,nSet,Label="Energy_s")
*                                                                      *
************************************************************************
*                                                                      *
*     Transform to the basis which diagonalizes the HMF Hessian.
*
      Call Trans_K(qInt(:,iS:iE),qInt_s,nInter,nRaw)

      Call Trans_K(dqInt(:,iS:iE),dqInt_s(:,:,1),nInter,nRaw)
      dqInt_s(:,:,1) = - dqInt_s(:,:,1)

      Energy_s(1:nRaw,1)=Energy(iS:iE)

      If (nSet>1) Then
         If (.NOT.Allocated(dqInt_Aux)) Call Abend()
         Call Trans_K(dqInt_Aux(:,iS:iE,1),dqInt_s(:,:,2),nInter,nRaw)
         dqInt_s(:,:,2) = - dqInt_s(:,:,2)
         Energy_s(1:nRaw,2)=Energy0(iS:iE)
         If (nSet==2) Then
            Energy_s(1:nRaw,2) = Energy_s(1:nRaw,1) - Energy_s(1:nRaw,2)
            dqInt_s(:,:,2)     = dqInt_s(:,:,1)     - dqInt_s(:,:,2)
         End If
      End If

      If (nSet>2) Then
         If (SIZE(dqInt_Aux,3)/=2) Call Abend()
         Call Trans_K(dqInt_Aux(:,iS:iE,2),dqInt_s(:,:,3),nInter,nRaw)
         dqInt_s(:,:,3) = - dqInt_s(:,:,3)

         Energy_s(1:nRaw,3)=Zero
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call RecPrt('Setup_kriging: qInt_s',' ',qInt_s,nInter,nRaw)
      Do i = 1, nSet
         Write (6,*) 'iSet=',i
         Call RecPrt('Setup_kriging: Energy_s',' ',Energy_s(:,i),1,nRaw)
         Call RecPrt('Setup_kriging: Grad_s',' ',dqInt_s(:,:,i),nInter,
     &                                                             nRaw)
      End Do
#endif
      Call Start_Kriging(nRaw,nInter,nSet,qInt_s,dqInt_s,Energy_s)
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deAllocate(Energy_s)
      Call mma_deAllocate(dqInt_s)
      Call mma_deAllocate(qInt_s)
*                                                                      *
************************************************************************
*                                                                      *
*     Pass the l-values to the GEK routine. This will initiate the
*     computation of the covariance matrix, and solve related GEK
*     equations.
*
      Call Set_l_Kriging(Array_l,nInter)
      Call mma_deAllocate(Array_l)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine Setup_Kriging
