************************************************************************
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
* Copyright (C) 1994,1997, Roland Lindh                                *
*               2014, Ignacio Fdez. Galvan                             *
************************************************************************
      Subroutine RS_P_RFO(H,g,nInter,dq,UpMeth,dqHdq,StepMax,
     &                    Step_Trunc)
************************************************************************
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             December '94                                             *
*                                                                      *
*                                                                      *
*             Solve      | H    g | | d |     | 1 0 | | d |            *
*                        |  T     | |   | = e |     | |   |            *
*                        | g    0 | | 1 |     | 0 1 | | 1 |            *
*                                                                      *
*             this corresponds to                                      *
*                                                                      *
*             H d + g = e d                                            *
*                                                                      *
*             and                                                      *
*                                                                      *
*               T                                                      *
*             g  d = e                                                 *
*                                                                      *
*             Modified from single negative eigenvalue to an arbitrary *
*             number, June '97, R. Lindh                               *
*                                                                      *
*             Removed full diagonalizations, April '14, I. Fdez. Galvan*
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
      Real*8 H(nInter,nInter), g(nInter), dq(nInter), Lambda
      Real*8, Allocatable:: Mat(:), Val(:), Vec(:,:), Tmp(:,:)
      Real*8, Allocatable:: MatN(:), ValN(:), VecN(:), TmpN(:),
     &        StepN(:), GradN(:)
      Real*8, Allocatable:: MatP(:), ValP(:), VecP(:), TmpP(:),
     &        StepP(:), GradP(:)
*
      Character*6 UpMeth
      Character*1 Step_Trunc
      Logical Found,Iterate
*
      iRout = 215
      Lu =6
      iPrint = nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt(' In RS_P_RFO: H','(10f10.6)',H,nInter,nInter)
         Call RecPrt(' In RS_P_RFO: g','(10f10.6)', g,nInter,1)
         Call RecPrt(' In RS_P_RFO:dq','(10f10.6)',dq,nInter,1)
      End If
*
      UpMeth='RSPRFO'
*
      NumVal=Min(2,nInter)
      nVStep=2
      Found=.False.
      Thr=1.0D-6
      Call mma_allocate(Vec,nInter,NumVal,Label='Vec')
      Call mma_allocate(Val,NumVal,Label='Val')
      Call mma_allocate(Mat,nInter*(nInter+1)/2,Label='Mat')
      Vec(:,:)=Zero
      Do i = 1, nInter
         Do j = 1, i
            ij = i*(i-1)/2+j
            Mat(ij)=H(i,j)
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Find the negative eigenvalue(s)
*     Stop when the highest eigenvalue found is larger than Thr
      Do While (.Not.Found)
        Call Davidson(Mat,nInter,NumVal,Val,Vec,iStatus)
        If (iStatus.gt.0) Then
          Call SysWarnMsg('RS_P_RFO',
     &      'Davidson procedure did not converge','')
        End If
        If ((Val(NumVal).gt.Thr).or.(NumVal.ge.nInter)) Then
          Found=.True.
        Else
*----     Increase the number of eigenpairs to compute
          Call mma_allocate(Tmp,nInter,NumVal,Label='Tmp')
          Tmp(:,:) = Vec(:,:)
          Call mma_deallocate(Vec)
          Call mma_deallocate(Val)

          NumVal=Min(NumVal+nVStep,nInter)

          Call mma_allocate(Vec,nInter,NumVal,Label='Vec')
          Vec(:,:)=Zero
          Call mma_allocate(Val,NumVal,Label='Val')
          Val(:)=Zero

          Vec(:,1:NumVal-nVStep) = Tmp(:,:)

          Call mma_deallocate(Tmp)
        End If
      End Do
      Call mma_deallocate(Mat)
*
      nNeg=0
      i=NumVal
      Do While ((i.ge.0).and.(nNeg.eq.0))
         If (Val(i).lt.Zero) nNeg=i
         i=i-1
      End Do
      If (iPrint.ge.99) Then
         Call RecPrt(' In RS_P_RFO: Eigenvalues',' ',Val,1,NumVal)
         Call RecPrt(' In RS_P_RFO: Eigenvectors',' ',Vec,nInter,NumVal)
         Write (Lu,*) ' nNeg=',nNeg
      End If
*
      If (iPrint.ge.6) Then
         Write (Lu,*)
         Write (Lu,*) 'RS-P-RF Optimization'
         Write (Lu,*) ' Iter   alpha        dqdq  StepMax'//
     &                '   EigVal_r  EigVal_t'
      End If
*
*     Write (Lu,*) 'Trust radius=',StepMax
      A_RFO=One   ! Initial seed of alpha
      IterMx=25
      Iter=0
      Iterate=.False.
      Thr=1.0D-7
      If (nNeg.gt.0) Then
         mInter=nNeg+1
         Call mma_allocate(StepN,nInter,Label='StepN')
         Call mma_allocate(GradN,nInter,Label='GradN')
         Call mma_allocate(VecN,mInter,Label='VecN')
         Call mma_allocate(ValN,1,Label='ValN')
         Call mma_allocate(MatN,mInter*(mInter+1)/2,Label='MatN')
         Call mma_allocate(TmpN,mInter,Label='TmpN')
         TmpN(:)=Zero
      End If
      mInter=nInter+1
      Call mma_allocate(StepP,nInter,Label='StepP')
      Call mma_allocate(GradP,nInter,Label='GradP')
      Call mma_allocate(VecP,mInter,Label='VecP')
      Call mma_allocate(ValP,1,Label='ValP')
      Call mma_allocate(MatP,mInter*(mInter+1)/2,Label='MatP')
      Call mma_allocate(TmpP,mInter,Label='TmpP')
      TmpP(:) = Zero
 998  Continue
         Iter=Iter+1
*        Write (Lu,*) 'Iter=',Iter
*        Write (Lu,*) 'A_RFO=',A_RFO
         Call FZero(dq,nInter)
         If (nNeg.gt.0) Then
*           write(Lu,*)
*           write(Lu,*) 'Process negative eigenvalues.'
*           write(Lu,*)
            mInter=nNeg+1
*                                                                      *
************************************************************************
*                                                                      *
*--------   Build the augmented matrix of the "negative" subspace
*            diagonal: negative eigenvalues (divided by alpha)
*            last row/column: components of the gradient along the
*                     negative eigenvectors (divided by sqrt(alpha))
*           Project the gradient in the negative subspace (but expressed
*           in the full space)
*           (Note that, since we are interested in the largest eigenvalue,
*            the whole matrix is multiplied by -1, so we actually find the
*            smallest eigenvalue and change its sign)
            GradN(:)=Zero
            MatN(:) =Zero
            j=mInter*(mInter-1)/2
            Do i=1,nNeg
               MatN(i*(i+1)/2)=-Val(i)/A_RFO
               gv=DDot_(nInter,g,1,Vec(:,i),1)
               MatN(j+i)=gv/Sqrt(A_RFO)
               Call DaXpY_(nInter,gv,Vec(:,i),1,GradN,1)
            End Do
*
*--------   Solve the partial RFO system for the negative subspace
            VecN(:) = TmpN(:)
            Call Davidson(MatN,mInter,1,ValN,VecN,iStatus)
            TmpN(:) = VecN(:)
            If (iStatus.gt.0) Then
               Call SysWarnMsg('RS_P_RFO',
     &              'Davidson procedure did not converge','')
            End If
            ValN(1) = - ValN(1)
*
*--------   Scale the eigenvector (combines eqs. (5) and (23))
*           Convert to full space and add to complete step
            Call DScal_(nNeg,One/(Sqrt(A_RFO)*VecN(1+nNeg)),VecN,1)
            Call dGeMV_('N',nInter,nNeg,One,Vec,nInter,
     &                                  VecN,1,Zero,StepN,1)
            Call DaXpY_(nInter,One,StepN,1,dq,1)
*           dqdq_max=Sqrt(DDot_(nInter,StepN,1,StepN,1))
*           write (Lu,*) 'dqdq_max=',dqdq_max
!           Sign
            EigVal_r=-DDot_(nInter,StepN,1,GradN,1)
            If (iPrint.ge.99) Then
               Call RecPrt('dq_r',' ',StepN,1,nInter)
               Call RecPrt(' g_r',' ',GradN,1,nInter)
               Write (Lu,*) 'Lambda=',EigVal_r
            End If
            If (EigVal_r.lt.-Thr) Then
               Write (Lu,*)
               Write (Lu,*) 'W A R N I N G !'
               Write (Lu,*) 'EigVal_r.lt.Zero',EigVal_r
               Write (Lu,*)
            End If
         Else
            EigVal_r=Zero
*           dqdq_max=Zero
         End If
*
*        write(Lu,*)
*        write(Lu,*) 'Process positive eigenvalues.'
*        write(Lu,*)
         mInter=nInter+1
*                                                                      *
************************************************************************
*                                                                      *
*----    Build the augmented matrix of the "positive" subspace
*        Instead of reducing the dimensions, the negative eigenvectors
*        are simply projected out from the gradient and the eigenvalues
*        are shifted to a large positive arbitrary value (10), to avoid
*        interferences
         GradP(:) = g(:)
         Do i=1,nNeg
           gv=DDot_(nInter,GradP(:),1,Vec(:,i),1)
           Call DaXpY_(nInter,-gv,Vec(:,i),1,GradP(:),1)
         End Do
         Do j=1,nInter
           call dcopy_(j,H(1,j),1,MatP(1+j*(j-1)/2),1)
           Do i=1,nNeg
             Do k=1,j
               jk=j*(j-1)/2+k
               MatP(jk)=MatP(jk)-(Val(i)-Ten)*Vec(j,i)*Vec(k,i)
             End Do
           End Do
           Call DScal_(j,One/A_RFO,MatP(1+j*(j-1)/2),1)
         End Do
         Call FZero(MatP(1+mInter*(mInter-1)/2),mInter)
         Call DaXpY_(nInter,-One/Sqrt(A_RFO),GradP(:),1,
     &                     MatP(1+mInter*(mInter-1)/2),1)
*
*----    Solve the partial RFO system for the positive subspace
         call dcopy_(mInter,TmpP(:),1,VecP(:),1)
         Call Davidson(MatP,mInter,1,ValP,VecP,iStatus)
         If (iStatus.gt.0) Then
           Call SysWarnMsg('RS_P_RFO',
     &          'Davidson procedure did not converge','')
         End If
         TmpP(1:mInter) = VecP(1:mInter)
         StepP(1:nInter) = VecP(1:nInter)
*
*----    Scale the eigenvector (combines eqs. (5) and (23))
*        Add to complete step
         Call DScal_(nInter,One/(Sqrt(A_RFO)*VecP(1+nInter)),StepP,1)
         Call DaXpY_(nInter,One,StepP,1,dq,1)
*        dqdq_min=Sqrt(DDot_(nInter,StepP,1,StepP,1))
*        write (Lu,*) 'dqdq_min=',dqdq_min
         EigVal_t=-DDot_(nInter,StepP,1,GradP,1) ! Sign
         If (iPrint.ge.99) Then
           Call RecPrt('dq_t',' ',StepP,1,nInter)
           Call RecPrt(' g_t',' ',GradP,1,nInter)
           Write (Lu,*) 'Lambda=',EigVal_t
         End If
         If (EigVal_t.gt.Thr) Then
           Write (Lu,*)
           Write (Lu,*) 'W A R N I N G !'
           Write (Lu,*) 'EigVal_t.gt.Zero',EigVal_t
           Write (Lu,*)
         End If
*
      Lambda = EigVal_t + EigVal_r
      dqdq=Sqrt(DDot_(nInter,dq,1,dq,1))
*
      If (iPrint.ge.6)
     &  Write (Lu,'(I5,5F10.5)') Iter,A_RFO,dqdq,StepMax,
     &                           EigVal_r,EigVal_t
*                                                                      *
************************************************************************
*                                                                      *
*------- Initialize data for iterative scheme (only at first iteration)
*
         If (.Not.Iterate) Then
            A_RFO_long=A_RFO
            dqdq_long=dqdq
            A_RFO_short=Zero
            dqdq_short=dqdq_long+One
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- RF with constraints. Start iteration scheme if computed step
*        is too long.
*
         If (Iter.eq.1.and.dqdq.gt.StepMax) Iterate=.True.
*                                                                      *
************************************************************************
*                                                                      *
*        Procedure if the step length is not equal to the trust radius
*
         If (Iterate.and.Abs(StepMax-dqdq).gt.Thr) Then
            Step_Trunc='*'
*           Write (Lu,*) 'StepMax-dqdq=',StepMax-dqdq
            Call Find_RFO_Root(A_RFO_long,dqdq_long,
     &                         A_RFO_short,dqdq_short,
     &                         A_RFO,dqdq,StepMax)
            If (Iter.gt.IterMx) Then
               Write (Lu,*) ' Too many iterations in RF'
               Go To 997
            End If
            Go To 998
         End If
*
 997  Continue

      Call mma_deallocate(Vec)
      Call mma_deallocate(Val)

      If (nNeg.gt.0) Then
         mInter=nNeg+1
         Call mma_deallocate(StepN)
         Call mma_deallocate(GradN)
         Call mma_deallocate(VecN)
         Call mma_deallocate(ValN)
         Call mma_deallocate(MatN)
         Call mma_deallocate(TmpN)
      End If
      mInter=nInter+1
      Call mma_deallocate(StepP)
      Call mma_deallocate(GradP)
      Call mma_deallocate(VecP)
      Call mma_deallocate(ValP)
      Call mma_deallocate(MatP)
      Call mma_deallocate(TmpP)
*
      If (iPrint.ge.6) Then
         Write (Lu,*)
         Write (Lu,*)
         Write (Lu,*) 'Rational Function Optimization: Lambda=',Lambda
         Write (Lu,*)
      End If
      dqHdq=dqHdq+Lambda*Half
*
      If (iPrint.ge.99) Then
         Write (Lu,*) 'EigVal,dqHdq=',Lambda,dqHdq
         Call RecPrt(' In RS_P_RFO: g','(10f10.6)', g,nInter,1)
         Call RecPrt(' In RS_P_RFO:dq','(10f10.6)',dq,nInter,1)
      End If
*
      Return
      End
