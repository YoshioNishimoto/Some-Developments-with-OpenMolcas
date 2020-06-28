!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! (c) Copyright, Roland Lindh  2020                                    *
!***********************************************************************
!                                                                      *
! For a set of d degenerate eigenvectors P with the length n, rotate   *
! the vectors to achieve optimal delacalization.                       *
!                                                                      *
!***********************************************************************
      Subroutine canonical_degenerate_eigenvectors(P_Old,n,d,Option)
      Use kriging_mod, only: blavAI
      Implicit None
      Integer, Intent(In):: n, d
      Real*8, Intent(InOut):: P_Old(n,d)
      Logical, Intent(In)::  Option
!
      External DDot_
      Real*8 DDot_
#include "stdalloc.fh"
      Integer k, kOpt, i, j, ij, nAngles, iAngles, jAngles,      &
              kAngles, iter, iSing, nRaw
      Real*8 D_k, Tmp, grad2, gradMax, Det, D_k_Sum
      Real*8 :: Thr1=1.0D-8, Thr2=1.0D-8, Hii
      Real*8, Allocatable:: Angles(:), QSub(:,:,:), PQ(:,:), Q(:,:)
      Real*8, Allocatable:: dAngles(:), dQSub(:,:,:), dPQ(:,:,:),      &
                            dQ(:,:,:), QTmp(:,:,:), P(:,:)
      Real*8, Allocatable:: Delta_Angles(:)
      Real*8, Parameter:: Zero=0.0D0, One=1.0D0
      Real*8, Allocatable:: HessInv(:,:), Hess(:,:)
      Integer, Parameter:: iter_max=100
      Real*8, Allocatable:: f_of_phi(:), phi(:,:), dfdphi(:,:)
      Real*8 :: Last_Value=0.0D0, Factor=0.0D0
      Logical :: Converged=.True.
!#define _DEBUG_
!#define _DEBUG2_
!#define _ONE_ROW_
!#define _NUMERICAL_TEST_
#ifdef _NUMERICAL_TEST_
      Real*8, Allocatable:: PQnum(:,:,:), Qnum(:,:,:), QSubNum(:,:,:)
      Real*8, Parameter:: Two=2.0D0
      Real*8 Delta
!
#endif
!
      Call mma_allocate(P,n,d,Label='deg:P')
      P(:,:)=P_Old(:,:)
      If (Option) Then
         Factor = -1.0D0
      Else
         Factor =  1.0D0
      End If
!
!     If (d.gt.2) Call Abend()
!
!     Find the row of the vectors P where the sum of the absolute values
!     of the coefficients are as large as possible. We will use this row
!     to compute the degree of the delocalization.
!
      kOpt=0
      D_k=Zero
      D_k_Sum = Zero
      Do k = 1, n
         Tmp=Zero
         Do i = 1, d
            Tmp=Tmp+Abs(P(k,i))
         End Do
         If (Tmp.gt.D_k) Then
            kOpt=k
            D_k=Tmp
         End If
         D_k_Sum = D_k_Sum + Tmp
      End Do
#ifdef _DEBUG_
      Call RecPrt('P(in)',' ',P,n,d)
#ifdef _ONE_ROW_
      Write (6,*)
      Write (6,*) 'kOpt=',kOpt
      Write (6,*) 'D_k=', Factor * D_k
      Write (6,*) 'P(kOpt,i)=',(P(kOpt,i),i=1,d)
      Write (6,*)
#else
      Write (6,*)
      Write (6,*) 'D_k_Sum=',Factor * D_k_Sum
      Write (6,*)
#endif
#endif
!
!     After this we will use the negative of D_k as the target function. In this
!     Respect we will look for a minimum and the Hessian will be positive definite.
!
!     The number of independent angles in the anti-symmetric rotational matrix.
      nAngles=d*(d-1)/2
!
!     Allocate memory for intermediate products
      Call mma_allocate(QTmp,d,d,2,Label='deg:QTmp')
      QTmp(:,:,:)=Zero
!
!     Allocate memory for the product PQ
      Call mma_allocate(PQ,n,d,Label='deg:PQ')
      Call mma_allocate(dPQ,n,d,nAngles,Label='deg:dPQ')
      PQ(:,:)=Zero
      dPQ(:,:,:)=Zero
!
!     Allocate memory for the rotational angles -- one for each pair in P
      Call mma_allocate(Angles,nAngles,Label='deg:Angles')
      Call mma_allocate(dAngles,nAngles,Label='deg:dAngles')
      Call mma_allocate(Delta_Angles,nAngles,Label='deg:Delta_Angles')
      Angles(:)=Zero
      dAngles(:)=Zero
      Delta_Angles(:)=Zero
      Call mma_allocate(f_of_phi,iter_max,Label='deg:f_of_phi')
      Call mma_allocate(phi,nAngles,iter_max,Label='deg:phi')
      Call mma_allocate(dfdphi,nAngles,iter_max,Label='deg:dfdphi')
!
!     The Q matrix is expressed as the product of the submatrices - QSub(1) x QSub(2) x ... x QSub(nAngles)
      Call mma_allocate(Q,d,d,Label='deg:Q')
      Call mma_allocate(dQ,d,d,nAngles,Label='deg:dQ')
      Q(:,:)=Zero
      Forall (i=1:d) Q(i,i)=One
!
!     Allocate memory for the individual submatrices based on the rotational angles
      Call mma_allocate(QSub,d,d,nAngles,Label='deg:QSub')
      Call mma_allocate(dQSub,d,d,nAngles,Label='deg:dQSub')  ! derivative
      QSub(:,:,:)=Zero
      dQSub(:,:,:)=Zero
!
!     Allocate memory for the appoximate inverse Hessian
      Call mma_allocate(HessInv,nAngles,nAngles,Label='deg:HessInv')
      Call mma_allocate(Hess,nAngles,nAngles,Label='deg:Hess')
      Hess(:,:)=Zero
      ForAll (i=1:nAngles) Hess(i,i)=1.0D0
      HessInv(:,:)=Hess(:,:)
!
      Call Init_Kriging()
#ifdef _ONE_ROW_
      blavAI= Abs(D_k)  !   set the base line.
#else
      blavAI= Abs(D_k_Sum)  !   set the base line.
#endif
!
!***************************************************************************************************
!***************************************************************************************************
!
!     Start of loop
!
      iter=0

      Do     ! inifinte loop
!
!***************************************************************************************************
!***************************************************************************************************
!
         iter=iter+1
         If (iter.gt.iter_max) Then
            Converged=.False.
            Exit
         End if
!
!     Given the values of all angles compute all the submatrices.
!
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) '***********************************************************************************'
      Write (6,*) '******* Iteration: ',iter
      Write (6,*) '***********************************************************************************'
#endif
#ifdef _NUMERICAL_TEST_
      Write (6,*)
      Write (6,*) '********************** Start numerical section ************************************'
      Write (6,*)
      Call mma_allocate(PQNum,n,d,2,Label='PQNum')
      Call mma_allocate(QNum,d,d,2,Label='QNum')
      Call mma_allocate(QSubNum,d,d,2,Label='QSubNum')
      dAngles(:)=Angles(:) ! temporary save
      Do kAngles = 1, nAngles
!
      Write (6,*)
      Write (6,*)
      Write (6,*) ' ***********************************************************'
      Write (6,*) ' ***************** kAngles: ',kAngles
      Write (6,*) ' ***********************************************************'
      Write (6,*)
      Write (6,*) ' **** x + d **** '
      Delta=1.0D-5m
      Angles(kAngles)=dAngles(kAngles)+Delta
!     Call RecPrt('Angles',' ',Angles,1,nAngles)
      Call Mk_Submatrices()
      Call Mk_Q()
      QNum(:,:,1)=Q(:,:)
      QSubNum(:,:,1)=QSub(:,:,kAngles)
      Call Mk_PQ()
      PQNum(:,:,1)=PQ(:,:)
      Call Compute_D_k
!     Write (6,*) 'D_k=',D_k
      tmp=D_k
!
      Write (6,*)
      Write (6,*) ' **** x - d **** '
      Angles(kAngles)=dAngles(kAngles)-Delta
      Call RecPrt('Angles',' ',Angles,1,nAngles)
      Call Mk_Submatrices()
      Call Mk_Q()
      QNum(:,:,2)=Q(:,:)
      QSubNum(:,:,2)=QSub(:,:,kAngles)
      Call Mk_PQ()
      PQNum(:,:,2)=PQ(:,:)
      Call Compute_D_k
!     Write (6,*) 'D_k=',D_k
!
      Write (6,*)
      Write (6,*) '***** numerical result ****'
      tmp = (tmp - D_k)/(Two*Delta)
      Write (6,*) 'Numerical gradient:',tmp
      Q(:,:)=(QNum(:,:,1)-QNum(:,:,2))/(Two*Delta)
      Call RecPrt('dQ(numerical)',' ',Q,d,d)
      Q(:,:)=(QSubNum(:,:,1)-QSubNum(:,:,2))/(Two*Delta)
      Call RecPrt('dQSub(numerical)',' ',Q,d,d)
      PQ(:,:)=(PQNum(:,:,1)-PQNum(:,:,2))/(Two*Delta)
#ifdef _ONE_ROW_
      Write (6,*) 'dPQ(kOpt,i)(numerical)=',(PQ(kOpt,i),i=1,d)
#endif
      Angles(:)=dAngles(:) ! Restore
      dfdphi(kAngles,iter)=tmp
      End Do
!
      Call mma_deallocate(QSubNum)
      Call mma_deallocate(QNum)
      Call mma_deallocate(PQNum)
!
      Write (6,*)
      Write (6,*) '************************ End numerical section ************************************'
      Write (6,*)
#endif
!
      phi(:,iter)=Angles(:)
      Call Mk_Submatrices()
      Call Mk_Q()
      Call Mk_PQ()
      Call Compute_D_k
      f_of_phi(iter)=D_k
#ifdef _DEBUG_
      Call RecPrt('Angles',' ',Angles,1,nAngles)
      Call RecPrt('Phi(:,iter)',' ',Phi(1,iter),1,nAngles)
      Write (6,*)
      Write (6,*) 'iter, D_k=',iter, D_k
      Write (6,*)
#endif
!
      Call Mk_derivatives
      dfdphi(:,iter)=dAngles(:)
#ifdef _DEBUG3_
      Call RecPrt('Phi(:,iter)',' ',Phi(1,iter),1,nAngles)
      Call RecPrt('dfdphi',' ',dfdphi(1,iter),1,nAngles)
#endif
!
!*******************************************************************************************
!
!     Check for convergence, if not do a simple Newton-Raphson update
!
!     Write (6,*) 'grad2=',grad2
!     Write (6,*) 'gradMax=',gradMax
!     Write (6,*) Sqrt(grad2/DBLE(nAngles)),Thr2
!     Write (6,*) gradMax,Thr1
      grad2=grad2/D_k**2
      gradMax=gradMax/Abs(D_k)


      If ((Sqrt(grad2/DBLE(nAngles)).lt.Thr2  .and. gradMax.lt.Thr1) .or.     &
          Abs(Last_value-D_k).le.Thr1)  Exit
      Last_value=D_k
!
!*******************************************************************************************
!
!     Use gradient-enhanced kriging to get the Hessian. Note that we will not micro iterate.
!
!     Always start from a positive definite matrix, if not the characteristic length will
!     be screwed.
      Hess(:,:)=Zero
      ForAll (i=1:nAngles) Hess(i,i)=1.0D0
#ifdef _DEBUG3_
      Call RecPrt('Hess(Unit)',' ',Hess,nAngles,nAngles)
      Call RecPrt('f_of_phi',' ',f_of_phi,1,iter)
      Call RecPrt('Phi(:,iter)',' ',Phi(1,iter),1,nAngles)
      Call RecPrt('dfdphi',' ',dfdphi(1,iter),1,nAngles)
#endif
#define _NEW_
#ifdef _NEW_
      nRaw=Min(2*d+1,iter)
      Call SetUp_Kriging(nRaw,nAngles,Phi(1,iter-nRaw+1),        &
                                   dfdphi(1,iter-nRaw+1),        &
                                   f_of_phi(iter-nRaw+1),Hess)
#else
      Call SetUp_Kriging(iter,nAngles,Phi,dfdphi,f_of_phi,Hess)
#endif
      Call Hessian_Kriging_Layer(phi(1,iter),Hess,nAngles)
      Call Close_Kriging()
#ifdef _DEBUG3_
      Call RecPrt('Hess(Kriging)',' ',Hess,nAngles,nAngles)
#endif
!     Make sure that the Hessian is positive definite -- we are going for a minimum!
      Call pdHess(Hess,nAngles)
#ifdef _DEBUG3_
      Call RecPrt('Hess(pdHess)',' ',Hess,nAngles,nAngles)
#endif
      Call Minv(Hess,HessInv,iSing,Det,nAngles)
#ifdef _DEBUG3_
      Call RecPrt('HessInv',' ',HessInv,nAngles,nAngles)
#endif
!
!*******************************************************************************************
!
!     For now, do a trivial Newton-Raphson update of the angles to find a maximum
#ifdef _DEBUG_
      Call RecPrt('dAngles',' ',dAngles,nAngles,1)
#endif
      Call DGEMM_('n','n',nAngles,1,nAngles,   &
                  One, HessInv,nAngles,       &
                       dAngles,nAngles,       &
                  Zero,Delta_Angles,nAngles)
#ifdef _DEBUG_
      Call RecPrt('Delta_Angles',' ',Delta_Angles,nAngles,1)
#endif
      Tmp=Sqrt(DDot_(nAngles,Delta_Angles,1,Delta_Angles,1))
!     Scale step if too large
      If (Tmp.gt.0.50D0/DBLE(iter)) Delta_Angles(:)=(0.50D0/(Tmp*DBLE(iter)))*Delta_Angles(:)
#ifdef _DEBUG_
      Call RecPrt('Delta_Angles',' ',Delta_Angles,nAngles,1)
#endif
!#define _INC_
#ifdef _INC_
      Angles(:)= - Delta_Angles(:)
      P(:,:)=PQ(:,:)
#else
      Angles(:)=Angles(:) - Delta_Angles(:)
#endif
!
!***************************************************************************************************
!***************************************************************************************************
!
      End Do
!
!***************************************************************************************************
!***************************************************************************************************
!
!     Here if converged. If not, skip procedure all together.
!
#ifdef _DEBUG_
      Call RecPrt('PQ(final)',' ',PQ,n,d)
#endif
      If (Converged) P_Old(:,:)=PQ(:,:)
!
!***************************************************************************************************
!     Release the resources
!
      Call mma_deallocate(Hess)
      Call mma_deallocate(HessInv)
      Call mma_deallocate(dQSub)
      Call mma_deallocate(QSub)
      Call mma_deallocate(dQ)
      Call mma_deallocate(Q)
      Call mma_deallocate(dfdphi)
      Call mma_deallocate(phi)
      Call mma_deallocate(f_of_phi)
      Call mma_deallocate(Delta_Angles)
      Call mma_deallocate(dAngles)
      Call mma_deallocate(Angles)
      Call mma_deallocate(dPQ)
      Call mma_deallocate(PQ)
      Call mma_deallocate(QTmp)
      Call mma_deallocate(P)
!
      Return
!
!***************************************************************************************************
!***************************************************************************************************
!
      Contains
!
!***************************************************************************************************
!***************************************************************************************************
!
      Subroutine Mk_Submatrices()
!
!***************************************************************************************************
!
      iAngles=0
      Do i = 1 , d
         Do j = i+1, d
            iAngles=iAngles+1
!           Write (6,*) 'iAngles=',iAngles
!
            ForAll (k=1:d) QSub(k,k,iAngles)=One
!           Write (6,*) 'Cos(phi),Sin(phi)=',Cos(Angles(iAngles)),Sin(Angles(iAngles))
            QSub(i,i,iAngles)= Cos(Angles(iAngles))
            QSub(i,j,iAngles)=-Sin(Angles(iAngles))
            QSub(j,i,iAngles)= Sin(Angles(iAngles))
            QSub(j,j,iAngles)= Cos(Angles(iAngles))
!
            dQSub(i,i,iAngles)=-Sin(Angles(iAngles))
            dQSub(i,j,iAngles)=-Cos(Angles(iAngles))
            dQSub(j,i,iAngles)= Cos(Angles(iAngles))
            dQSub(j,j,iAngles)=-Sin(Angles(iAngles))
#ifdef _DEBUG2_
            Call RecPrt('QSub_i(analytic)',' ',  QSub(1,1,iAngles),d,d)
            Call RecPrt('dQSub_i(analytic)',' ',dQSub(1,1,iAngles),d,d)
#endif
         End Do
      End Do
      Return
      End Subroutine Mk_Submatrices
!
!***************************************************************************************************
!***************************************************************************************************
!
      Subroutine Mk_Q()
!
!***************************************************************************************************
!
!     Compute Q =  Q_1*Q_2*...*Q_N
!
      i=1
      j=2
      QTmp(:,:,i)=QSub(:,:,1)
      Do iAngles = 2, nAngles
!
         Call DGEMM_('N','N',d,d,d,            &
                     One, QTmp(1,1,i),d,       &
                          QSub(1,1,iAngles),d, &
                     Zero,QTmp(1,1,j),d)
!
         i=MOD(i,2)+1
         j=MOD(j,2)+1
      End Do
      Q(:,:)=QTmp(:,:,i)
#ifdef _DEBUG2_
      Call RecPrt('Q',' ',Q,d,d)
#endif
      Return
      End Subroutine Mk_Q
!
!***************************************************************************************************
!***************************************************************************************************
!
      Subroutine Mk_PQ()
!
!***************************************************************************************************
!
!     Compute P*Q
!
      Call DGEMM_('N','N',n,d,d,    &
                  One ,P,n,         &
                       Q,d,         &
                  Zero,PQ,n)
!     Call RecPrt('PQ',' ',PQ,n,d)
#ifdef _DEBUG2_
#ifdef _ONE_ROW_
      Write (6,*) 'PQ(kOpt,i)=',(PQ(kOpt,i),i=1,d)
#endif
#endif
      Return
      End Subroutine Mk_PQ
!
!***************************************************************************************************
!***************************************************************************************************
!
      Subroutine Compute_D_k()
!
!***************************************************************************************************
!
!     Compute the value of the function to optimize
!
      D_k=Zero
#ifdef _ONE_ROW_
      Do i = 1, d
         D_k=D_k + Factor * Abs(PQ(kOpt,i))
      End Do
#else
      Do i = 1, d
         Do k = 1, n
             D_k=D_k + Factor * Abs(PQ(k,i))
         End Do
      End Do
#endif
      Return
      End Subroutine Compute_D_k
!
!***************************************************************************************************
!***************************************************************************************************
!
      Subroutine Mk_derivatives()
!
!***************************************************************************************************
!
!     Write (6,*)
!     Write (6,*) '************* Enter Mk_derivatives******************'
!     Write (6,*)
!
!     Let us now do the derivatives
!
      grad2=Zero
      gradMax=Zero
!
!     Loop over all the rotational angles.
!
      Do iAngles = 1, nAngles
!
!        Form: Q_1*Q*_2*...*dQ_iAngles*...*Q_nAngles
!        Write (6,*)
!        Write (6,*) ' *************** iAngles= ',iAngles
!        Write (6,*)
         i=1
         j=2
         If (iAngles.eq.1) Then
            QTmp(:,:,i)=dQSub(:,:,1)
         Else
            QTmp(:,:,i)= QSub(:,:,1)
         End If
         Do jAngles = 2, nAngles
!
!           Call RecPrt('A',' ',QTmp(1,1,i),d,d)
            If (jAngles.eq.iAngles) Then
!              Call RecPrt('B',' ',dQSub(1,1,jAngles),d,d)
               Call DGEMM_('N','N',d,d,d,            &
                           One, QTmp(1,1,i),d,       &
                               dQSub(1,1,jAngles),d, &
                           Zero,QTmp(1,1,j),d)
            Else
!              Call RecPrt('B',' ', QSub(1,1,jAngles),d,d)
               Call DGEMM_('N','N',d,d,d,            &
                           One, QTmp(1,1,i),d,       &
                                QSub(1,1,jAngles),d, &
                           Zero,QTmp(1,1,j),d)
            End If
!           Call RecPrt('C',' ',QTmp(1,1,j),d,d)
!
!           Effectively interchange the values of i and j
            i=MOD(i,2)+1
            j=MOD(j,2)+1
         End Do  ! jAngles
!        Write (*,*) 'i=',i
!
!        The result of the last multiplication is in the matrix indexed i!
!
         dQ(:,:,iAngles)=QTmp(:,:,i)
#ifdef _DEBUG2_
         Call RecPrt('dQ',' ',dQ(1,1,iAngles),d,d)
#endif
!
!        Compute P*dQ
!
         Call DGEMM_('N','N',n,d,d,                  &
                     One ,P,n,                       &
                          dQ(1,1,iAngles),d,         &
                     Zero,dPQ(1,1,iAngles),n)
#ifdef _DEBUG2_
         Call RecPrt('dPQ(Analytic)',' ',dPQ(1,1,iAngles),n,d)
#ifdef _ONE_ROW_
         Write (6,*) 'dPQ(kOpt,i)(Analytic)=',(dPQ(kOpt,i,iAngles),i=1,d)
#endif
#endif
!
!        Form the derivative
!
         Tmp = Zero
#ifdef _ONE_ROW_
         Do i = 1, d
            Tmp = Tmp + Factor * SIGN(One,PQ(kOpt,i)) * dPQ(kOpt,i,iAngles)
         End Do
#else
         Do i = 1, d
            Do k = 1, n
               Tmp = Tmp + Factor * SIGN(One,PQ(k,i)) * dPQ(k,i,iAngles)
            End Do
         End Do
#endif
         dAngles(iAngles)=Tmp
!
         grad2=grad2+dAngles(iAngles)**2
         gradMax=Max(gradMax,Abs(dAngles(iAngles)))
      End Do
#ifdef _DEBUG2_
      Call RecPrt('dAngles',' ',dAngles,1,nAngles)
#endif
!     Write (6,*)
!     Write (6,*) '************* Exit Mk_derivatives ******************'
!     Write (6,*)
      Return
      End Subroutine Mk_derivatives
!
!***************************************************************************************************
!
      End Subroutine canonical_degenerate_eigenvectors
