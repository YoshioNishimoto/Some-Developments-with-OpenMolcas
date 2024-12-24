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
! Copyright (C) 2003, Per-Olof Widmark                                 *
!               2017, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine Optim2(E_Pred,G,H,D,C,n,nDim,n0,n1,r2)
!***********************************************************************
!                                                                      *
!    purpose: Solve a set of non-linear equations to get interpolation *
!             coefficients for damping                                 *
!                                                                      *
!     input:                                                           *
!       G       : linear energy terms                                  *
!       H       : bilinear energy terms                                *
!       D       : bilinear density terms                               *
!       n       : size of matrices                                     *
!       nDim    : leading dimension for H                              *
!       n0      : index number of the lowest energy                    *
!       n1      : index number of the second lowest energy             *
!       r2      : square of density change constraint                  *
!                                                                      *
!     output:                                                          *
!       C       : interpolation coefficients                           *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! The KISS method is employed here. Do a scan of the energy surface    *
! along lines parallel to the edges. Reduce the steplength gradually.  *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark                                                     *
!     University of Lund, Sweden, 2003                                 *
!                                                                      *
!     hacked for optimization on a sphere by:                          *
!     Roland Lindh                                                     *
!     Uppsala University, Uppsala, Sweden, 2017.                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
      Use Constants, only: Zero, Half, One, Two
      Implicit None
!----------------------------------------------------------------------*
! Dummy arguments.                                                     *
!----------------------------------------------------------------------*
      Integer n,nDim,n0,n1
      Real*8  G(nDim),H(nDim,nDim),C(nDim),D(nDim,nDim)
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
      Real*8  Step,Eref,E_Pred
      Real*8  Step_pi,Step_mi,Step_pj,Step_mj,Step_pk,Step_mk
      Real*8  Optim_E
#ifdef _DEBUGPRINT_
      Real*8  sum
#endif
      Logical DidChange
      Integer Iter
      Integer i,j,k,l,nLow
      Real*8  A, B, r2
      Real*8  AI, AJ1, AK1, AJK, AJ2, AK2
      Real*8  BI, BJ1, BJ2
      Real*8  Steps(3,4), Es(4)
!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
!----------------------------------------------------------------------*
! Make first guess.                                                    *
!----------------------------------------------------------------------*
!
!     The initial guessed values of the Cs are set by finding the
!     coeffcients for the densities with the two lowest energies,
!     n0 and n1. All other coefficients are assumed to be zero.
!     The coefficient are found by solving the two simultaneous
!     equations,
!
!     Sum_i C(i) = 0
!
!     Sum_i C(i) C(j) D(i,j) = r2
!
      Do i=1,n
         C(i)=Zero
      End Do
      A=Two*(D(n0,n1)-D(n1,n1))/(D(n0,n0)-Two*D(n0,n1)+D(n1,n1))
      B=(D(n1,n1)-r2)/(D(n0,n0)-Two*D(n0,n1)+D(n1,n1))
      C(n0) = -(A/Two) + Sqrt((A/Two)**2-B)
      If (C(n0).gt.One) C(n0) = -(A/Two) - Sqrt((A/Two)**2-B)
      If (C(n0).lt.Zero) Then
         Write (6,*) 'C(n0).lt.Zero'
         Call Abend()
      End If
      C(n1)=One-C(n0)
      If (C(n1).gt.One) Then
         Write (6,*) 'C(n1).gt.One'
         Call Abend()
      End If
      If (C(n1).lt.Zero) Then
         Write (6,*) 'C(n1).lt.Zero'
         Call Abend()
      End If
!
!     At this point we have a set of Cs which fulfil the two contraints.
#ifdef _DEBUGPRINT_
      Write(6,'(a)') 'Start C:'
      Write(6,'(6F15.6)') (C(i),i=1,n)
      Write(6,'(a)') 'Start G:'
      Write(6,'(6F15.6)') (G(i),i=1,n)
      Call RecPrt('H',' ',H,nDim,nDim)
      Call RecPrt('D',' ',D,nDim,nDim)
#endif
!----------------------------------------------------------------------*
! Compute start energy.                                                *
!----------------------------------------------------------------------*
      Eref=Zero
      Do i=1,n
         Eref=Eref+C(i)*G(i)
         Do j=1,n
            Eref=Eref+Half*C(i)*C(j)*H(i,j)
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Write(6,'(a,F15.6)') 'Eref = ',Eref
#endif
!----------------------------------------------------------------------*
! Do a scan in all tripletwise directions.                             *
!----------------------------------------------------------------------*
      Iter=0
      DidChange=.true.
      Step=0.10d0
      Call Abend()
!
!     Below we explore changes to all triplets that are consistent with
!     that the constrains are not broken.
!
100   If(Iter.lt.500 .and. DidChange) Then
!
         Iter=Iter+1
         DidChange=.false.
!
         Do i=1,n-2
!
            Do j=i+1,n-1
!
               Do k=j+1, n
!
!                 Make sure that step trivially is within the range.
!                 The first constraint. Now we vary d(j) and d(k)
!                 such that the constraints are not broken.
!
                  Step_pi= Min(Step,One-C(i))
                  Steps(1,1) = Step_pi
                  Steps(1,2) = Step_pi
!
!                 Given the change d(i) compute d(k) and d(l),
!                 subjects to the simultaneous equations.
!
!                 Sum_i d(i) = 0
!
!                 Sum_ij (2(C(i)+d(i))*d(j)*D(i,j) = 0
!
                  AI=Two*(C(i)+    Step_pi)*Step_pi*D(i,i)+C(k)*Step_pi*D(k,i)+C(j)*Step_pi*D(j,i)
                  AJ1=   (C(i)+Two*Step_pi)*        D(i,j)+C(j)*        D(j,j)+C(k)*        D(k,j)
                  AK1=   (C(i)+Two*Step_pi)*        D(i,k)+C(k)*        D(k,k)+C(j)*        D(j,k)
                  AJK=Two*D(j,k)
                  AJ2=D(j,j)
                  AK2=D(k,k)
!
                  BI = AI - AJ1*Step_pi+AJ2*Step_pi**2
                  BJ1= -AJ1 + AK1 - (AJK - Two*AJ2)*Step_pi
                  BJ2= -AJK + AJ2 + AK2
!
                  A = BJ1/(Two*BJ2)
                  B = BI/BJ2
!
                  Step_pk = - A + Sqrt(A**2-B)
                  Step_mk = - A - Sqrt(A**2-B)
!
                  Step_pj=-(Step_pi+Step_pk)
                  Step_mj=-(Step_pi+Step_mk)
!
                  C(i)=C(i) + Step_pi
!
                  C(j)=C(j) + Step_mj
                  C(k)=C(k) + Step_mk
                  Steps(2,1)=Step_mj
                  Steps(3,1)=Step_mk
!
                  Es(1)=optim_E(C,G,H,n,nDim)
!
                  C(j)=C(j) - Step_mj
                  C(k)=C(k) - Step_mk
!
                  C(j)=C(j) + Step_pj
                  C(k)=C(k) + Step_pk
                  Steps(2,2)=Step_pj
                  Steps(3,2)=Step_pk
!
                  Es(2)=optim_E(C,G,H,n,nDim)
!
                  C(j)=C(j) - Step_pj
                  C(k)=C(k) - Step_pk
!
                  C(i)=C(i) - Step_pi
!
                  Step_mi= Min(C(i),Step)
                  Steps(1,3) = Step_mi
                  Steps(1,4) = Step_mi
!
                  AI =Two*(C(i)+    Step_mi)*Step_mi*D(i,i)+C(k)*Step_mi*D(k,i)+C(j)*Step_mi*D(j,i)
                  AJ1=    (C(i)+Two*Step_mi)*        D(i,j)+C(j)*        D(j,j)+C(k)        *D(k,j)
                  AK1=    (C(i)+Two*Step_mi)*        D(i,k)+C(k)*        D(k,k)+C(j)        *D(j,k)
                  AJK=Two*D(j,k)
                  AJ2=D(j,j)
                  AK2=D(k,k)
!
                  BI = AI - AJ1*Step_mi+AJ2*Step_mi**2
                  BJ1= -AJ1 + AK1 - (AJK - Two*AJ2)*Step_mi
                  BJ2= -AJK + AJ2 + AK2
!
                  A = BJ1/(Two*BJ2)
                  B = BI/BJ2
!
                  Step_pk = - A + Sqrt(A**2-B)
                  Step_mk = - A - Sqrt(A**2-B)
!
                  Step_pj=-(Step_mi+Step_pk)
                  Step_mj=-(Step_mi+Step_mk)
!
                  C(i)=C(i) + Step_mi
!
                  C(j)=C(j) + Step_mj
                  C(k)=C(k) + Step_mk
                  Steps(2,3)=Step_mj
                  Steps(3,3)=Step_mk
!
                  Es(3)=optim_E(C,G,H,n,nDim)
!
                  C(j)=C(j) - Step_mj
                  C(k)=C(k) - Step_mk
!
                  C(j)=C(j) + Step_pj
                  C(k)=C(k) + Step_pk
                  Steps(2,4)=Step_pj
                  Steps(3,4)=Step_pk
!
                  Es(4)=optim_E(C,G,H,n,nDim)
!
                  C(j)=C(j) - Step_pj
                  C(k)=C(k) - Step_pk
!
                  C(i)=C(i) - Step_mi
!
#ifdef _DEBUGPRINT_
                  Write(6,'(a,3I3,4F20.10)')
     &                  'i,j,k,Es(i),i=1,4)=',
     &                   i,j,k,(Es(l),l=1,4)
#endif
                  nLow=0
                  Do l = 1, 4
                     If (Es(l).lt.Eref) Then
                        nLow=l
                        Eref=Es(l)
                        DidChange=.true.
                     End If
                  End Do
                  If (nLow.ne.0) Then
                     C(i) = C(i) + Steps(1,nLow)
                     C(j) = C(j) + Steps(2,nLow)
                     C(k) = C(k) + Steps(3,nLow)
                  End If
!
               End Do
            End Do
         End Do
!
         If (.not.DidChange) Then
            If (Step.gt.0.9d-4) Then
               Step=0.1d0*Step
               DidChange=.true.
#ifdef _DEBUGPRINT_
               Write(6,*) 'Step is',Step
#endif
            End If
         End If
#ifdef _DEBUGPRINT_
!
!        Check that the constraint is fullfilled.
!
         sum=Zero
         Do i=1,n
            sum=sum+C(i)
         End Do
         Write(6,*) 'optim: sum-1',sum-One
!
         sum=Zero
         Do i = 1, n
            Do j = 1, n
               sum=sum+C(i)*C(j)*D(i,j)
            End Do
         End Do
         Write(6,*) 'optim: sum-r2',sum-r2
#endif
         Go To 100
      End If
#ifdef _DEBUGPRINT_
      Write (6,*) 'ERef=',ERef
#endif
      E_Pred=ERef

      End SubRoutine Optim2
