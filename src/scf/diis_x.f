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
* Copyright (C) 1994,2017, Roland Lindh                                *
*               1995, Per-Olof Widmark                                 *
*               1995, Markus P. Fuelscher                              *
*               1995, Piotr Borowski                                   *
*               1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine DIIS_x(nD,CInter,nCI,QNRStp,Ind)
************************************************************************
*                                                                      *
*     purpose: Accelerate convergence using DIIS method                *
*                                                                      *
*     input:                                                           *
*                                                                      *
*     output:                                                          *
*       CInter  : Interpolation coefficients of length nCI             *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*               uses SubRoutines and Functions from Module lnklst.f    *
*               -linked list implementation to store series of vectors *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     modified by:                                                     *
*     P.O. Widmark, M.P. Fuelscher, P. Borowski & M.Schuetz            *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*     Derived from code for c1- and c2-DIIS as implemented by          *
*     R. Lindh in Slapaf and SCF in 1994.                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
*#define _DEBUGPRINT_
*#define _NEW_CODE_
*#define _NEW_
      use InfSO, only: IterSO, Energy
      use InfSCF, only: TimFld, mOV, kOptim, Iter, C1DIIS, AccCon,
     &                  Iter_Start
      use Constants, only: One, Two, Zero
#ifdef _NEW_
      use Constants, only: Half
#endif
#ifdef _NEW_CODE_
      Use InfSCF, only: Iter_Start
#endif
#ifdef _DEBUGPRINT_
      use InfSCF, only: kOV
#endif
      use MxDM, only: MxOptm
      Implicit None
*
#include "stdalloc.fh"
#include "file.fh"
*
#ifdef _NEW_CODE_
      Integer k
      Real*8 E_tmp
#endif
      Integer nCI, nD
      Real*8 CInter(nCI,nD)
      Real*8, Dimension(:,:), Allocatable:: EVector, Bij
      Real*8, Dimension(:), Allocatable:: EValue, Err1, Err2, Scratch
*     Real*8, Dimension(:), Allocatable:: Err3, Err4
*
*---- Define local variables
      Integer Ind(MxOptm)
      Real*8 GDiis(MxOptm + 1),BijTri(MxOptm*(MxOptm + 1)/2)
      Real*8 EMax, Fact, ee2, ee1, E_Min_G, Dummy, Alpha, B11
      Logical QNRstp, Case2
      Integer iVec, jVec, kVec, nBij, nFound
      Integer :: i, j
*     Integer :: iPos
      Integer :: ipBst, ij, iErr, iDiag, iDum
      Real*8 :: tim1, tim2, tim3
      Real*8 :: thrld=1.0D-15
      Real*8 :: ThrCff=4.0D2
      Real*8 :: cpu1, cpu2, c2, Bii_Min
      Real*8, External:: DDot_
#ifdef _NEW_CODE_
      Logical Ignore
#endif
      Character*80 Text,Fmt
#ifdef _DEBUGPRINT_
      Real*8 cDotV
#endif
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
 100  Continue
*
#ifdef _NEW_CODE_
!
!     Select from iterations with the lowest energies.
!
      Do i = kOptim, 1, -1
         Ind(i)=0
*
         E_Min_G= Zero
         Do j = Iter_Start, iter
*
            Ignore=.False.
            Do k = kOptim, i+1, -1
               Ignore = Ignore .or. Ind(k).eq.j
            End Do
            If (Ignore) Cycle
*
            E_tmp = Energy(j)
*
            If (E_tmp.lt.E_Mi_Gn) Then
               E_Min_G=E_tmp
               Ind(i)=j
            End If
*
         End Do
      End Do
#else
!
!     Select from the kOptim last iterations
!
      Do i = 1, kOptim
         Ind(i) = iter-kOptim+i
      End Do
#endif
#ifdef _DEBUGPRINT_
*     Write (6,*) 'Iter, Iter_Start=', Iter, Iter_Start
      Write (6,*) 'kOptim=',kOptim
      Write (6,*) 'Ind(i):',(Ind(i),i=1,kOptim)
#endif
*
*-----The following piece of code computes the DIIS coeffs
*     with error vectors chosen as the grds (DIIS only)
*     or as delta=-Hinv*grd (QNR/DIIS)
*
*     Allocate memory for error vectors (gradient or delta)
*
      Call mma_allocate(Err1,mOV,Label='Err1')
      Call mma_allocate(Err2,mOV,Label='Err2')
*     Call mma_allocate(Err3,mOV,Label='Err3')
*     Call mma_allocate(Err4,mOV,Label='Err4')
      nBij=kOptim+1
      Call mma_allocate(Bij,nBij,nBij)
      Call FZero(Bij,nBij**2)
*
*---- Compute norms, <e_i|e_j>
*
#ifdef _DEBUGPRINT_
      Write (6,*) 'kOV(:)=',kOV
      Write (6,*) 'mOV   =',mOV
      Call RecPrt('Energy',' ',Energy,1,iter)
#endif
      E_Min_G=Zero
      Bii_min=1.0D+99
      Do i=1,kOptim
         Call ErrV(mOV,Ind(i),QNRStp,Err1)
*        If (QNRStp) Call ErrV(mOV,Ind(i),.False.,Err3)
#ifdef _DEBUGPRINT_
         Call NrmClc(Err1,mOV,'Diis  ','Err(i) ')
#endif
         Do j=1,i-1

            Call ErrV(mOV,Ind(j),QNRStp,Err2)
*           If (QNRStp) Call ErrV(mOV,Ind(j),.False.,Err4)
#ifdef _DEBUGPRINT_
            Call NrmClc(Err2,mOV,'Diis  ','Err(j)  ')
#endif
            If (QNRStp) Then
#ifdef _NEW_
*              Bij(i,j) = Half*DBLE(nD)*(
*    &                    DDot_(mOV,Err1,1,Err4,1)
*    &                  + DDot_(mOV,Err3,1,Err2,1)
*    &                                  )
#else
               Bij(i,j) = DBLE(nD)*DDot_(mOV,Err1,1,Err2,1)
#endif
            Else
               Bij(i,j) = DBLE(nD)*DDot_(mOV,Err1,1,Err2,1)
            End If
            Bij(j,i) = Bij(i,j)
         End Do
         If (QNRStp) Then
#ifdef _NEW_
*           Bij(i,i) = DBLE(nD)*DDot_(mOV,Err1,1,Err3,1)
#else
            Bij(i,i) = DBLE(nD)*DDot_(mOV,Err1,1,Err1,1)
#endif
         Else
            Bij(i,i) = DBLE(nD)*DDot_(mOV,Err1,1,Err1,1)
         End If
         If (Bij(i,i).lt.Bii_Min) Then
            E_min_G=Energy(Ind(i))
            Bii_Min=Bij(i,i)
         End If
      End Do

      i = kOptim
!
!     Monitor if the sequence of norms of the error vectors and their
!     corresponding energies are consistent with a single concave
!     potential energy minimum.

!     Case 2
!     Check if we are sliding off a shoulder, that is, we have a
!     lowering of the energy while the norm of the error vector
!     increase.
      Case2 = Bij(i,i)>Bii_Min .and. Energy(Ind(i))+1.0D-4<E_Min_G
     &        .and. kOptim>1

      If ( qNRStp .and. Case2 ) Then
#ifdef _DEBUGPRINT_
         Write(6,*)'   RESETTING kOptim!!!!'
         Write(6,*)'   Calculation of the norms in Diis :'
         Fmt  = '(6f16.8)'
         Text = 'B-matrix squared in Diis :'
         Call RecPrt(Text,Fmt,Bij,nBij,nBij)
         Write (6,'(A,2F16.6)') 'Bij(i,i),      Bii_Min=',Bij(i,i),
     &                                                    Bii_Min
         Write (6,'(A,2F16.6)') 'Energy(Ind(i)),E_Min_G=',
     &                Energy(Ind(i)),E_Min_G
         Write (6,*) Energy(Ind(i)),E_Min_g
         Write (6,*)

#endif
!        Rest the depth of the DIIS and the BFGS update.
         Write(6,*) 'DIIS_X: Resetting kOptim!'
         Write(6,*) '        Caused by energies and gradients which are'
     &            //' inconsistent with a concave energy functional.'
         kOptim=1
         Iter_Start = Iter
         IterSO=1
*        Call mma_deallocate(Err4)
*        Call mma_deallocate(Err3)
         Call mma_deallocate(Err2)
         Call mma_deallocate(Err1)
         Call mma_deallocate(Bij)
         Go To 100
      End If
*
*---- Deallocate memory for error vectors & gradient
*     Call mma_deallocate(Err4)
*     Call mma_deallocate(Err3)
      Call mma_deallocate(Err2)
      Call mma_deallocate(Err1)
*
#ifdef _DEBUGPRINT_
      Write(6,*)'   Calculation of the norms in Diis :'
      Fmt  = '(6f16.8)'
      Text = 'B-matrix squared in Diis :'
      Call RecPrt(Text,Fmt,Bij,nBij,nBij)
      Write(6,*)
      Write(6,*)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Here, the DIIS coeffs for DIIS only or QNR/DIIS are
*     computed, either with C1DIIS or C2DIIS
*     (-> stored in vector CInter)
*                                                                      *
************************************************************************
*                                                                      *
      If (.not.c1Diis) Then
*                                                                      *
*-------   C2DIIS case                                                 *
*                                                                      *
*         References:                                                  *
*         H. Sellers, Int. J. Quantum Chem. 45, 31-41(1993).           *
*         doi:10.1002/qua.560450106                                    *
*                                                                      *
************************************************************************
*                                                                      *
         If (QNRStp) Then
            AccCon = 'QNRc2DIIS'
         Else
            AccCon = 'c2DIIS   '
         End If
*
*------- Form a unit eigenvector matrix
         Call mma_allocate(EVector,kOptim,kOptim,Label='EVector')
         Call mma_allocate(EValue,kOptim,Label='EValue')
*
         call dcopy_(kOptim**2,[Zero],0,EVector,       1)
         call dcopy_(kOptim,   [One], 0,EVector,kOptim+1)
*
*------- Form a triangular B-matrix
*
         ij = 1
         Do i = 1, kOptim
            call dcopy_(i,Bij(i,1),nBij,BijTri(ij),1)
            ij = ij + i
         End Do
*
#ifdef _DEBUGPRINT_
         Fmt  = '(5g25.15)'
         Text = 'B-matrix before Jacobi :'
         Call TriPrt(Text,Fmt,BijTri,kOptim)
         Text = 'EigenVectors before Jacobi :'
         Call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
         Write(6,*)
         Write(6,*)
#endif
*
*------- Diagonalize B-matrix
*
         EMax=Zero
         Do i = 1, kOptim*(kOptim+1)/2
            EMax=Max(EMax,Abs(BijTri(i)))
         End Do
         Do i = 1, kOptim*(kOptim+1)/2
            If (Abs(BijTri(i)).lt.EMax*1.0D-14) BijTri(i)=Zero
         End Do
*
         Call mma_allocate(Scratch,kOptim**2,Label='Scratch')
*
         Dummy=Zero
         iDum=0
         Call Diag_Driver('V','A','L',kOptim,BijTri,
     &                    Scratch,kOptim,Dummy,Dummy,iDum,iDum,
     &                    EValue,EVector,kOptim,1,0,'J',
     &                    nFound,iErr)
*
         Call mma_deallocate(Scratch)
         Call dCopy_(kOptim*(kOptim+1)/2,[Zero],0,BijTri,1)
*
         iDiag = 0
         Do i = 1,kOptim
            iDiag = iDiag + i
            BijTri(iDiag) = EValue(i)
         End Do

*
#ifdef _DEBUGPRINT_
         Fmt  = '(5g25.15)'
         Text = 'B-matrix after Jacobi :'
         Call TriPrt(Text,Fmt,BijTri,kOptim)
         Text = 'EigenValues :'
         Call RecPrt(Text,Fmt,EValue,1,kOptim)
         Text = 'EigenVectors :'
         Call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
         Write(6,*)
         Write(6,*)
#endif
*
*------  Renormalize the eigenvectors to the C1-DIIS format
*
#ifdef _DEBUGPRINT_
         Write(6,*)' Normalization constants :'
#endif
*
         Do iVec = 1, kOptim
            Alpha = Zero
            Do i = 1, kOptim
               Alpha = Alpha + EVector(i,iVec)
            End Do
*
#ifdef _DEBUGPRINT_
            Fmt = '(A7,i2,A4,f16.8)'
            Write(6,Fmt)' Alpha(',iVec,') = ',Alpha
#endif
*
            EVector(:,iVec) = EVector(:,iVec)/Alpha
         End Do

         Do kVec = 1, kOptim
            ee1 = Zero
            Do iVec = 1, kOptim
               Do jVec = 1, kOptim
                  ee1 = ee1 +
     &            EVector(iVec,kVec)*EVector(jVec,kVec)*Bij(iVec,jVec)
               End Do
            End Do
*           Write (6,*)'EValue(kVec),ee1:',EValue(kVec),ee1
            EValue(kVec)=ee1
         End Do
*
#ifdef _DEBUGPRINT_
         Fmt  = '(6e16.8)'
         Text = 'B-matrix after scaling :'
         Call TriPrt(Text,Fmt,BijTri,kOptim)
         Text = 'EigenValues after scaling :'
         Call RecPrt(Text,Fmt,EValue,1,kOptim)
         Text = 'EigenVectors after scaling :'
         Call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
         Write(6,*)
         Write(6,*)
#endif
*
*------  Select a vector.
         ee1   = 1.0D+72
#ifdef _DEBUGPRINT_
         cDotV = 1.0D+72
#endif
         ipBst =-99999999
         Do iVec = 1, kOptim
*
*           Pick up eigenvalue (ee2) and the norm of the
*           eigenvector (c2).
*
            ee2 = EValue(iVec)
            c2 = DDot_(kOptim,EVector(1,iVec),1,EVector(1,iVec),1)
#ifdef _DEBUGPRINT_
            Write (6,*) '<e|e>=',ee2
            Write (6,*) 'c**2=',c2
#endif
*
*---------  Reject if <e|e> is too low (round-off),
*           analys further.
*
            If (ee2.lt.Thrld) Then
#ifdef _DEBUGPRINT_
               Fmt  = '(A,i2,5x,g12.6)'
               Text = '<e|e> is low,         iVec, <e|e> = '
               Write(6,Fmt)Text(1:36),iVec,ee2
#endif
*
*------------  Reject if coefficients are too large (linear dep.).
*
               If (Sqrt(c2).gt.ThrCff) Then
#ifdef _DEBUGPRINT_
                  Fmt  = '(A,i2,5x,g12.6)'
                  Text = 'c**2 is too large,     iVec, c**2 = '
                  Write(6,Fmt)Text(1:36),iVec,c2
#endif
*                 Go To 520
                  Cycle
               End If
            End If
*
*---------  Reject if coefficients are too large (linear dep.).
*
            If (Sqrt(c2).gt.ThrCff*Two) Then
#ifdef _DEBUGPRINT_
               Fmt  = '(A,i2,5x,g12.6)'
               Text = 'c**2 is too large,     iVec, c**2 = '
               Write(6,Fmt)Text(1:36),iVec,c2
#endif
*              Go To 520
               Cycle
            End If
*
*-----------Keep the best candidate
*
*
*---------  Reject if coefficients are too large (linear dep.).
*
            If (Sqrt(c2).gt.ThrCff) Then
#ifdef _DEBUGPRINT_
               Fmt  = '(A,i2,5x,g12.6)'
               Text = 'c**2 is too large,     iVec, c**2 = '
               Write(6,Fmt)Text(1:36),iVec,c2
#endif
*              Go To 520
               Cycle
            End If
            If (ee2<ee1) Then
*-----------   New vector lower eigenvalue.
               ee1   = ee2
               ipBst = iVec
#ifdef _DEBUGPRINT_
               cDotV = c2
#endif
            End If
*
*520        Continue
*
         End Do
*
         If (ipBst.lt.1 .or. ipBst.gt.kOptim) Then
            Write(6,*)' No proper solution found in C2-DIIS !'
            Fmt  = '(6e16.8)'
            Text = 'EigenValues :'
            Call RecPrt(Text,Fmt,EValue,1,kOptim)
            Text = 'EigenVectors :'
            Call RecPrt(Text,Fmt,EVector,kOptim,kOptim)
            Call Quit_OnConvError()
         End If
         call dcopy_(kOptim,EVector(1,ipBst),1,CInter(1,1),1)
*
#ifdef _DEBUGPRINT_
         Write(6,*)
         Write(6,*)' Selected root :',ipBst
         Write(6,'(A,f16.8)')'  c**2 =         ',cDotV
#endif
*
         Call mma_deallocate(EValue)
         Call mma_deallocate(EVector)
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
*------- C1DIIS                                                        *
*                                                                      *
*         References:                                                  *
*         P. Csaszar and P. Pulay, J. Mol. Struc., 114, 31-34 (1984).  *
*         doi:10.1016/S0022-2860(84)87198-7                            *
*                                                                      *
************************************************************************
*                                                                      *
         If (QNRStp) Then
            AccCon = 'QNRc1DIIS'
         Else
            AccCon = 'c1DIIS   '
         End If
*
*        Set up the missing part of the matrix in Eq. (5) and the
*        vector on the RHS in the same equation. Note the sign change!
*
         Do i = 1, kOptim
            Bij(kOptim + 1,i) = - One ! note sign change
            Bij(i,kOptim + 1) = - One ! note sign change
            GDiis(i)          =   Zero
         End Do
         Bij(kOptim + 1,kOptim + 1) =   Zero
         GDiis(kOptim + 1)          = - One  ! note sign change
*
#ifdef _DEBUGPRINT_
         Write(6,*)' B matrix in DIIS_e:'
         Do i = 1, kOptim + 1
            Write(6,'(7f16.8)')(Bij(i,j),j = 1, kOptim + 1),
     &                          GDiis(i)
         End Do
#endif
*
*------- Condition the B matrix
*
         B11 = Sqrt(Bij(1,1)*Bij(kOptim,kOptim))
         Do i = 1, kOptim
            Do j = 1, kOptim
               Bij(i,j) = Bij(i,j)/B11
            End Do
         End Do
*
*------- Solve for the coefficients, solve the equations.
*
         Call Gauss(kOptim + 1,nBij,Bij,CInter(1,1),GDiis)
*
*------- Normalize sum of interpolation coefficients
*
         Fact = Zero
         Do i = 1, kOptim
            Fact = Fact + CInter(i,1)
         End Do
*
         Fact = One/Fact
         Do i = 1, kOptim
            CInter(i,1) = Fact*CInter(i,1)
         End Do
*
*------- Make sure new density gets a weight
         Call C_Adjust(CInter(1,1),kOptim,0.05D0)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*
      Call mma_deallocate(Bij)
*
*     Temporary fix for UHF.
*
      If (nD.eq.2) Call DCopy_(nCI,CInter(1,1),1,CInter(1,2),1)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Fmt  = '(6e16.8)'
      Text = 'The solution vector :'
      Call RecPrt(Text,Fmt,CInter(1,1),1,kOptim)
#endif
*
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 6) = TimFld( 6) + (Cpu2 - Cpu1)
      Return
      End
