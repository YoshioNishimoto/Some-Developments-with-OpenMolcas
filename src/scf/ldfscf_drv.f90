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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************
      Subroutine LDFSCF_Drv(nD,nSym,nBas,DSQ,DLT,DSQ_ab,DLT_ab,FLT,FLT_ab,nFLT,ExFac,nOcc,nOcc_ab)
      use LDFSCF, only: LDF_IntegralMode, LDF_Timing, LDF_IntegralPrescreening, LDF_ContributionPrescreening, &
                        LDF_UseConventionalIntegrals, LDF_UseLDFIntegrals, LDF_UsePSDIntegrals,                          &
                        LDF_UseExactIntegralDiagonal, LDF_IntegralCheck, LDF_FitVerification,                            &
                        LDF_CoefficientCheck, LDF_UBCNorm, LDF_CoulombCheck, LDF_OverlapCheck,                           &
                        LDF_ChargeCheck, LDF_ChargePrint, LDF_ModeCheck, LDF_IntegralPSDCheck, LDF_UseExactIntegrals

      use stdalloc, only: mma_allocate, mma_deallocate
!
!     Thomas Bondo Pedersen, September 2010.
!
!     Purpose: compute two-electron contributions to the SCF Fock matrix
!              using Local Density Fitting (LDF) coefficients.
!
      Use k2_arrays, only: DeDe
      use Constants, only: Zero, Half, One, Two
      Implicit None
      Integer nD, nSym, nFLT
      Integer nBas(nSym), nOcc(nSym), nOcc_ab(nSym)
      Real*8  DSQ(*), DLT(*)
      Real*8  DSQ_ab(*), DLT_ab(*)
      Real*8  FLT(*), FLT_ab(*)
      Real*8  ExFac
#include "localdf.fh"

      Real*8, External::   LDF_Charge, LDF_FittedCharge, dDot_

      Integer, External::  ip_of_Work, iPrintLevel
#if defined (_DEBUGPRINT_)
      Logical, External:: LDF_X_IsSet
#endif

      Character(LEN=6)  RegNam
      Character(LEN=10) SecNam
      Parameter (RegNam='LDFSCF')
      Parameter (SecNam=RegNam//'_Drv')

      Real*8, Parameter:: Tol2C=1.0d-14
      Real*8, Parameter:: Tol_ModeCheck=1.0d-10

      Logical Timing

      Integer Mode, IntegralOption

      Integer nDen_Max
      Parameter (nDen_Max=3)

      Real*8  ThrPS(2)
      Real*8  FactC(nDen_Max)
      Integer ip_D(nDen_Max)

      Integer nDen
      Integer u, v, uv
      Integer iPrint
      Integer irc
      Integer lF
      Integer iDen
      Integer AB_MAE, AB_MRNrm
      Integer lMode
      Integer n

      Logical Add
      Logical PackedD
      Logical PackedF
      Logical ComputeF

      Logical QuitOnError
      Logical Silent

      Real*8  MAE, MRNrm, Q, deltaQ
      Real*8  Q_LDF
      Real*8  MaxErr, MaxAbsErr, AverageErr, RMSErr
      Real*8  MaxErr_OffD, MaxAbsErr_OffD
      Real*8  DiffNorm, DiffSum
      Real*8  ConvNorm, LDFNorm
      Real*8  ConvSum, LDFSum
      Real*8  FNorm, factor

      Integer i, j
      Integer iTri
      Real*8, Allocatable:: UBFNorm(:), DrvF(:,:)

      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

#if defined (_DEBUGPRINT_)
      ! Enter
#endif
      Call Untested('LDF-SCF')

      ! Get print level
#if defined (_DEBUGPRINT_)
      iPrint=max(iPrintLevel(-1),4)
#else
      iPrint=iPrintLevel(-1)
#endif

      ! Print entry message
      If (iPrint.ge.3) Then
         Call Cho_Head('Enter '//SecNam,'-',80,6)
      End If

      ! Symmetry not implemented for LDF
      If (nSym.ne.1) Then
         Call WarningMessage(2,SecNam//': Symmetry not implemented for LDF')
         Call LDF_Quit(1)
      End If
!
      Call Set_Basis_Mode('WithAuxiliaryr')
      Call Setup_iSD()

#if defined (_DEBUGPRINT_)
      If (.not.LDF_X_IsSet()) Then
         Call WarningMessage(2,SecNam//': LDF info not set!')
         Call LDF_Quit(1)
      End If
#endif

      ! Debug option: print charge
      If (LDF_ChargePrint) Then
         PackedD=.False.
         ip_D(1)=ip_of_Work(DSQ(1))
         Q=LDF_Charge(PackedD,ip_D(1))
         Q_LDF=LDF_FittedCharge(PackedD,ip_D(1))
         Write(6,'(/,A,A,1P,D20.10,2(2X,A,D20.10))')SecNam,': Q=',Q,'Q_LDF=',Q_LDF,'Error=',Q_LDF-Q
         Call xFlush(6)
         If (Q.lt.Zero .or. Q_LDF.lt.Zero) Then
            Call WarningMessage(2,SecNam//': this is unphysical!')
            Call LDF_Quit(1)
         End If
      End If

      Call mma_allocate(DeDe,[-1,-1],label='Dede') ! Dummy allocation

      ! Integral representation?
      Mode=LDF_IntegralMode
      ! Time the Fock build?
      Timing=LDF_Timing .or. iPrint.ge.3
      ! Set integral handling (options)
      If (LDF_UseLDFIntegrals) Then
         IntegralOption=111
      Else If (LDF_UseConventionalIntegrals) Then
         IntegralOption=222
      Else If (LDF_UsePSDIntegrals) Then
         IntegralOption=333
      Else If (LDF_UseExactIntegralDiagonal) Then
         IntegralOption=444
      Else
         IntegralOption=0
      End If
      ! Determine integral and contribution prescreening thresholds
      If (LDF_IntegralPrescreening.lt.Zero) Then
         LDF_IntegralPrescreening=min(abs(Thr_Accuracy)*1.0d-2,1.0d-8)
         If (iPrint.ge.4) Then
            Write(6,'(A,1P,D15.6)')'Setting default Integral Prescreeening threshold:    ',LDF_IntegralPrescreening
            Call xFlush(6)
         End If
      End If
      If (LDF_ContributionPrescreening.lt.Zero) Then
         LDF_ContributionPrescreening=min(abs(Thr_Accuracy)*1.0d-1,1.0d-6)
         If (iPrint.ge.4) Then
            Write(6,'(A,1P,D15.6)')'Setting default Contribution Prescreeening threshold:',LDF_ContributionPrescreening
            Call xFlush(6)
         End If
      End If
      ! Set prescreening thresholds to pass to Fock calculator
      ThrPS(1)=LDF_IntegralPrescreening
      ThrPS(2)=LDF_ContributionPrescreening

      ! Debug: check that full integral matrix is PSD (if requested)
      ! Note: this is only done in first call to this routine
      !       (unless flag is reset outside)
      If (LDF_IntegralPSDCheck.gt.0) Then
         Call WarningMessage(0,SecNam//': Checking full integral matrix PSD')
         If (LDF_UseExactIntegrals.eq.1) Then
            Write(6,'(A)')'Using exact diagonal blocks of the integral matrix'
         Else If (LDF_UseExactIntegrals.eq.2) Then
            Write(6,'(A)')'Using exact off-diagonal blocks of the integral matrix'
         End If
         Call xFlush(6)
         Call LDF_CheckPSD_Full(LDF_IntegralPSDCheck.eq.1,Mode,LDF_UseExactIntegrals,ThrPS(1),irc)
         If (irc.ne.0) Then
            Call WarningMessage(0,SecNam//': full integral matrix NOT PSD')
         Else
            Call WarningMessage(0,SecNam//': full integral matrix PSD')
         End If
         Call xFlush(6)
         LDF_IntegralPSDCheck=-1
      End If

      ! Debug: check fitting coefficients (if requested)
      ! Note: this is only done in first call to this routine
      !       (unless flag is reset outside)
      If (LDF_CoefficientCheck) Then
         Call WarningMessage(0,SecNam//': Checking coefficients')
         Call xFlush(6)
         Call LDF_CheckAllC_wLD(iPrint.ge.3)
         Write(6,'(/,A)') 'Coefficient check completed: all OK!'
         Call xFlush(6)
         LDF_CoefficientCheck=.False.
      End If

      ! Debug: check overlap integrals (if requested)
      ! Note: this is only done in first call to this routine
      !       (unless flag is reset outside)
      If (LDF_OverlapCheck) Then
         Call WarningMessage(0,SecNam//': Checking overlap integrals')
         Call LDF_CheckAllOverlapIntegrals(iPrint.ge.3,Tol2C,MAE,AB_MAE,MRNrm,AB_MRNrm)
         If (iPrint.lt.3) Then
            Write(6,'(/,2X,A,1P,D20.10,2X,A)')'Tolerance for 2C errors..',Tol2C,'(all 2C passed)'
            Write(6,'(2X,A,1P,D20.10,2X,A,I10)')'Max abs error............',MAE,'@AB=',AB_MAE
            Write(6,'(2X,A,1P,D20.10,2X,A,I10)')'Max relative norm error..',MRNrm,'@AB=',AB_MRNrm
         End If
         Write(6,'(/,A)') 'Overlap check completed: 2C OK!'
         Call xFlush(6)
         LDF_OverlapCheck=.False.
      End If

      ! Debug: verify fit for each atom pair (if requested)
      ! Note: this is only done in first call to this routine
      !       (unless flag is reset outside)
      If (LDF_FitVerification) Then
         Call WarningMessage(0,SecNam//': Verifying fit')
         Call xFlush(6)
         Call LDF_VerifyFit_Drv(irc)
         If (irc.ne.0) Then
            Call WarningMessage(2,SecNam//': Fit verification failed!')
            Write(6,'(A,I10)')'LDF_VerifyFit_Drv returned code',irc
            Call LDF_Quit(1)
         Else
            Write(6,'(/,A)') 'Fit verification completed: all OK!'
            Call xFlush(6)
         End If
         LDF_FitVerification=.False.
      End If

      ! Debug: check all integrals (if requested)
      ! Note: this is only done in first call to this routine
      !       (unless flag is reset outside)
      If (LDF_IntegralCheck) Then
         Call WarningMessage(0,SecNam//': Checking all integrals - this may take a while!')
         Call xFlush(6)
         QuitOnError=.True.
         Silent=iPrint.lt.3
         Call LDF_CheckAllIntegrals(QuitOnError,Silent,Mode,ThrPS(1),MaxErr,MaxAbsErr,AverageErr,RMSErr,    &
                                    MaxErr_OffD, MaxAbsErr_OffD,DiffNorm,DiffSum,ConvNorm,LDFNorm,          &
                                    ConvSum,LDFSum)
         If (Silent) Then
            Call Cho_Head('LDF Integral Error Analysis','-',80,6)
            Write(6,'(3X,A,A)')'(Increase print level to get more details on ','individual atom pair blocks)'
            If (Mode.eq.1) Then
               Write(6,'(3X,A)') 'LDF integral representation: ROBUST'
            Else If (Mode.eq.2) Then
               Write(6,'(3X,A)') 'LDF integral representation: NONROBUST'
            Else If (Mode.eq.3) Then
               Write(6,'(3X,A)') 'LDF integral representation: HALF-AND-HALF'
            Else
               Write(6,'(3X,A)') 'LDF integral representation: UNKNOWN'
            End If
            Write(6,'(3X,A,1P,D20.10)') 'LDF integral prescreening threshold:',ThrPS(1)
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)') 'Conventional norm...',ConvNorm, 'Conventional sum....',ConvSum
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)') 'LDF norm............',LDFNorm, 'LDF sum.............',LDFSum
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)') 'Difference norm.....',DiffNorm, 'Difference sum......',DiffSum
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)') 'Max Error...........',MaxErr, 'Max Abs Error.......',MaxAbsErr
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)') 'Max OffD Error......',MaxErr_OffD, 'Max Abs OffD Error..',MaxAbsErr_OffD
            Write(6,'(3X,A,1P,D20.10,3X,A,D20.10)') 'Average Error.......',AverageErr, 'RMS Error...........',RMSErr
         End If
         If (MaxAbsErr.gt.Thr_Accuracy) Then
            Write(6,'(/,A)') 'Integral check completed:'
            Call WarningMessage(1,'Integral check completed:Max Abs Error is greater than LDF target accuracy!')
         Else
            Write(6,'(/,A)') 'Integral check completed: all OK!'
            Call xFlush(6)
         End If
         Call xFlush(6)
         LDF_IntegralCheck=.False. ! Done, so do not do it again
      End If

!--------------------------------------------------------------
!     Branch according to Coulomb-and-exchange or Coulomb-only.
!--------------------------------------------------------------

      If (ExFac.ne.Zero) Then ! Coulomb-and-exchange
         Write(6,'(//,A,A)') SecNam,': Exchange not implemented yet!'
         Call LDF_NotImplemented()
      Else ! Coulomb-only
         ! Set number of density matrices
         nDen=1
         ! Defensive coding: make sure nDen_Max has not changed
         ! below nDen
         If (nDen.gt.nDen_Max) Then
            Call WarningMessage(2,SecNam//': nDen>nDen_Max')
            Call LDF_Quit(1)
         End If
         ! Set scaling factors
         FactC(1)=One
         If (nD==1) Then ! spin-restricted Coulomb-only
            ! Get pointers to D and F
            ! Off-diagonal elements of DLT are scaled by 2 by the SCF
            ! program. This is incompatible with the LDF implementation.
            ! So, use instead the quadratic DSQ.
            ! Use packed F (FLT) for result.
            ip_D(1)=ip_of_Work(DSQ(1))
            ! Set flags for quadratic (SQ) density matrix and packed
            ! Fock matrix
            PackedD=.False.
            PackedF=.True.
            ! Do not add result to F, simply put it there
            Add=.False.
            ! Debug: check charge (if requested)
            If (LDF_ChargeCheck) Then
               Call WarningMessage(0,SecNam//': Checking charge')
               Call xFlush(6)
               Call LDF_CheckCharge(.True.,PackedD,ip_D(1),MAE,AB_MAE,Q,deltaQ)
            End If
            ! Debug: compute norm of upper bound Fock matrix error for
            !        the Coulomb contribution (if requested).
            If (LDF_UBCNorm) Then
               Call WarningMessage(0,SecNam//': Computing norm of upper bound to Coulomb error')
               Call xFlush(6)
               Call mma_Allocate(UBFNorm,nDen,Label='UBFNorm')
               Call LDF_Fock_CoulombUpperBoundNorm_Full(.True.,PackedD,nDen,FactC,ip_D,UBFNorm)
               Call mma_deAllocate(UBFNorm)
            End If
            ! Print args
            If (iPrint.ge.3) Then
               Write(6,'(A)') 'Calling RHF Fock matrix calculator with args'
               Write(6,'(2X,A,I15,2X,A,14X,L1,2X,A,I15)') 'IntegralOption...',IntegralOption,'Timing...........',Timing,  &
               'Mode.............',Mode
               Write(6,'(2X,A,1P,2D15.6)') 'ThrPS............',ThrPS(1),ThrPS(2)
               Write(6,'(2X,A,14X,L1,2X,A,14X,L1,2X,A,14X,L1)')'Add..............',Add, 'PackedD..........',PackedD,      &
               'PackedF..........',PackedF
               Write(6,'(2X,A,I15)') 'nDen.............',nDen
               Write(6,'(2X,A,1P,3D15.6)') 'FactC............',(FactC(i),i=1,nDen)
               If (iPrint.ge.5) Then
                  Write(6,'(2X,A,3I15)') 'ip_D.............',(ip_D(i),i=1,nDen)
               End If
            End If
            ! Compute two-electron contributions to Fock matrix
            Call LDF_Fock_CoulombOnly(IntegralOption,Timing,Mode,ThrPS,Add,PackedD,PackedF,nDen,FactC,ip_D,FLT)
            ! Debug: check Coulomb error
            If (LDF_CoulombCheck) Then
               Call WarningMessage(0,SecNam//': Analysis of Coulomb error')
               Call xFlush(6)
               If (PackedF) Then
                  lF=nBas(1)*(nBas(1)+1)/2
               Else
                  lF=nBas(1)**2
               End If
               Call mma_allocate(DrvF,lF,nDen,Label='DrvF')
               Do iDen=1,nDen
                  Call dCopy_(lF,FLT(iDen),1,DrvF(:,iDen),1)
               End Do
               ComputeF=.False.
               Call LDF_Fock_CoulombErrorAnalysis(ComputeF,Mode,PackedD,PackedF,nDen,FactC,ip_D,DrvF)
               Call mma_deallocate(DrvF)
            End If
            ! Debug: check mode consistency
            If (LDF_ModeCheck) Then
               LDF_ModeCheck=.False. ! check only once
               Call WarningMessage(0,SecNam//': Mode Check')
               Call xFlush(6)
               If (PackedF) Then
                  lF=nBas(1)*(nBas(1)+1)/2
               Else
                  lF=nBas(1)**2
               End If
               Call mma_allocate(DrvF,lF,2*nDen,Label='DrvF')
               Do iDen=1,nDen
                  Call dCopy_(lF,FLT(iDen),1,DrvF(:,nDen+iDen),1)
               End Do
               If (Mode.eq.1) Then
                  lMode=3
                  factor=-Two
               Else If (Mode.eq.2) Then
                  lMode=3
                  factor=-Two
               Else If (Mode.eq.3) Then
                  Do iDen=1,nDen
                     Call dScal_(lF,Two,DrvF(:,nDen+iDen),1)
                  End Do
                  lMode=1
                  factor=-One
               Else
                  Call WarningMessage(2,SecNam//': unknown Mode')
                  Call LDF_Quit(1)
                  lMode=0
                  factor=Zero
               End If
               Call LDF_Fock_CoulombOnly(IntegralOption,Timing,lMode,ThrPS,Add,PackedD,PackedF,nDen,FactC,ip_D,DrvF)
               Do iDen=1,nDen
                  Call dAXPY_(lF,factor,DrvF(:,     iDen),1,DrvF(:,nDen+iDen),1)
               End Do
               If (Mode.eq.1) Then
                  lMode=2
                  factor=One
               Else If (Mode.eq.2) Then
                  lMode=1
                  factor=One
               Else If (Mode.eq.3) Then
                  lMode=2
                  factor=-One
               Else
                  Call WarningMessage(2,SecNam//': unknown Mode')
                  Call LDF_Quit(1)
                  lMode=0
                  factor=Zero
               End If
               Call LDF_Fock_CoulombOnly(IntegralOption,Timing,lMode,ThrPS,Add,PackedD,PackedF,nDen,FactC,ip_D,DrvF)
               Do iDen=1,nDen
                  Call dAXPY_(lF,factor,DrvF(:,     iDen),1,DrvF(:,nDen+iDen),1)
               End Do
               Call Cho_Head(SecNam//': Mode Check','=',80,6)
               n=0
               Do iDen=1,nDen
                  FNorm=sqrt(dDot_(lF,DrvF(:,nDen+iDen),1,DrvF(:,nDen+iDen),1))
                  If (FNorm.gt.Tol_ModeCheck) Then
                     Write(6,'(3X,A,I3,A,1P,D20.10,A)')'Density no.',iDen,' Check norm=',Fnorm,'  (FAIL)'
                     n=n+1
                  Else
                     Write(6,'(3X,A,I3,A,1P,D20.10,A)')'Density no.',iDen,' Check norm=',Fnorm,'  (pass)'
                  End If
               End Do
               If (n.ne.0) Then
                  Call WarningMessage(2,SecNam//': Mode Check: Failed')
                  Call LDF_Quit(1)
               End If
               Call xFlush(6)
               Call mma_deallocate(DrvF)
            End If
         Else ! spin-unrestricted Coulomb-only
            ! Add alpha and beta parts of density matrix
            Call dAXPY_(nBas(1)**2,One,DSQ_ab,1,DSQ,1)
            ! Get pointers to D and F
            ! Off-diagonal elements of DLT are scaled by 2 by the SCF
            ! program. This is incompatible with the LDF implementation.
            ! So, use instead the quadratic DSQ.
            ! Use packed F (FLT) for result.
            ip_D(1)=ip_of_Work(DSQ(1))
            ! Set flags for quadratic (SQ) density matrix and packed
            ! Fock matrix
            PackedD=.False.
            PackedF=.True.
            ! Do not add result to F, simply put it there
            Add=.False.
            ! Debug: check charge (if requested)
            If (LDF_ChargeCheck) Then
               Call WarningMessage(0,SecNam//': Checking charge')
               Call xFlush(6)
               Call LDF_CheckCharge(.True.,PackedD,ip_D(1),MAE,AB_MAE,Q,deltaQ)
            End If
            ! Debug: compute norm of upper bound Fock matrix error for
            !        the Coulomb contribution (if requested).
            If (LDF_UBCNorm) Then
               Call WarningMessage(0,SecNam//': Computing norm of upper bound to Coulomb error')
               Call xFlush(6)
               Call mma_Allocate(UBFNorm,nDen,Label='UBFNorm')
               Call LDF_Fock_CoulombUpperBoundNorm_Full(.True.,PackedD,nDen,FactC,ip_D,UBFNorm)
               Call mma_deAllocate(UBFNorm)
            End If
            ! Print args
            If (iPrint.ge.3) Then
               Write(6,'(A)')'Calling UHF Fock matrix calculator with args'
               Write(6,'(2X,A,I15,2X,A,14X,L1,2X,A,I15)')'IntegralOption...',IntegralOption,'Timing...........',Timing,     &
               'Mode.............',Mode
               Write(6,'(2X,A,1P,2D15.6)') 'ThrPS............',ThrPS(1),ThrPS(2)
               Write(6,'(2X,A,14X,L1,2X,A,14X,L1,2X,A,14X,L1)')'Add..............',Add,'PackedD..........',PackedD,         &
               'PackedF..........',PackedF
               Write(6,'(2X,A,I15)') 'nDen.............',nDen
               Write(6,'(2X,A,1P,3D15.6)') 'FactC............',(FactC(i),i=1,nDen)
               If (iPrint.ge.5) Then
                  Write(6,'(2X,A,3I15)') 'ip_D.............',(ip_D(i),i=1,nDen)
               End If
            End If
            ! Compute two-electron contributions to Fock matrix
            Call LDF_Fock_CoulombOnly(IntegralOption,Timing,Mode,ThrPS,Add,PackedD,PackedF,nDen,FactC,ip_D,FLT)
            ! Debug: check Coulomb error
            If (LDF_CoulombCheck) Then
               Call WarningMessage(2,SecNam//': Analysis of Coulomb error')
               Call xFlush(6)
               If (PackedF) Then
                  lF=nBas(1)*(nBas(1)+1)/2
               Else
                  lF=nBas(1)**2
               End If
               Call mma_Allocate(DrvF,lF,nDen,Label='DrvF')
               Do iDen=1,nDen
                  Call dCopy_(lF,FLT(iDen),1,DrvF(:,iDen),1)
               End Do
               ComputeF=.False.
               Call LDF_Fock_CoulombErrorAnalysis(ComputeF,Mode,PackedD,PackedF,nDen,FactC,ip_D,DrvF)
               Call mma_deallocate(DrvF)
            End If
            ! Debug: check mode consistency
            If (LDF_ModeCheck) Then
               LDF_ModeCheck=.False. ! check only once
               Call WarningMessage(0,SecNam//': Mode Check')
               Call xFlush(6)
               If (PackedF) Then
                  lF=nBas(1)*(nBas(1)+1)/2
               Else
                  lF=nBas(1)**2
               End If
               Call mma_Allocate(DrvF,lF,2*nDen,Label='DrvF')
               Do iDen=1,nDen
                  Call dCopy_(lF,FLT(iDen),1,DrvF(:,nDen+iDen),1)
               End Do
               If (Mode.eq.1) Then
                  lMode=3
                  factor=-Two
               Else If (Mode.eq.2) Then
                  lMode=3
                  factor=-Two
               Else If (Mode.eq.3) Then
                  Do iDen=1,nDen
                     Call dScal_(lF,Two,DrvF(:,nDen+iDen),1)
                  End Do
                  lMode=1
                  factor=-One
               Else
                  Call WarningMessage(2,SecNam//': unknown Mode')
                  Call LDF_Quit(1)
                  lMode=0
                  factor=Zero
               End If
               Call LDF_Fock_CoulombOnly(IntegralOption,Timing,lMode,ThrPS,Add,PackedD,PackedF,nDen,FactC,ip_D,DrvF)
               Do iDen=1,nDen
                  Call dAXPY_(lF,factor,DrvF(:,     iDen),1,DrvF(:,nDen+iDen),1)
               End Do
               If (Mode.eq.1) Then
                  lMode=2
                  factor=One
               Else If (Mode.eq.2) Then
                  lMode=1
                  factor=One
               Else If (Mode.eq.3) Then
                  lMode=2
                  factor=-One
               Else
                  Call WarningMessage(2,SecNam//': unknown Mode')
                  Call LDF_Quit(1)
                  lMode=0
                  factor=Zero
               End If
               Call LDF_Fock_CoulombOnly(IntegralOption,Timing,lMode,ThrPS,Add,PackedD,PackedF,nDen,FactC,ip_D,DrvF)
               Do iDen=1,nDen
                  Call dAXPY_(lF,factor,DrvF(:,     iDen),1,DrvF(:,nDen+iDen),1)
               End Do
               Call Cho_Head(SecNam//': Mode Check','=',80,6)
               n=0
               Do iDen=1,nDen
                  FNorm=sqrt(dDot_(lF,DrvF(:,nDen+iDen),1,DrvF(:,nDen+iDen),1))
                  If (FNorm.gt.Tol_ModeCheck) Then
                     Write(6,'(3X,A,I3,A,1P,D20.10,A)')'Density no.',iDen,' Check norm=',Fnorm,'  (FAIL)'
                     n=n+1
                  Else
                     Write(6,'(3X,A,I3,A,1P,D20.10,A)')'Density no.',iDen,' Check norm=',Fnorm,'  (pass)'
                  End If
               End Do
               If (n.ne.0) Then
                  Call WarningMessage(2,SecNam//': Mode Check: Failed')
                  Call LDF_Quit(1)
               End If
               Call xFlush(6)
               Call mma_deallocate(DrvF)
            End If
            ! Copy result to FLT_ab
            Call dCopy_(nFLT,FLT,1,FLT_ab,1)
            ! Restore DSQ from DLT (which is scaled by two on the
            ! off-diagonal)
            uv=0
            Do v=1,nBas(1)
               Do u=1,v-1
                  uv=uv+1
                  DSQ(uv)=half*DLT(iTri(u,v))
               End Do
               uv=uv+1
               DSQ(uv)=DLT(iTri(v,v))
               Do u=v+1,nBas(1)
                  uv=uv+1
                  DSQ(uv)=half*DLT(iTri(u,v))
               End Do
            End Do
         End If ! spin-restricted or spin-unrestricted
      End If ! Coulomb-and-exchange or Coulomb-only

      Call mma_deallocate(DeDe)

      ! Print exit message
      If (iPrint.ge.3) Then
         Write(6,'(2X,A,A)') 'Exit ',SecNam
         Call xFlush(6)
      End If
      Call Free_iSD()

#if defined (_DEBUGPRINT_)
      ! Exit
#endif
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(DLT_ab)
         Call Unused_integer_array(nOcc)
         Call Unused_integer_array(nOcc_ab)
      End If
      End Subroutine LDFSCF_Drv
