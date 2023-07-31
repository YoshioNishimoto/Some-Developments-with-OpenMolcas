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
      Subroutine Convrg(iter,kIter, nInter,iStop,MxItr,
     &                  mIntEff,mTtAtm,GoOn,Step_Trunc,
     &                  Just_Frequencies)
      Use Chkpnt
      Use Slapaf_Info, only: Cx, Gx, Coor, GNrm, Energy, Shift, qInt,
     &                       dqInt, Lbl
      use Slapaf_Parameters, only: HUpMet, FindTS, Analytic_Hessian,
     &                             MaxItr, Numerical, iNeg, GrdMax,
     &                             E_Delta, ThrEne, ThrGrd, nLambda,
     &                             iOptC, ThrCons, ThrMEP, Baker,
     &                             eMEPTest, rMEP, MEP, nMEP, Stop,
     &                             NADC, EDiffZero, ApproxNADC
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "warnings.h"
      Integer:: IRC=0
      Real*8 Maxed, MaxErr
      Character(LEN=5) ConLbl(5)
      Character(LEN=1) Step_Trunc
      Character(LEN=16) StdIn
      Character(LEN=80) Point_Desc
      Character(LEN=16) MEP_Text
      Logical Conv1, GoOn, Found, Terminate, Last_Energy,
     &        Just_Frequencies, Saddle, eTest,
     &        IRCRestart, Conv2, ConvTmp, TSReg, BadConstraint,
     &        TurnBack
      Character(LEN=8) Temp
      Real*8, Allocatable:: Coor1(:,:), Coor2(:,:)
      Real*8, Allocatable:: E_IRC(:), C_IRC(:,:,:), G_IRC(:,:,:)
      Real*8, Allocatable:: E_S(:), C_S(:,:,:), G_S(:,:,:)
      Real*8, Allocatable:: E_R(:), C_R(:,:), G_R(:,:)
      Real*8, Allocatable:: E_P(:), C_P(:,:), G_P(:,:)
      Real*8, Allocatable:: E_MEP(:), G_MEP(:,:,:)
      Real*8, Allocatable, Target:: C_MEP(:,:,:)
      Real*8, Allocatable:: L_MEP(:), Cu_MEP(:)
      Integer, Allocatable:: Information(:)
      Real*8, Allocatable:: Tmp(:)
      Real*8, Allocatable, Target:: Not_Allocated(:,:), OfRef(:,:)
      Real*8 rDum(1,1,1,1)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      Subroutine SphInt(xyz,nCent,OfRef,RR0,Bf,l_Write,Label,dBf,ldB)
      Integer nCent
      Real*8  xyz(3,nCent)
      Real*8, Allocatable, Target:: OfRef(:,:)
      Real*8  RR0
      Real*8  Bf(3,nCent)
      Logical l_Write
      Character(LEN=8) Label
      Real*8  dBf(3,nCent,3,nCent)
      Logical ldB
      End Subroutine SphInt
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      nAtom=SIZE(Cx,2)
      TSReg = iAnd(iOptC,8192).eq.8192
*                                                                      *
************************************************************************
*                                                                      *
      Lu=6
      nSaddle_Max=100
      iRout=116
      iPrint=nPrint(iRout)
      If (iPrint.ge.99) Then
         Call RecPrt('Convrg: Energy',' ',Energy,1,iter)
         Call RecPrt('Convrg: Shift',' ',Shift,nInter,iter)
         Call RecPrt('Convrg: qInt',' ',qInt,nInter,iter+1)
         Call RecPrt('Convrg: dqInt',' ',dqInt,nInter,iter)
         Call RecPrt('Convrg: Cx',' ',Cx,3*nAtom,iter+1)
         Call RecPrt('Convrg: Gx',' ',Gx,3*nAtom,iter+1)
      End If
*
      Call Get_iScalar('Saddle Iter',iter_S)
      If (iter_S.eq.0) Then
         iter_S=1
         Call Put_iScalar('Saddle Iter',iter_S)
         Call f_Inquire('RUNFILE2',Found)
         If (Found) Then
            Call NameRun('RUNFILE2')
            Call Put_iScalar('Saddle Iter',iter_S)
            Call NameRun('#Pop')
         End If
      End If
      Temp=' '
      If (Analytic_Hessian) Then
         If (HUPMET.eq.'  No  '.or.HUPMET.eq.' None ') Then
*           Temp='Analytic'
            Temp='Computed'
         Else
            Temp(1:6)= HUPMET(1:6)
         End if
      Else
         Temp(1:6)= HUPMET(1:6)
      End If
*
      Call mma_allocate(Coor1,3,mTtAtm,Label='Coor1')
      Call mma_allocate(Coor2,3,mTtAtm,Label='Coor2')
      Call AtmLst(Cx(:,:,iter  ),nAtom,Coor1,mTtAtm)
      Call AtmLst(Cx(:,:,iter+1),nAtom,Coor2,mTtAtm)
      Call OptRMS_Slapaf(Coor1,Coor2,mTtAtm,RMS,RMSMax)
      Call mma_deallocate(Coor1)
      Call mma_deallocate(Coor2)
*
      If (kIter.ne.iter .and. kIter.eq.1) Then
         Fabs   = GNrm(kiter)
         E = Energy(kiter)
      Else
         Fabs   = GNrm(iter)
         E = Energy(iter)
      End If
      Fabs = Max(Zero,Fabs)
      E0 = E + E_delta
      Energy(iter+1)=E0
      Gx(:,:,iter+1)=Zero
      If (kiter.eq.1) Then
         eChng=Zero
      Else
         If (kIter.ne.iter .and. kIter.eq.2) Then
            eChng=Energy(iter)-Energy(1)
         Else
            eChng=Energy(iter)-Energy(iter-1)
         End If
      End If
*
      eDiffMEP=Zero
      If (MEP.or.rMEP) Then
         Saddle=.False.
         iMEP=0
         Call Qpg_iScalar('nMEP',Found)
         If (Found) Call Get_iScalar('nMEP',iMEP)
         If (iMEP.eq.0) Then
            iOff_Iter=0
            Call Put_iScalar('iOff_Iter',iOff_Iter)
         Else
            Call Get_iScalar('iOff_Iter',iOff_Iter)
         End If
      Else
         iOff_Iter=0
         iSaddle=0
         Call qpg_dArray('Saddle',Saddle,nSaddle)
         Saddle=Saddle.and..NOT.Just_Frequencies
         If (Saddle) Then
            Call Qpg_iScalar('nMEP',Found)
            If (Found) Call Get_iScalar('nMEP',iSaddle)
            Call Get_iScalar('iOff_Iter',iOff_Iter)
         End If
      End If
*
*----- Convergence criteria
*
*      Too many iterations
*
*      or
*
* 1)   a la Baker
*      Abs(GrdMax).lt.ThrGrd
*      and
*      Abs(eChng).lt.ThrEne or RMSMax.lt.ThrGrd
*
* 2)   a la Gaussian
*      Abs(Fabs/Sqrt(mIntEff)).lt.ThrGrd
*      and
*      Abs(GrdMax).lt.ThrGrd*1.5
*      and
*      ((RMS.lt.ThrGrd*4
*        and
*        RMSMax.lt.ThrGrd*6)
*       or
*       Abs(eChng).lt.ThrEne)
*
      If (Baker) Then
         Val1= Abs(eChng)
         Thr1= ThrEne
         If (kIter.le.1) Then
            ConLbl(1)=' --- '
         Else If (Val1.lt.Thr1) Then
            ConLbl(1)=' Yes '
         Else
            ConLbl(1)=' No  '
         End If
         Val2= RMSMax
         Thr2= ThrGrd
         If (Val2.lt.Thr2) Then
            If (Step_Trunc.ne.' ') Then
               ConLbl(2)=' No *'
            Else
               ConLbl(2)=' Yes '
            End If
         Else
            ConLbl(2)=' No  '
         End If
         Val3= Abs(GrdMax)
         Thr3= ThrGrd
         If (Val3.lt.Thr3) Then
            ConLbl(3)=' Yes '
         Else
            ConLbl(3)=' No  '
         End If
         Conv1= Val1.lt.Thr1.and.kIter.gt.1
         Conv1= Conv1.or. (Val2.lt.Thr2 .and. Step_Trunc.eq.' ')
         Conv1= Conv1.and. Val3.lt.Thr3
      Else
         Val2=Abs(Fabs/Sqrt(DBLE(mIntEff)))
         Thr2=ThrGrd
         Conv1=Val2.lt.Thr2
         If (Conv1) Then
            ConLbl(2)=' Yes '
         Else
            ConLbl(2)=' No  '
         End If
         Conv1= Conv1.and.Abs(GrdMax).lt.ThrGrd*1.5D0
         Val4=Abs(GrdMax)
         Thr4=ThrGrd*1.5D0
         ConvTmp=Val4.lt.Thr4
         Conv1=Conv1.and.ConvTmp
         If (ConvTmp) Then
            ConLbl(4)=' Yes '
         Else
            ConLbl(4)=' No  '
         End If
         Conv2= RMS.lt.ThrGrd*4.D0 .and. Step_Trunc.eq.' '
         Val1=RMS
         Thr1=ThrGrd*4.0D0
         ConvTmp=Val1.lt.Thr1
         Conv2=ConvTmp .and. Step_Trunc.eq.' '
         If (ConvTmp) Then
            If (Step_Trunc.ne.' ') Then
               ConLbl(1)=' No *'
            Else
               ConLbl(1)=' Yes '
            End If
         Else
            ConLbl(1)=' No  '
         End If
         Val3=RMSMax
         Thr3=ThrGrd*6.0D0
         ConvTmp=Val3.lt.Thr3
         Conv2=Conv2.and.ConvTmp
         If (ConvTmp) Then
            If (Step_Trunc.ne.' ') Then
               ConLbl(3)=' No *'
            Else
               ConLbl(3)=' Yes '
            End If
         Else
            ConLbl(3)=' No  '
         End If
         Val5=Abs(eChng)
         Thr5=ThrEne
         ConvTmp=Val5.lt.Thr5 .and. kIter.gt.1
         If (ConvTmp) Then
            ConLbl(5)=' Yes '
         Else
            If (kIter.gt.1) Then
               ConLbl(5)=' No  '
            Else
               ConLbl(5)=' --- '
            End If
         End If
         Conv2=Conv2.or.ConvTmp
         Conv1=Conv1.and.Conv2
      End If
*
      Stop = Conv1 .or. (kIter-iOff_Iter).ge.MxItr ! CGG
      iStop=1
      If (kIter-iOff_Iter.ge.MxItr) iStop=16       ! CGG
      If (Conv1.or.Just_Frequencies)  iStop= 0
*
      If (GoOn) Then
         Stop=.False.
         iStop=1
      Else
         If (Just_Frequencies) Stop=.True.
         nPrint(52)=nPrint(52)+1
         nPrint(54)=nPrint(54)+1
         iPrint    =iPrint+1
         nPrint(53)=nPrint(53)+1
      End If
      If (.Not.Just_Frequencies)
     &   Call Status(kIter-iOff_Iter,E,Fabs,E0,MaxItr-1,eChng,Temp,
     &               Step_Trunc,.NOT.Numerical)
*
      If (Baker) Then
         If (iPrint.ge.5) Then
            Write (Lu,'(A)')
     &        '                +----------------------------------+'
            Write (Lu,'(A)')
     &        '                +  Value      Threshold Converged? +'
            Write (Lu,'(A)')
     &        '+---------------+----------------------------------+'
            Write (Lu,4)
     &        '+ Max. gradient +',Val3,' ',Thr3,'    ',ConLbl(3),'  +'
            Write (Lu,'(A)')
     &        '+---------------+----------------------------------+'
            Write (Lu,4)
     &        '+ Max. disp.    +',Val2,' ',Thr2,'    ',ConLbl(2),'  +'
            Write (Lu,'(A)')
     &        '+---------------+----------------------------------+'
            Write (Lu,4)
     &        '+ Energy diff.  +',Val1,' ',Thr1,'    ',ConLbl(1),'  +'
            Write (Lu,'(A)')
     &        '+---------------+----------------------------------+'
 4          Format(A,2(ES11.4,A),A,A)
         End If
      Else
         If (iPrint.ge.5) Then
            Write (Lu,'(A)')
     &        '       +----------------------------------+'
     &      //'----------------------------------+'
            Write (Lu,'(A)')
     &        '       +    Cartesian Displacements       +'
     &      //'    Gradient in internals         +'
            Write (Lu,'(A)')
     &        '       +  Value      Threshold Converged? +'
     &      //'  Value      Threshold Converged? +'
            Write (Lu,'(A)')
     &        ' +-----+----------------------------------+'
     &      //'----------------------------------+'
            Write (Lu,5)
     &        ' + RMS +',Val1,   ' ',Thr1,'    ',ConLbl(1),  '  +',
     &                Val2,   ' ',Thr2,'    ',ConLbl(2),  '  +'
            Write (Lu,'(A)')
     &       ' +-----+----------------------------------+'
     &          //'----------------------------------+'
            Write (Lu,5)
     &        ' + Max +',Val3,   ' ',Thr3,'    ',ConLbl(3),  '  +',
     &                Val4,   ' ',Thr4,'    ',ConLbl(4),  '  +'
            Write (Lu,'(A)')
     &       ' +-----+----------------------------------+'
     &      //'----------------------------------+'
            If (ThrEne.gt.Zero) Then
               Write (Lu,5)
     &           ' + dE  +',Val5,   ' ',Thr5,'    ',ConLbl(5),  '  +'
               Write (Lu,'(A)')
     &          ' +-----+----------------------------------+'
            End If
            Write (Lu,*)
 5          Format(A,2(2(ES11.4,A),A,A))
         End If
      End If
*
      If (Stop.and.Conv1) Then
         Call Qpg_dScalar('Max error',Found)
         If (Found) Call Get_dScalar('Max error',MaxErr)
         If (MaxErr.gt.ThrCons) Then
            iStop=1
            Conv1=.False.
            Stop=.False.
            Write(Lu,'(A,ES11.4)') 'Maximum constraint error: ',MaxErr
            Write(Lu,*)
         End If
      End If
*
      nConst=0
      If (iNeg(1).eq.0) Then
         If (EDiffZero) Then
            If (NADC) Then
               nConst=2
            Else
               nConst=1
            End If
            Point_Desc='Minimum Energy Crossing Point Structure'
         Else
            Point_Desc='Minimum Structure'
         End If
      Else If (iNeg(1).eq.1) Then
         Point_Desc='Transition State Structure'
      Else
         Point_Desc='Higher Order Saddle Point Structure'
      End If
      If (nLambda.gt.nConst) Point_Desc='Constrained '//Trim(Point_Desc)
      If (iPrint.ge.5) Then
         If (Stop) Then
            If (Conv1) Then
               Write (Lu,'(A,I3,A)') ' Geometry is converged in ',
     &            kIter-iOff_iter,' iterations to a '//Trim(Point_Desc)
            Else
               Write (Lu,'(A)') ' No convergence after max iterations'
               If (Lu.ne.6) Write (6,'(/A)')
     &                          ' No convergence after max iterations'
            End If
         Else
            Write (Lu,'(A)') ' Convergence not reached yet!'
         End If
      End If
      If (FindTS.and.Stop.and.Conv1) Then
         If (.Not.TSReg) Then
            If (iPrint.ge.5) Then
               Write (Lu,*)
               Write (Lu,'(A)')
     &' FindTS was requested, but the TS regime was not reached.'
               Write (Lu,'(A)')
     &' The converged structure is probably not the desired TS.'
            End If
            iStop=16
         End If
      End If
c      If (iPrint.eq.7) Then
c         Write (Lu,*)
c         Write (Lu,'(A)') '*********************************'//
c     &      ' Geometry Statistics for Geometry Optimization '//
c     &                    '*********************************'
c         nPrint(118) = 7
c         Call List(' Internal coordinates ',Lbl,qInt,nInter,iter+1)
c         Call List(' Internal forces    ',Lbl,dqInt,nInter,iter)
c      End If
*
*     The energy change should not be too large
      Maxed=1.0d2
      If (Abs(E_Delta).gt.Maxed) Then
         Write (6,*) 'The predicted energy change is too large: ',
     &                E_Delta
         Write (6,'(A)') ' This can''t be right!'
         Write (6,'(A)') ' This job will be terminated.'
         iStop=8
         Stop=.True.
      End If
      If (iPrint.ge.5) Then
         Write (Lu,*)
         Write (Lu,'(A)') '*********************'//
     &      '*******************************************************'
     &      //'*************************************'
         Write (Lu,'(A)') '*********************'//
     &      '*******************************************************'
     &      //'*************************************'
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Terminate=.False.
      IRCRestart=.False.
*                                                                      *
************************************************************************
*                                                                      *
*-----Write summary of conical intersection characterization data
*
      If (Conv1.and.NADC.and.EDiffZero) Then
         If (iPrint.ge.5) Then
            If (.Not.ApproxNADC) Call CI_Summary(Lu)
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Book keeping for Saddle optimization for a TS. Note that this will
*     have to be done on both of the runfiles!
*
      If (Conv1.and.Saddle) Then
*
*        Here if a macro iteration in the Saddle TS optimization is
*        completed.
*
C        ENew=Energy(iter)+E_Delta
         Call mma_allocate(Tmp,nSaddle,Label='Tmp')
*
*        Store the info for later generation of MOLDEN formated files
*
         Call Get_dArray('Saddle',Tmp,nSaddle)
         E_Reac=Tmp(6*nAtom+1)
         E_Prod=Tmp(6*nAtom+2)
*
         Call mma_allocate(E_S,nSaddle_Max,Label='E_S')
         Call mma_allocate(C_S,3,nAtom,nSaddle_Max,Label='C_S')
         Call mma_allocate(G_S,3,nAtom,nSaddle_Max,Label='G_S')
         If (iSaddle.eq.0) Then
*
*           Initiate with data from the starting points
*
            E_S(:)=Zero
            C_S(:,:,:)=Zero
            G_S(:,:,:)=Zero
            iSaddle=1
            If (E_Reac.le.E_Prod) Then
               E_S(iSaddle)=E_Reac
               Call DCopy_(3*nAtom,Tmp(1:3*nAtom),1,C_S(:,:,iSaddle),1)
            Else
               E_S(iSaddle)=E_Prod
               Call DCopy_(3*nAtom,Tmp(3*nAtom+1:6*nAtom),1,
     &                     C_S(:,:,iSaddle),1)
            End If
*
         Else
*
            Call Get_dArray('MEP-Energies',E_S,nSaddle_Max)
            Call Get_dArray('MEP-Coor',C_S,3*nAtom*nSaddle_Max)
            Call Get_dArray('MEP-Grad',G_S,3*nAtom*nSaddle_Max)
*
         End If
*
*        Add the new data
*
         iSaddle=iSaddle+1
         E_S(iSaddle)=Energy(iter)
         C_S(:,:,iSaddle) = Cx(:,:,iter)
         G_S(:,:,iSaddle) = Gx(:,:,iter)
*
*        Put data on RUNFILE
*
         Call Put_dArray('MEP-Energies',E_S,nSaddle_Max)
         Call Put_dArray('MEP-Coor',C_S,3*nAtom*nSaddle_Max)
         Call Put_dArray('MEP-Grad',G_S,3*nAtom*nSaddle_Max)
         Call Put_iScalar('nMEP',iSaddle)
*
         Call mma_deallocate(G_S)
         Call mma_deallocate(C_S)
         Call mma_deallocate(E_S)
*                                                                      *
************************************************************************
*                                                                      *
*        Now update the "Saddle" field on both runfiles.
*
         Do iFile = 1, 2
*
         If (iFile.eq.1) Then
            Call NameRun('RUNREAC')
         Else
            Call NameRun('RUNPROD')
         End If
*
*        Update info on the runfile.
*
         Call Get_dArray('Saddle',Tmp,nSaddle)
         E1=Tmp(6*nAtom+1)
         E2=Tmp(6*nAtom+2)
C        Write (6,*) 'ENew=',ENew
C        Write (6,*) 'E1,E2=',E1,E2
         If (E1.le.E2) Then
C           Write (6,*) 'Update reactant'
            Tmp(6*nAtom+1)=Energy(iter)
            E1=Energy(iter)
            Call DCopy_(3*nAtom,Cx(:,:,iter),1,Tmp(1:3*nAtom),1)
         Else
C           Write (6,*) 'Update product'
            Tmp(6*nAtom+2)=Energy(iter)
            E2=Energy(iter)
            Call DCopy_(3*nAtom,Cx(:,:,iter),1,Tmp(3*nAtom+1:6*nAtom),1)
         End If
*        Set flag that seward should process the info! This should not
*        be done for the final macro iteration.
         If (.Not.FindTS) Tmp(6*nAtom+5)=One
         Call Put_dArray('Saddle',Tmp,nSaddle)
*
         End Do
*
*                                                                      *
************************************************************************
*                                                                      *
*        Converged or not, create the saddle.molden file
*        after each macro iteration
*
         Call NameRun('RUNREAC')
         Call mma_allocate(E_r,nSaddle_Max,Label='E_r')
         Call mma_allocate(C_r,3*nAtom,nSaddle_Max,Label='C_r')
         Call mma_allocate(G_r,3*nAtom,nSaddle_Max,Label='G_r')
         Call Qpg_iScalar('nMEP',Found)
         if(Found) Then
            Call Get_dArray('MEP-Energies',E_r,nSaddle_Max)
            Call Get_dArray('MEP-Coor',C_r,3*nAtom*nSaddle_Max)
            Call Get_dArray('MEP-Grad',G_r,3*nAtom*nSaddle_Max)
            Call Get_iScalar('nMEP',iSaddle_r)
         Else
            E_r(1)=E_Reac
            C_r(:,1) = Tmp(1:3*nAtom)
            G_r(:,:) = Zero
            iSaddle_r=1
         End If
*
         Call NameRun('RUNPROD')
         Call mma_allocate(E_p,nSaddle_Max,Label='E_p')
         Call mma_allocate(C_p,3*nAtom,nSaddle_Max,Label='C_p')
         Call mma_allocate(G_p,3*nAtom,nSaddle_Max,Label='G_p')
         Call Qpg_iScalar('nMEP',Found)
         if(Found) Then
            Call Get_dArray('MEP-Energies',E_p,nSaddle_Max)
            Call Get_dArray('MEP-Coor',C_p,3*nAtom*nSaddle_Max)
            Call Get_dArray('MEP-Grad',G_p,3*nAtom*nSaddle_Max)
            Call Get_iScalar('nMEP',iSaddle_p)
         Else
            E_p(1)=E_Prod
            C_p(:,1) = Tmp(3*nAtom+1:6*nAtom)
            G_p(:,:) = Zero
            iSaddle_p=1
         End If
*
*        Merge the two lists
*
         jSaddle=iSaddle_r
         Do iSaddle=iSaddle_p, 1, -1
            jSaddle=jSaddle+1
            E_r(jSaddle) = E_p(iSaddle)
            C_r(:,jSaddle) = C_p(:,iSaddle)
            G_r(:,jSaddle) = G_p(:,iSaddle)
         End Do
         Call mma_deallocate(G_p)
         Call mma_deallocate(C_p)
         Call mma_deallocate(E_p)
*
*        Align the structures sequentially, only for visualization
*        (gradients are not changed, though)
* TODO   Rotate the gradients too
*
         Do iSaddle=1,iSaddle_r+iSaddle_p-1
            Call Align(C_r(:,iSaddle+1),C_r(:,iSaddle),nAtom)
         End Do
*
         Call Intergeo('MD_SADDLE',E_r,C_r,
     &                 G_r,nAtom,iSaddle_r+iSaddle_p)
*
         Call mma_deallocate(G_r)
         Call mma_deallocate(C_r)
         Call mma_deallocate(E_r)
*
*                                                                      *
************************************************************************
*                                                                      *
*        If the Saddle TS optimization is not yet completed set up
*        data for the next macro iteration.
*
         If (.Not.FindTS) Then
            Call mma_deallocate(Tmp)
C           Write (*,*) 'Reset $SubProject'
*
*           Reset $SubProject for the next macro iteration.
*
            LuInput=11
            LuInput=IsFreeUnit(LuInput)
            Call StdIn_Name(StdIn)
            Call Molcas_Open(LuInput,StdIn)
            If (E1.le.E2) Then
               Write (LuInput,'(A)') '> EXPORT SubProject=.Reac'
C              Write (6,*) 'SubProject=.Reac'
               Call NameRun('RUNREAC')
            Else
               Write (LuInput,'(A)') '> EXPORT SubProject=.Prod'
C              Write (6,*) 'SubProject=.Prod'
               Call NameRun('RUNPROD')
            End If
*
*           Signal whether next iteration will be the first in the branch
*
            Call Qpg_iScalar('nMEP',Found)
            If (.Not.Found) Then
               Write (LuInput,'(A)') '> EXPORT SADDLE_FIRST=1'
            End If
            Write (LuInput,'(A,I3)') '> EXIT ',_RC_CONTINUE_LOOP_
            Close(LuInput)
*
*           Set flags to request yet another macro iteration.
*
            Terminate=.False.
            iStop=6
            Stop=.False.
         Else
            Call NameRun('RUNFILE')
            Call mma_deallocate(Tmp)
            nSaddle=0
            Call Put_dArray('Saddle',[Zero],nSaddle)
            Call Put_iScalar('nMEP',nSaddle)
         End If
*
*        Update the active runfile wrt the total number of micro
*        iterations done in all macro iterations of this branch.
*
         Call NameRun('RUNFILE')
         Call Put_iScalar('iOff_Iter',iter)
      End If
*
*     Disable first iteration signal right after the first iteration
*     (in each branch)
*
      If ((.Not.Conv1).and.Saddle) Then
         If ((iter.Eq.1).and.(iStop.eq.1)) Then
            LuInput=11
            LuInput=IsFreeUnit(LuInput)
            Call StdIn_Name(StdIn)
            Call Molcas_Open(LuInput,StdIn)
            Write (LuInput,'(A)') '> EXPORT SADDLE_FIRST=0'
            Write (LuInput,'(A,I3)') '> EXIT ',_RC_CONTINUE_LOOP_
            Close(LuInput)
            iStop=6
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Book keeping for minimum energy path search
*
      Call Qpg_iScalar('IRC',Found)
      If (Found) Call Get_iScalar('IRC',IRC)
*
      TurnBack=.False.
      If (MEP.or.rMEP) Then
       If (Conv1) Then
*
*        Is this the first iteration or not?
*
         iMEP=iMEP+1
*
*        Save information for the current step
*
         Call mma_allocate(E_MEP,nMEP+1,Label='E_MEP')
         Call mma_allocate(C_MEP,3,nAtom,nMEP+1,Label='C_MEP')
         Call mma_allocate(G_MEP,3,nAtom,nMEP+1,Label='G_MEP')
         If (iMEP.gt.1) Then
            Call Get_dArray('MEP-Energies',E_MEP,nMEP+1)
            Call Get_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
            Call Get_dArray('MEP-Grad',G_MEP,3*nAtom*(nMEP+1))
         Else
            E_MEP(:)=Zero
            C_MEP(:,:,:)=Zero
            G_MEP(:,:,:)=Zero
            E_MEP(iMEP)=Energy(iOff_iter+1)
            C_MEP(:,:,iMEP)= Cx(:,:,iOff_iter+1)
            G_MEP(:,:,iMEP)= Gx(:,:,iOff_iter+1)
         End If
*
         E_MEP(iMEP+1)=Energy(iter)
         C_MEP(:,:,iMEP+1)= Cx(:,:,iter)
         G_MEP(:,:,iMEP+1)= Gx(:,:,iter)
         Call Put_dArray('MEP-Energies',E_MEP,nMEP+1)
         Call Put_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
         Call Put_dArray('MEP-Grad',G_MEP,3*nAtom*(nMEP+1))
         Call Put_iScalar('nMEP',iMEP)
*
*        Save the path so far (energies, coordinates and forces)
*
         Call Intergeo('MD_MEP',E_MEP,C_MEP,G_MEP,nAtom,iMEP+1)
*
*        Compute energy difference and RMS between last two structures
*
         eDiffMEP=E_MEP(iMEP+1)-E_MEP(iMEP)
         Call mma_allocate(Coor1,3,mTtAtm,Label='Coor1')
         Call mma_allocate(Coor2,3,mTtAtm,Label='Coor2')
         Call AtmLst(C_MEP(:,:,iMEP  ),nAtom,Coor1,mTtAtm)
         Call AtmLst(C_MEP(:,:,iMEP+1),nAtom,Coor2,mTtAtm)
         Call OptRMS_Slapaf(Coor1,Coor2,mTtAtm,RMS,RMSMax)
         Call mma_deallocate(Coor1)
         Call mma_deallocate(Coor2)
*
         Call mma_deallocate(G_MEP)
         Call mma_deallocate(C_MEP)
         Call mma_deallocate(E_MEP)
*
       Else
*
*        Test for "turn back", i.e. when the trial structure
*        is getting too close to the previous converged structure,
*        this may be an indication of an ill-behaved constraint
         If (iMEP.ge.1) Then
            Call mma_allocate(C_MEP,3,nAtom,nMEP+1,Label='C_MEP')
            Call mma_allocate(OfRef,3,nAtom,Label='OfRef')
            Call mma_Allocate(Tmp,3*nAtom,Label='Tmp')
            Call Get_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
            OfRef(:,:)=C_MEP(:,:,iMEP+1)

*           Using hypersphere measure, even with "transverse" MEPs,
*           this should not be a problem
            Call SphInt(Cx(:,:,iter),nAtom,Not_Allocated,refDist,Tmp,
     &         .False.,'dummy   ',rDum,.False.)
            Call SphInt(Cx(:,:,iter),nAtom,OfRef,prevDist,
     &                  Tmp,.False.,'dummy   ',rDUm,.False.)
            If (prevDist.lt.Half*refDist) Then
               TurnBack=.True.
               Conv1=.True.
               Stop=.True.
               iStop=0
               Terminate=.True.
            End If
            Call mma_deallocate(Tmp)
            Call mma_deallocate(OfRef)
            Call mma_deallocate(C_MEP)
         End If
*
       End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----List internal coordinates and gradients
*
      kkIter=iter+1
      If (iPrint.ge.8) Then
         Call List(' Internal coordinates ',Lbl,qInt,nInter,kkIter)
         Call List(' Internal forces    ',Lbl,dqInt,nInter,iter)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Put out the new reference structure and the new starting
*     structure to be used for the next MEP point.
*     For rMEP keep the reference structure!
*     Note that this is done in weighted Cartesian coordinates!
*
      If ((Conv1.or.(iter.eq.1)).and.(MEP.or.rMEP)) Then
         If ((iMEP.ge.1).and.(iPrint.ge.5)) Then
            Write (6,*)
            Call CollapseOutput(1,'IRC/Minimum Energy Path Information')
         End If
*
         ResGrad=Huge(ResGrad)
         If (.Not.Terminate) Then
            BadConstraint=.False.
            Call MEP_Dir(Cx,Gx,nAtom,iMEP,iOff_iter,iPrint,IRCRestart,
     &                   ResGrad,BadConstraint)
            Call dCopy_(3*nAtom,Cx(:,:,iter+1),1,Coor,1)
            Call Put_iScalar('iOff_Iter',iter)
         End If
*
         If (MEP) Then
            If (IRC.eq.0) Then
               MEP_Text='MEP'
            Else If (IRC.eq.1) Then
               MEP_Text='IRC(forward)'
            Else
               MEP_Text='IRC(backward)'
            End If
         Else If (rMEP) Then
            MEP_Text='rMEP'
         Else
            MEP_Text=''
         End If
*
*        Should we terminate or not? Not done on the first iteration.
*
         If (iMEP.gt.0) Then
*
*           Test for energy increase (optionally disabled).
            eTest=eMEPTest.and.(eDiffMEP.gt.Zero)
            If ((MEP.and.eTest).and.(.not.Terminate)) Then
               Terminate=.True.
               If (iPrint.ge.5) Then
                  Write (6,*)
                  Write (6,'(A)') ' '//Trim(MEP_Text)//'-search'//
     &               ' terminated due to energy increase!'
                  Write (6,*)
               End If
            End If
*
*           Test for energy decrease (optionally disabled).
            eTest=eMEPTest.and.(eDiffMEP.lt.Zero)
            If ((rMEP.and.eTest).and.(.not.Terminate)) Then
               Terminate=.True.
               If (iPrint.ge.5) Then
                  Write (6,*)
                  Write (6,'(A)') ' '//Trim(MEP_Text)//'-search'//
     &               ' terminated due to energy decrease!'
                  Write (6,*)
               End If
            End If
*
*           Test for small gradient.
            If ((iMEP.gt.1).or.(IRC.eq.0)) Then
              If ((ResGrad.lt.ThrMEP).and.(.not.Terminate)) Then
                 Terminate=.True.
                 If (iPrint.ge.5) Then
                    Write (6,*)
                    Write (6,'(A)') ' '//Trim(MEP_Text)//'-search'//
     &               ' terminated due to small gradient!'
                    Write (6,*)
                 End If
              End If
            End If
*
*           Test for small step.
            If ((RMS.lt.ThrGrd*4.D0).and.(.not.Terminate)) Then
               Terminate=.True.
               If (iPrint.ge.5) Then
                  Write (6,*)
                  Write (6,'(A)') ' '//Trim(MEP_Text)//'-search'//
     &               ' terminated due to small geometry change!'
                  Write (6,*)
               End If
            End If
*
*           Test for max number of points.
            If ((iMEP.ge.nMEP).and.(.not.Terminate)) Then
               Terminate=.True.
               If (iPrint.ge.5) Then
                  Write (6,*)
                  Write (6,'(A)') ' '//Trim(MEP_Text)//'-search'//
     &               ' terminated due to max number of path points!'
                  Write (6,*)
               End If
            End If
*
*           Test for constraint misbehavior.
            If ((BadConstraint.and.(.not.Terminate)).or.TurnBack) Then
               Terminate=.True.
               If (iPrint.ge.5) Then
                  Write (6,*)
                  Write (6,'(A)') ' '//Trim(MEP_Text)//'-search'//
     &               ' terminated due to problematic constraint!'
                  Write (6,*)
               End If
            End If
*
*           If IRC reset for backward IRC search.
*
            If (Terminate) Then
               If (IRC.eq.1) Then
                  IRCRestart=.True.
               End If
            End If
*
         End If
*
         If (Conv1) Then
            Call Chkpnt_update_MEP(.Not.TurnBack,IRCRestart)
         End If
*
         If (Conv1.and.Terminate) Then
            If (IRC.ne.0) Then
               Call mma_allocate(E_MEP,nMEP+1,Label='E_MEP')
               Call mma_allocate(C_MEP,3,nAtom,nMEP+1,Label='C_MEP')
               Call mma_allocate(G_MEP,3,nAtom,nMEP+1,Label='G_MEP')
               Call Get_dArray('MEP-Energies',E_MEP,nMEP+1)
               Call Get_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
               Call Get_dArray('MEP-Grad',G_MEP,3*nAtom*(nMEP+1))
               If (IRC.eq.1) Then
                  IRCRestart=.True.
                  IRC=-1
                  Call Put_iScalar('IRC',IRC)
*
*                 Store away data for IRC molden file. Forward part.
*
                  Call Put_dArray('IRC-Energies',E_MEP,iMEP+1)
                  Call Put_dArray('IRC-Coor',C_MEP,3*nAtom*(iMEP+1))
                  Call Put_dArray('IRC-Grad',G_MEP,3*nAtom*(iMEP+1))
                  Call Put_dArray('Ref_Geom',Cx,3*nAtom)
*
*                 Write a temporary file
*                 (will be overwritten when the backward part is done)
*
                  Call Intergeo('MD_IRC',E_MEP,C_MEP,G_MEP,nAtom,iMEP+1)
*
                  Terminate=.False.
*
               Else If (IRC.eq.-1) Then
*
*                 Assemble molden file for IRC
*
                  nBackward=iMEP+1
                  Call qpg_dArray('IRC-Energies',Found,nForward)
                  nIRC=nForward+nBackward-1
                  Call mma_allocate(E_IRC,nIRC,Label='E_IRC')
                  Call mma_allocate(C_IRC,3,nAtom,nIRC,Label='C_IRC')
                  Call mma_allocate(G_IRC,3,nAtom,nIRC,Label='G_IRC')
*
                  j=0
                  Do i = nBackward, 1, -1
                     j = j+1
                     E_IRC(j)=E_MEP(i)
                     C_IRC(:,:,j) = C_MEP(:,:,i)
                     G_IRC(:,:,j) = G_MEP(:,:,i)
                  End Do
*
                  Call Get_dArray('IRC-Energies',E_IRC(nBackward),
     &                            nForward)
                  Call Get_dArray('IRC-Coor',C_IRC(:,:,nBackward:nIRC),
     &                            nForward*3*nAtom)
                  Call Get_dArray('IRC-Grad',G_IRC(:,:,nBackward:nIRC),
     &                            nForward*3*nAtom)
*
                  Call Intergeo('MD_IRC',E_IRC,C_IRC,G_IRC,nAtom,nIRC)
*
                  Call mma_deallocate(G_IRC)
                  Call mma_deallocate(C_IRC)
                  Call mma_deallocate(E_IRC)
               End If

               Call mma_deallocate(G_MEP)
               Call mma_deallocate(C_MEP)
               Call mma_deallocate(E_MEP)
            End If
         End If
*
         If (.Not.Terminate) Then
           iStop=1
           Stop=.False.
         End If
*
*        Print out the path so far
*
         If ((iMEP.ge.1).and.(iPrint.ge.5)) Then
            Call mma_allocate(E_MEP,nMEP+1,Label='E_MEP')
            Call mma_allocate(C_MEP,3,nAtom,nMEP+1,Label='C_MEP')
            Call mma_allocate(L_MEP,nMEP+1,Label='L_MEP')
            Call mma_allocate(Cu_MEP,nMEP+1,Label='Cu_MEP')
            Call Get_dArray('MEP-Energies',E_MEP,nMEP+1)
            Call Get_dArray('MEP-Coor',C_MEP,3*nAtom*(nMEP+1))
            Call Get_dArray('MEP-Lengths',L_MEP,nMEP+1)
            Call Get_dArray('MEP-Curvatures',Cu_MEP,nMEP+1)
            Write (6,*)
            CumLen=Zero
            If (Cu_MEP(1+iMEP).ge.Zero) Then
               Write(6,*) '         Cumul.'
               Write(6,*) 'Point  Length (bohr)       Energy  Curvature'
               Write(6,*) '--------------------------------------------'
               Do i = 0, iMEP
                  CumLen=CumLen+L_MEP(1+i)
                  Write (6,200) i,CumLen,E_MEP(1+i),Cu_MEP(1+i)
               End Do
            Else
               Write(6,*) '         Cumul.'
               Write(6,*) 'Point  Length (bohr)       Energy'
               Write(6,*) '---------------------------------'
               Do i = 0, iMEP
                  CumLen=CumLen+L_MEP(1+i)
                  Write (6,200) i,CumLen,E_MEP(1+i)
               End Do
            End If
200         Format (1X,I5,1X,F10.6,1X,F16.8,1X,F10.6)
            If (iPrint.gt.6) Then
               Write (6,*)
               Do i = 0, iMEP
                  Call RecPrt(' Coordinates',' ',C_MEP(:,:,i+1),3,nAtom)
               End Do
            End If
            Call CollapseOutput(0,'IRC/Minimum Energy Path Information')
            Write(6,*)
            Call mma_deallocate(E_MEP)
            Call mma_deallocate(C_MEP)
            Call mma_deallocate(L_MEP)
            Call mma_deallocate(Cu_MEP)
         End If
*
      End If
*
      If (IRCRestart) Then
*
*         Prepare the runfile to start from Scratch
*
          iMEP=0
          Call Put_iScalar('nMEP',iMEP)
          Call mma_allocate(Information,7,Label='Information')
          Information(:)=0
          Information(1)=-99     ! Deactivate the record
          Call Put_iArray('Slapaf Info 1',Information,7)
          Call mma_deallocate(Information)
          iOff_Iter=0
          Call Put_iScalar('iOff_Iter',iOff_Iter)
*
*         Restore data
*
          Call Put_dScalar('Last Energy',Energy(1))
          Call Put_dArray('GRAD',Gx,3*nAtom)
          Call Put_dArray('Unique Coordinates',Cx,3*nAtom)
          Call Put_Coord_New(Cx,nAtom)
          call dcopy_(3*nAtom,Cx,1,Coor,1)
          call dcopy_(3*nAtom,Cx,1,Cx(:,:,iter+1),1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Figure out if the last energy should be computed!
*
      Last_Energy = Stop .and. iStop.ne.16 .and. iStop.ne.8
      Last_Energy = Last_Energy .and. .Not.MEP .and. .Not.rMEP
      Last_Energy = Last_Energy .and.
     &             .Not. (Numerical .and. kIter.eq.1)
      If (Last_Energy) iStop = 2
*                                                                      *
************************************************************************
*                                                                      *
*     Write geometry file for MOLDEN. Note that this file should not be
*     generated if Slapaf is running new geometries for any numerical
*     procedure!
*
      If (
     &     .NOT. Just_Frequencies  .AND.
     &     .NOT. Numerical
     &   ) Then
         Call Write_QMMM(Cx,nAtom,iter)
         Call Intergeo('MD_GEO',Energy,Cx,Gx,nAtom,iter+1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
