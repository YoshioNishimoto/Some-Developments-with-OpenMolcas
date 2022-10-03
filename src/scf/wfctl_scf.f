************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it ana/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               1995,1996, Martin Schuetz                              *
*               2003, Valera Veryazov                                  *
*               2016,2017, Roland Lindh                                *
************************************************************************
      SubRoutine WfCtl_SCF(iTerm,Meth,FstItr,SIntTh)
      use SCF_Arrays
      use InfSCF
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
#include "mxdm.fh"
      Logical FstItr
      Character(LEN=*) Meth
*
      nD = iUHF + 1
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate memory for some arrays
*
      nTr=MxIter
      Call mma_allocate(TrDh,nTR**2,nD,Label='TrDh')
      Call mma_allocate(TrDP,nTR**2,nD,Label='TrDP')
      Call mma_allocate(TrDD,nTR**2,nD,Label='TrDD')
      nCI = MxOptm + 1
      Call mma_allocate(CInter,nCI,nD,Label='CInter')
*                                                                      *
************************************************************************
*                                                                      *
      Call WfCtl_SCF_Internal(iTerm,Meth,FstItr,SIntTh,OneHam,TwoHam,
     &                        Dens,Ovrlp,TrDh,TrDP,TrDD,
     &                        CInter,OccNo,Vxc,TrM,nBT,
     &                        nDens,nD,nTr,nBB,nCI,nnB
     &                       )
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(CInter)
      Call mma_deallocate(TrDD)
      Call mma_deallocate(TrDP)
      Call mma_deallocate(TrDh)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine WFCtl_SCF
      SubRoutine WfCtl_SCF_Internal(
     &                      iTerm,Meth,FstItr,SIntTh,
     &                      OneHam,TwoHam,Dens,Ovrlp,
     &                      TrDh,TrDP,TrDD,CInter,OccNo,
     &                      Vxc,TrM,mBT,mDens,nD,nTr,mBB,nCI,mmB
     &                             )
************************************************************************
*                                                                      *
*     purpose: Optimize SCF wavefunction.                              *
*                                                                      *
*     called from: SCF                                                 *
*                                                                      *
*     calls to: PrBeg,Aufbau,DMat,PMat,EneClc,SOIniH,UpdFck,           *
*               TraFck,DIIS_x,DIIS_i,NewOrb,MODens,PrIte               *
*               uses SubRoutines and Functions from Module cycbuf.f    *
*               -cyclic buffer implementation                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*     QNR/DIIS, M. Schuetz, 1995                                       *
*     NDDO start orbitals, M. Schuetz, 1996                            *
*     University of Lund, Sweden, 1992,95,96                           *
*     UHF, V.Veryazov, 2003                                            *
*     Cleanup, R. Lindh, 2016                                          *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
#ifdef _MSYM_
      Use, Intrinsic :: iso_c_binding, only: c_ptr
#endif
      Use Interfaces_SCF, Only: TraClc_i
      use LnkLst, only: SCF_V
      use LnkLst, only: LLGrad,LLDelt,LLx
      use InfSO
      use SCF_Arrays, only: EOrb, CMO, Fock
      use InfSCF
      Implicit Real*8 (a-h,o-z)
      External Seconds
      Real*8 Seconds
      Real*8 OneHam(mBT), TwoHam(mBT,nD,mDens), Dens(mBT,nD,mDens),
     &       Ovrlp(mBT), TrDD(nTr*nTr,nD),
     &       TrDh(nTr*nTr,nD), TrDP(nTr*nTr,nD),
     &       CInter(nCI,nD), Vxc(mBT,nD,mDens), TrM(mBB,nD),
     &       OccNo(mmB,nD)
#include "real.fh"
#include "stdalloc.fh"
#include "file.fh"
#include "twoswi.fh"
#include "ldfscf.fh"
#include "warnings.h"
#include "mxdm.fh"
      Real*8, Dimension(:),   Allocatable:: D1Sao
      Real*8, Dimension(:), Allocatable:: Grd1, Disp, Xnp1, Xn

*---  Tolerance for negative two-electron energy
      Real*8 E2VTolerance
      Parameter (E2VTolerance=-1.0d-8)

*---  Define local variables
      Real*8 ::  StepMax=0.100D0
      Logical QNR1st,FstItr
      Logical :: Converged=.False.
      Character Meth*(*), Meth_*10
      Character*72 Note
      Logical AufBau_Done, Diis_Save, Reset, Reset_Thresh, AllowFlip
      Logical ScramNeworb
      Integer iAufOK, Ind(MxOptm)
      Character*128 OrbName
#ifdef _MSYM_
      Type(c_ptr) msym_ctx
#endif
      Dimension Dummy(1),iDummy(7,8)
      External DNRM2_
      Integer iSym
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call CWTime(TCpu1,TWall1)
*---  Put empty spin density on the runfile.
      Call mma_allocate(D1Sao,nBT,Label='D1Sao')
      Call FZero(D1Sao,nBT)
      Call Put_D1Sao(D1Sao,nBT)
      Call mma_deallocate(D1Sao)
*
      If (Len_Trim(Meth) > Len(Meth_)) Then
        Meth_ = '[...]'
      Else
        Meth_ = Trim(Meth)
      End If
      iTerm=0
      iDMin = - 1
*---  Choose between normal and minimized differences
      MinDMx = 0
      If(MiniDn) MinDMx = Max(0,nIter(nIterP)-1)
*
*
*     Optimization options
*
*     iOpt=0: DIIS interpolation on density or/and Fock matrices
*     iOpt=1: DIIS extrapolation on gradients w.r.t orbital rot.
*     iOpt=2: DIIS extrapolation on the anti-symmetrix X matrix.
*     iOpt=3: RS-RFO in the space of the anti-symmetric X matrix.
*
      iOpt=0
      QNR1st=.TRUE.
*
*     START INITIALIZATION
*
*---
*---  Initialize Cholesky information if requested
*
      If (DoCholesky) Then
         If (DoLDF) Then ! Local DF
            Call LDF_X_Init(.True.,0,LDFracMem,irc)
            If (irc.ne.0) Then
               Call WarningMessage(2,
     &                    'WfCtl: non-zero return code from LDF_X_Init')
               Call LDF_Quit(1)
            End If
            Call LDF_SetIntegralPrescreeningInfo()
            Call Free_iSD()
         Else ! Cholesky or RI/DF
            Call Cho_X_init(irc,ChFracMem)
            If (irc.ne.0) Then
               Call WarningMessage(2,
     &                              'WfCtl. Non-zero rc in Cho_X_init.')
               Call Abend
            End If
         End If
      End If
*                                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
      iterSt=iter0
      iterso=0        ! number of second order steps.
      kOptim=1
      Iter_no_Diis=2
      Converged=.False.
*                                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*     END INITIALIZATION
*                                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
      DiisTh=Max(DiisTh,QNRTh)
*
*     If DIIS is turned off make threhold for activation impossible
*
      If (.NOT. DIIS) Then
         DiisTh=Zero
         Iter_no_Diis=1000000
      End If
*
*     If no damping to preceed the DIIS make threshold being fullfiled
*     at once.
*
      If (.NOT. Damping) Then
         DiisTh=DiisTh*1.0D99
         Iter_no_Diis=1
      End If
*
*---  turn temporarily off DIIS & QNR/DIIS, if Aufbau is active...
*
      DiisTh_save=DiisTh
      Diis_Save=Diis
      If(Aufb) Then
         DiisTh=Zero
         Diis=.False.
      End If
*
*---  Print header to iterations
*
      If(KSDFT.eq.'SCF'.or.One_Grid) Call PrBeg(Meth_)
      AufBau_Done=.False.
*                                                                      *
*======================================================================*
*======================================================================*
*                                                                      *
*---  If direct SCF modify thresholds temporarily
*
      Reset=.False.
      Reset_Thresh=.False.
      EThr_New=EThr*100.0D0
*
*--- pow: temporary disabling of threshold switching
*
      If((DSCF.or.KSDFT.ne.'SCF').and.nIter(nIterP).gt.10) Then
*
         If(DSCF.and.KSDFT.eq.'SCF'.and.Two_Thresholds) Then
            Reset=.True.
            Reset_Thresh=.True.
            Call Reduce_Thresholds(EThr_New,SIntTh)
         End If
*
         If(KSDFT.ne.'SCF'.and..Not.One_Grid) Then
            Reset=.True.
            Call Modify_NQ_Grid()
            Call PrBeg(Meth_)
         End If
*
      End If
*
*======================================================================*
*======================================================================*
*                                                                      *
*                     I  t  e  r  a  t  i  o  n  s                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* nIterP=0: NDDO                                                       *
* nIterP=1: SCF                                                        *
*                                                                      *
*======================================================================*
*
*     Set some paramters to starting defaults.
*
      AllowFlip=.true.
      iAufOK=0
      IterX=0
      If(Scrmbl) IterX=-2
      Iter_DIIS=0
      EDiff=0.0D0
      DMOMax=0.0D0
      FMOMax=0.0D0
      DltNrm=0.0D0

      If(MSYMON) Then
#ifdef _MSYM_
         Write(6,*) 'Symmetrizing orbitals during SCF calculation'
         Call fmsym_create_context(msym_ctx)
         Call fmsym_set_elements(msym_ctx)
         Call fmsym_find_symmetry(msym_ctx)
#else
         Write(6,*) 'No msym support, skipping symmetrization '
     $           //'of SCF orbitals...'
#endif
      End If
*                                                                      *
*======================================================================*
*                                                                      *
*     Start of iteration loop
*                                                                      *
*======================================================================*
*                                                                      *
      Do iter_ = iterSt+1, iterSt+nIter(nIterP)
         iter = iter_
         IterX=IterX+1
         WarnCfg=.false.
*
         If(.not.Aufb.and.iter.gt.MaxFlip) AllowFlip=.false.

         TCP1=seconds()

         iDMin = iDMin + 1
         If(iDMin.gt.MinDMx) iDMin = MinDMx

#ifdef _MSYM_
         If (MSYMON) Then
            Do iD = 1, nD
               Call fmsym_symmetrize_orbitals(msym_ctx,CMO(1,iD))
               Call ChkOrt(CMO(1,iD),nBO,Ovrlp,nBT,Whatever)
               Call Ortho(CMO(1,iD),nBO,Ovrlp,nBT)
               Call ChkOrt(CMO(1,iD),nBO,Ovrlp,nBT,Whatever)
            End Do
         End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Do Aufbau procedure, if active...
*
         If (Aufb) Call Aufbau(nAufb,OccNo,mmB,iAufOK,nD)
*                                                                      *
************************************************************************
*                                                                      *
*        Save the variational energy of the previous step
*
         EnVold=EneV
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Call SCF_Energy(FstItr,E1V,E2V,EneV)
         Energy(iter)=EneV
         If(iter.eq.1) Then
            EDiff  = Zero
         Else
            EDiff = EneV-EnVold
         End If
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!        Select on the fly the current optimization method
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
!        Test if this is a DIIS extrapolation iteration, alternatively
!        the the iteration is a DIIS interpolation iteration. The former
!        is activated if the DMOMax is lower than the threshold
!        after a specific number of iteration, or if the condition
!        has already been achived.
!
!        2017-02-03: add energy criterion to make sure that the DIIS
!                    gets some decent to work with.
!
         If ((DMOMax.lt.DiisTh .AND. IterX.gt.Iter_no_Diis
     &             .AND. EDiff.lt.1.0D-1)) Then
*
*           Reset kOptim such that the extraploation scheme is not
*           corrupted by iterations with too high energies. Those can
*           not be used in the scheme.
*
            If (iOpt.eq.0) kOptim=2
            iOpt=Max(1,iOpt)
            Iter_DIIS = Iter_DIIS + 1
         End If
*
*        Test if the DIIS scheme will be operating in an orbital
*        rotation mode or linear combination of density matrices. This
*        option is avaliable only in the extrapolation mode of DIIS.
*
*        2017-02-03: Make sure that the density based DIIS is in
*                    action for at least 2 iterations such that
*                    when the orbital rotation DIIS is turned on
*                    we are firmly in the NR region.
*
         If (iOpt.ge.2 .OR.
     &      (iOpt.eq.1 .AND. DMOMax.lt.QNRTh .AND.  Iter_DIIS.ge.2))
     &      Then
            If (RSRFO.or.RGEK) Then
               If (RSRFO) Then
                  iOpt=3
               Else
                  iOpt=4
               End If
               kOptim=2
            Else
               iOpt=2
            End If

         End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Select Case(iOpt)

         Case(0)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*           Interpolation DIIS
*
*           Compute traces: TrDh, TrDP and TrDD.
*
            Call TraClc_i(OneHam,Dens,TwoHam,Vxc,nBT,nDens,iter,
     &                    TrDh,TrDP,TrDD,MxIter,nD)
*
            If (kOptim.eq.1) Then
*
*              If we only have one density, then nothing much to intra-
*              polate over.
*
               AccCon = 'None     '
*
            Else
*
*              DIIS interpolation optimization: EDIIS, ADIIS, LDIIS
*
               iOpt_DIIS=1 ! EDIIS option
*
               Call DIIS_i(CInter,nCI,TrDh,TrDP,TrDD,MxIter,nD,
     &                     iOpt_DIIS,Ind)
*
*----          Compute optimal density, dft potentials, and TwoHam
*
               Call OptClc(Dens,TwoHam,Vxc,nBT,nDens,CInter,nCI,nD,Ind)
*
            End If
*
*---        Update Fock Matrix from OneHam and extrapolated TwoHam & Vxc
*
            Call UpdFck(OneHam,TwoHam,Vxc,nBT,nDens,Fock,nIter(nIterP),
     &                  nD)
*---        Diagonalize Fock matrix and obtain new orbitals
*
            ScramNeworb=Scrmbl.and.iter.eq.1
            Call NewOrb_SCF(Fock,nBT,CMO,nBO,FMOMax,EOrb,nnO,Ovrlp,nFO,
     &                  AllowFlip,ScramNeworb,nD)
*
*---        Transform density matrix to MO basis
*
            Call MODens(Dens,Ovrlp,nBT,nDens,CMO,nBB,nD)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Case(1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*           Extrapolation DIIS
*
*           The section of the code operates on the derivatives of
*           the energy w.r.t the elements of the anti-symmetric X
*           matrix.
*
*           The optimal density matrix is found with the DIIS and
*           canonical CMOs are generated by the diagonalization of the
*           Fock matrix.
*
            If (FrstDs) iter_ref=1
            Call GrdClc(FrstDs)
*
            Call DIIS_x(nD,CInter,nCI,iOpt.eq.2,Ind)
*
*----       Compute optimal density, dft potentials, and TwoHam
*
            Call OptClc(Dens,TwoHam,Vxc,nBT,nDens,CInter,nCI,nD,Ind)
*
*---        Update Fock Matrix from OneHam and extrapolated TwoHam & Vxc
*
            Call UpdFck(OneHam,TwoHam,Vxc,nBT,nDens,Fock,nIter(nIterP),
     &                  nD)
*---        Diagonalize Fock matrix and obtain new orbitals
*
            ScramNeworb=Scrmbl.and.iter.eq.1
            Call NewOrb_SCF(Fock,nBT,CMO,nBO,FMOMax,EOrb,nnO,Ovrlp,nFO,
     &                  AllowFlip,ScramNeworb,nD)
*
*---        Transform density matrix to MO basis
*
            Call MODens(Dens,Ovrlp,nBT,nDens,CMO,nBB,nD)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Case(2,3,4)
*                                                                      *
*           Extrapolation DIIS & QNR
*
*           or
*
*           Quasi-2nd order scheme (rational function optimization)
*           with restricted.step.
*
*           In this section we operate directly on the anti-symmetric X
*           matrix.
*
*           Initially the energy is determined through a line seach,
*           followed by the Fock matrix etc. In this respect, this
*           approach is doing in this iteration what was already done
*           by the two other approaches in the end of the previous
*           iteration.
*
*           Note also, that since we work directly with the X matrix we
*           also rotate the orbitals with this matrix (note the call to
*           RotMOs right before the call to SCF_Energy in the line
*           search routine, linser). Thus, the orbitals here are not
*           canonical and would require at the termination of the
*           optimization that these are computed.
*
*           Initiate if the first QNR step
*
            If (QNR1st) Then
               iter_Ref = iter
!------        1st QNR step, reset kOptim to 1

               kOptim = 1
               CInter(1,1) = One
               CInter(1,nD) = One

!              init 1st orb rot parameter X1 (set it to zero)
               Call mma_allocate(Xn,mOV,Label='Xn')
               Xn(:)=Zero
!              and store it on appropriate LList
               Call PutVec(Xn,mOV,iter,'NOOP',LLx)
               Call mma_deallocate(Xn)
*
*---           compute initial diagonal Hessian, Hdiag
*
               Call SOIniH()
               AccCon(8:8)='H'
            Else
               AccCon(8:8)=' '
            End If

            Call GrdClc(QNR1st)
*
            Call dGrd()
*
*---        Update the Fock Matrix from actual OneHam, Vxc & TwoHam
*           AO basis
*
            Call UpdFck(OneHam,TwoHam,Vxc,nBT,nDens,Fock,
     &                  nIter(nIterP),nD)
*
*---        and transform the Fock matrix into the new MO space,
*
            Call TraFck(Fock,nBT,CMO,nBO,.FALSE.,FMOMax,
     &                  EOrb,nnO,Ovrlp,nD)
*
*---        update QNR iteration counter
*
            iterso=iterso+1

*           Allocate memory for the current gradient and
*           displacement vector.
*
            Call mma_allocate(Grd1,mOV,Label='Grd1')
            Call mma_allocate(Disp,mOV,Label='Disp')
            Call mma_allocate(Xnp1,mOV,Label='Xnp1')

         Select Case(iOpt)

         Case(2)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*
*----       Compute extrapolated g_x(n) and X_x(n)
*
            Call DIIS_x(nD,CInter,nCI,iOpt.eq.2,Ind)
            Call OptClc_QNR(CInter,nCI,nD,Grd1,Xnp1,mOV,Ind,MxOptm,
     &                      kOptim,kOV)
*
*-------    compute new displacement vector delta
*           dX_x(n) = -H(-1)*g_x(n) ! Temporary storage in Disp
*
            Call SOrUpV(Grd1,mOV,Disp,'DISP','BFGS')
            Disp(:)=-Disp(:)
!
!           from this, compute new orb rot parameter X(n+1)
!
!           X(n+1) = X_x(n) - H(-1)g_x(n)
!           X(n+1) = X_x(n) + dX_x(n)
!
            Xnp1(:)=Xnp1(:)+Disp(:)
*
*           get address of actual X(n) in corresponding LList
*
            jpXn=LstPtr(iter,LLx)
*
*           and compute actual displacement dX(n)=X(n+1)-X(n)
*
            Disp(:)=Xnp1(:)-SCF_V(jpXn)%A(:)
            dqdq = DNRM2_(mOV,Disp,1)
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
         Case(3,4)

*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
*           Get g(n)
*
            Call GetVec(iter,LLGrad,inode,Grd1,mOV)
*
*           Use rs-rfo to compute dX(n)
*
*
            If (iter>1     .AND.
     &          EDiff<Zero .AND.
     &          Abs(EDiff)<1.0D-1
     &         ) Then
               Write (6,*) 'Increase step restriction parameter.'

*              Increase steplength if there was an energy decrease.

               StepMax=Min(StepMax*Two,0.45D0)

             Else If (iter>1 .AND.
     &                EDiff>Zero ) Then

*              If last step represented an energy increase reset the
*              threshold parameter for the step restriction
               Write (6,*) 'Decrease step restriction parameter.'

               StepMax=Max(StepMax*0.75D0,0.8D-3)
            End If

************************************************************************
*                                                                      *
            Select Case(iOpt)

            Case(3)

               dqHdq=Zero
               Call rs_rfo_scf(Grd1,mOV,Disp,AccCon(1:6),dqdq,
     &                         dqHdq,StepMax,AccCon(9:9))

            Case(4)
*                                                                      *
               Call DIIS_GEK_Optimizer(Disp,mOV,dqdq,AccCon(1:6),
     &                                               AccCon(9:9))
            End Select
*                                                                      *
************************************************************************

*           Pick up X(n) and compute X(n+1)=X(n)+dX(n)

            Call GetVec(iter,LLx,inode,Xnp1,mOV)

            Xnp1(:)=Xnp1(:)+Disp(:)

         End Select
*                                                                      *
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
*           Store X(n+1) and dX(n)
*
            Call PutVec(Xnp1,mOV,iter+1,'NOOP',LLx)
            Call PutVec(Disp,mOV,iter,'NOOP',LLDelt)
*
*           compute Norm of dX(n)
*
            DltNrm=DBLE(nD)*dqdq
*
*           Generate the CMOs, rotate MOs accordingly to new point
*
            Call RotMOs(Disp,mOV,CMO,nBO,nD,Ovrlp,mBT)
*
*           and release memory...
            Call mma_deallocate(Xnp1)
            Call mma_deallocate(Disp)
            Call mma_deallocate(Grd1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Case Default
            Write (6,*) 'WfCtl_SCF: Illegal option'
            Call Abend()
         End Select
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*----    Update DIIS interpolation depth kOptim
*
         If (iOpt<=3) Then
            If (idKeep.eq.0) Then
               kOptim = 1
            Else If (idKeep.eq.1) Then
               kOptim = 2
            Else
               kOptim = Min(kOptim_Max,kOptim + 1)
            End If
         End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*        End optimization section
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
         Call Save_Orbitals()
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*======================================================================*
*                                                                      *
*        A U F B A U  section                                          *
*                                                                      *
*======================================================================*
*                                                                      *
         If (Aufb.and..Not.Teee) Then
*
            If(iter.ne.1 .and. Abs(EDiff).lt.0.002D0 .and.
     &         .Not.AufBau_Done) Then
               Aufbau_Done=.True.
               Diis=Diis_Save
               DiisTh=DiisTh_save
            End If
*
         Else If (Aufb.and.Teee) Then
*
* If one iteration has been performed with the end temperature and we
* have a stable aufbau Slater determinant, terminate.
* The temperature is scaled down each iteration until it reaches the end
* temperature
*
* !!! continue here
            If(Abs(RTemp-TStop).le.1.0d-3
     &         .and. iAufOK.eq.1
     &         .and. Abs(Ediff).lt.1.0d-2 ) Then
               Aufbau_Done=.True.
               Diis=Diis_Save
               DiisTh=DiisTh_save
            End If
            RTemp=TemFac*RTemp
*           RTemp=Max(TStop,RTemp)
            TSTop=Min(TStop,RTemp)
         End If
*                                                                      *
*======================================================================*
*======================================================================*
*                                                                      *
*------- Print out necessary things
*
         TCP2=seconds()
         CpuItr = TCP2 - TCP1
         Call PrIte(iOpt.ge.2,CMO,mBB,nD,Ovrlp,mBT,OccNo,mmB)
*
*----------------------------------------------------------------------*
         Call Scf_Mcontrol(iter)
*----------------------------------------------------------------------*
*
*------- Check that two-electron energy is positive.
*        The LDF integrals may not be positive definite, leading to
*        negative eigenvalues of the ERI matrix - i.e. attractive forces
*        between the electrons! When that occurs the SCF iterations
*        still might find a valid SCF solution, but it will obviously be
*        unphysical and we might as well just stop here.
*
         If (Neg2_Action.ne.'CONT') Then
            If (E2V.lt.E2VTolerance) Then
               Call WarningMessage(2,
     &                        'WfCtl_SCF: negative two-electron energy')
               Write(6,'(A,1P,D25.10)') 'Two-electron energy E2V=',E2V
               Call xFlush(6)
               If (Neg2_Action.eq.'STOP') Call Abend()
            End If
         End If
*                                                                      *
************************************************************************
************************************************************************
!                                                                      *
!------- Perform another iteration if necessary
!
!        Convergence criteria are
!
!        Either,
!
!        1) it is not the first iteration
!        2) the absolute energy change is smaller than EThr
!        3) the absolute Fock matrix change is smaller than FThr
!        4) the absolute density matrix change is smaller than DThr,
!           and
!        5) step is NOT a Quasi NR step
!
!        or
!
!        1) step is a Quasi NR step, and
!        2) DltNrm.le.DltNth
!
*                                                                      *
************************************************************************
*                                                                      *
         If (EDiff>0.0.and..Not.Reset) EDiff=Ten*EThr
         If (iter.ne.1             .AND.
     &       (Abs(EDiff).le.EThr)  .AND.
     &       (Abs(FMOMax).le.FThr) .AND.
     &       (((Abs(DMOMax).le.DThr).AND.(iOpt.lt.2))
     &       .OR.
     &       ((DltNrm.le.DltNTh).AND.iOpt.ge.2))
     &      ) Then
*                                                                      *
************************************************************************
*                                                                      *
            If(Aufb) Then
               WarnPocc=.true.
            End If
*
*           If calculation with two different sets of parameters reset
*           the convergence parameters to default values. This is
*           possibly used for direct SCF and DFT.
*
            If (Reset) Then
               Call Reset_some_stuff()
               Cycle
            End If
*
*           Here if we converged!
*
*           Branch out of the iterative loop! Done!!!
*
            Converged=.True.
            Exit
*                                                                      *
************************************************************************
*                                                                      *
         Else If(Aufb.and.AufBau_Done) Then
*                                                                      *
************************************************************************
*                                                                      *
*           Here if Aufbau stage is finished and the calculation will
*           continue normally.
*
            Aufb=.FALSE.
            If (jPrint.ge.2) Then
               Write(6,*)
               Write(6,'(6X,A)') ' Fermi aufbau procedure completed!'
            End If
            IterX=-2
            If(iUHF.eq.0) Then
               iOffOcc=0
               Do iSym = 1, nSym
                  Do iBas = 1, nOrb(iSym)
                     iOffOcc=iOffOcc+1
                     If(iBas.le.nOcc(iSym,1)) Then
                        OccNo(iOffOcc,1)=2.0d0
                     Else
                        OccNo(iOffOcc,1)=0.0d0
                     End If
                  End Do
               End Do
               If (jPrint.ge.2) Then
                Write (6,'(6X,A,8I5)') 'nOcc=',
     &                (nOcc(iSym,1),iSym=1,nSym)
                Write (6,*)
               End If
            Else
               iOffOcc=0
               Do iSym = 1, nSym
                  Do iBas = 1, nOrb(iSym)
                     iOffOcc=iOffOcc+1
                     If(iBas.le.nOcc(iSym,1)) Then
                        OccNo(iOffOcc,1)=1.0d0
                     Else
                        OccNo(iOffOcc,1)=0.0d0
                     End If
                     If(iBas.le.nOcc(iSym,2)) Then
                        OccNo(iOffOcc,2)=1.0d0
                     Else
                        OccNo(iOffOcc,2)=0.0d0
                     End If
                  End Do
               End Do
               If (jPrint.ge.2) Then
                  Write (6,'(6X,A,8I5)')
     *              'nOcc(alpha)=',(nOcc(iSym,1),iSym=1,nSym)
                  Write (6,'(6X,A,8I5)')
     *              'nOcc(beta) =',(nOcc(iSym,2),iSym=1,nSym)
                  Write (6,*)
               End If

            End If
*
*---------- Now when nOcc is known compute standard sizes of arrays.
*
            Call Setup()
*
         End If
*                                                                      *
*======================================================================*
*                                                                      *
*                                                                      *
*     End of iteration loop
*
      End Do ! iter_
*                                                                      *
*======================================================================*
*                                                                      *

      If (Converged) Then

         If (jPrint.ge.2) Then
            Write (6,*)
            Write (6,'(6X,A,I3,A)')' Convergence after',
     &                iter, ' Macro Iterations'
         End If

      Else

!        Here if we didn't converge or if this was a forced one
!        iteration  calculation.

         iter=iter-1
         If(nIter(nIterP).gt.1) Then
            If (jPrint.ge.1) Then
               Write(6,*)
               Write(6,'(6X,A,I3,A)')
     &              ' No convergence after',iter,' Iterations'
            End If
            iTerm = _RC_NOT_CONVERGED_
         Else
            If (jPrint.ge.2) Then
               Write(6,*)
               Write(6,'(6X,A)') ' Single iteration finished!'
            End If
         End If
*
         If (jPrint.ge.2) Write(6,*)
*
         If(Reset) Then
            If(DSCF.and.KSDFT.eq.'SCF') Call Reset_Thresholds()
            If(KSDFT.ne.'SCF'.and..Not.One_Grid) then
              Call Reset_NQ_grid()
            End If
         End If

      End If
*
***********************************************************
*                      S   T   O   P                      *
***********************************************************
*
      If (jPrint.ge.2) Then
         Call CollapseOutput(0,'Convergence information')
         Write(6,*)
      End If
*---  Compute total spin (if UHF)
      If(iUHF.eq.0) Then
         s2uhf=0.0d0
      Else
         Call s2calc(CMO(1,1),CMO(1,2),Ovrlp,
     &               nOcc(1,1),nOcc(1,2),nBas,nOrb, nSym,
     &               s2uhf)
      End If
*---
      Call Add_Info('SCF_ITER',[DBLE(Iter)],1,8)
c     Call Scf_XML(0)
*
      Call KiLLs()
*
*     If the orbitals are generated with orbital rotations in
*     RotMOs we need to generate the canonical orbitals.
*
      If (iOpt.ge.2) Then
*
*---    Generate canonical orbitals
*
         Call TraFck(Fock,nBT,CMO,nBO,.TRUE.,FMOMax,EOrb,nnO,
     &               Ovrlp,nD)
*
*        Transform density matrix to MO basis
*
         Call MODens(Dens,Ovrlp,nBT,nDens,CMO,nBB,nD)
*
      End If
*
*---- Compute correct orbital energies
*
      Call MkEorb(Fock,nBT,CMO,nBB,EOrb,nnB,nSym,nBas,nOrb,nD)
*
*     Put orbital coefficients and energies on the runfile.
*
      Call Put_darray('SCF orbitals',   CMO(1,1),nBB)
      Call Put_darray('OrbE',   Eorb(1,1),nnB)
      If (nD.eq.2) Then
         Call Put_darray('SCF orbitals_ab',   CMO(1,2),nBB)
         Call Put_darray('OrbE_ab',   Eorb(1,2),nnB)
      End If
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*---  Finalize Cholesky information if initialized
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
      If (DoCholesky) Then
         If (DoLDF) Then
            Call LDF_UnsetIntegralPrescreeningInfo()
            Call LDF_X_Final(.True.,irc)
            If (irc.ne.0) Then
               Call WarningMessage(2,
     &                   'WfCtl: non-zero return code from LDF_X_Final')
               Call LDF_Quit(1)
            End If
         Else
            Call Cho_X_Final(irc)
            If (irc.ne.0) Then
               Call WarningMessage(2,
     &                             'WfCtl. Non-zero rc in Cho_X_Final.')
               Call Abend
            End If
         End If
      End If
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
#ifdef _MSYM_
      If (MSYMON) Then
         Call fmsym_release_context(msym_ctx)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Exit                                                             *
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
      Call SavTim(3,TCpu2-TCpu1,TWall2-TWall1)
      TimFld( 2) = TimFld( 2) + (TCpu2 - TCpu1)

*                                                                      *
************************************************************************
*                                                                      *
      Contains
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Save_Orbitals()
      Integer iSym
*
*---  Save the new orbitals in case the SCF program aborts

      iTrM = 1
      iCMO = 1
      Do iSym = 1, nSym
         nBs = nBas(iSym)
         nOr = nOrb(iSym)
         lth = nBs*nOr
         Do iD = 1, nD
            Call DCopy_(lth,CMO(iCMO,iD),1,TrM(iTrM,iD),1)
            Call FZero(TrM(iTrm+nBs*nOr,iD),nBs*(nBs-nOr))
         End Do
         iTrM = iTrM + nBs*nBs
         iCMO = iCMO + nBs*nOr
      End Do
*
      If(iUHF.eq.0) Then
         OrbName='SCFORB'
         Note='*  intermediate SCF orbitals'

         Call WrVec_(OrbName,LuOut,'CO',iUHF,nSym,nBas,nBas,
     &               TrM(1,1), Dummy,OccNo(1,1), Dummy,
     &               Dummy,Dummy, iDummy,Note,2)
         Call Put_darray('SCF orbitals',TrM(1,1),nBB)
         Call Put_darray('OrbE',Eorb(1,1),nnB)
         If(.not.Aufb) Then
            Call Put_iarray('SCF nOcc',nOcc(1,1),nSym)
         End If
      Else
         OrbName='UHFORB'
         Note='*  intermediate UHF orbitals'
         Call WrVec_(OrbName,LuOut,'CO',iUHF,nSym,nBas,nBas,
     &               TrM(1,1), TrM(1,2),OccNo(1,1),OccNo(1,2),
     &               Dummy,Dummy, iDummy,Note,3)
         Call Put_darray('SCF orbitals',   TrM(1,1),nBB)
         Call Put_darray('SCF orbitals_ab',TrM(1,2),nBB)
         Call Put_darray('OrbE',   Eorb(1,1),nnB)
         Call Put_darray('OrbE_ab',Eorb(1,2),nnB)
         If(.not.Aufb) Then
            Call Put_iarray('SCF nOcc',   nOcc(1,1),nSym)
            Call Put_iarray('SCF nOcc_ab',nOcc(1,2),nSym)
         End If
      End If
      End Subroutine Save_Orbitals
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Reset_some_Stuff()
*                                                                      *
*-------       Reset thresholds for direct SCF procedure
*
               Reset=.False.
*
*---------------
               If(iOpt.eq.2) Then
                  iOpt = 1        ! True if step is QNR
                  QNR1st=.TRUE.
               End If
               iterso=0
               If(Reset_Thresh) Call Reset_Thresholds()
               If(KSDFT.ne.'SCF') Then
                  If (.Not.One_Grid) Then
                     iterX=0
                     Call Reset_NQ_grid()
*                    Call PrBeg(Meth_)
                  End If
                  If ( iOpt.eq.0 ) kOptim=1
               End If
      End Subroutine Reset_some_Stuff
*                                                                      *
************************************************************************
*                                                                      *
      End SubRoutine WfCtl_SCF_Internal
