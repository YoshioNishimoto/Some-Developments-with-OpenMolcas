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
* Copyright (C) Naoki Nakatani                                         *
************************************************************************
*  DMRGCtl
*
*> @brief
*>   DMRG Control
*> @author N. Nakatani
*>
*> @details
*> Many functionality is same as ::CICtl (copy most part of ::CICtl to here)
*> To use DMRG as CAS-CI solver.
*> \p IRst = ``0`` is used for the first and the last (DMRG-) CI calculations
*> to fully optimize the DMRG wavefunction.
*> \p IRst = ``1`` enables to restart DMRG calculation from last iteration
*>
*> @param[in]     CMO    MO coefficients
*> @param[out]    D      Average 1-dens matrix
*> @param[out]    DS     Average spin 1-dens matrix
*> @param[out]    P      Average symm. 2-dens matrix
*> @param[out]    PA     Average antisymm. 2-dens matrix
*> @param[out]    FI     Fock matrix from inactive density
*> @param[in,out] D1I    Inactive 1-dens matrix
*> @param[in,out] D1A    Active 1-dens matrix
*> @param[in]     TUVX   Active 2-el integrals
*> @param[in]     IFINAL Calculation status switch
*> @param[in]     IRst   DMRG restart status switch
************************************************************************
#if defined (_ENABLE_BLOCK_DMRG_) || defined (_ENABLE_CHEMPS2_DMRG_) || defined (_ENABLE_DICE_SHCI_)
      Subroutine DMRGCtl(CMO,D,DS,P,PA,FI,D1I,D1A,TUVX,IFINAL,IRst)

      use wadr, only: FMO
      use stdalloc, only: mma_allocate, mma_deallocate
      use rctfld_module
      Use casvb_global, Only: ifvb
      use rasscf_lucia, only: PAtmp, Pscr, Ptmp, DStmp, Dtmp
!     use sxci, only: IDXSX

      Implicit Real* 8 (A-H,O-Z)
      Dimension CMO(*),D(*),DS(*),P(*),PA(*),FI(*),D1I(*),D1A(*),
     &          TUVX(*)
c     Logical Exist
      Logical Do_ESPF

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "output_ras.fh"
      Character*16 ROUTINE
      Parameter (ROUTINE='DMRGCTL ')
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "timers.fh"
#include "gas.fh"
#include "pamint.fh"
*PAM05      SymProd(i,j)=1+iEor(i-1,j-1)
*
C Local print level (if any)
      IPRLEV=IPRLOC(3)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

* set up flag 'IFCAS' for GAS option, which is set up in gugatcl originally.
* IFCAS = 0: This is a CAS calculation
* IFCAS = 1: This is a RAS calculation - This might cause an error in DMRG-CASSCF.
*
      if(iDoGas) call setsxci
      If(IPRLEV.gt.DEBUG ) then
        Write(LF,*)
        Write(LF,*) ' Enter DMRG section'
        Write(LF,*) ' =================='
        Write(LF,*)
        Write(LF,*) ' iteration count =',ITER
      End If
      if(ifinal.ne.0) PamGen1=.True.
*
* SOME DIRTY SETUPS
*
      S=0.5D0*DBLE(ISPIN-1)
*
* COMPUTE ONE ELECTRON INTEGRALS IN MO BASIS
* AND ADD CORE INTERACTION
*
* FMO FOCK MATRIX IN MO-BASIS
* LW2: 1-PARTICLE DENSITY MATRIX ALSO USED IN MO/AO TRANSFORMATION
*
      CALL mma_allocate(FMO,NACPAR,Label='FMO')
      Call DecideOnESPF(Do_ESPF)
      If ( lRf .or. KSDFT.ne.'SCF' .or. Do_ESPF) THEN
*
* In case of a reaction field in combination with an average CAS
* select the potential of the appropriate state.
*
        jDisk = IADR15(3)
        Do i=1,IPCMROOT-1
          Call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
          Call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
          Call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
          Call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
        End Do
*
        CALL GETMEM('D1A_FULL','ALLO','REAL',LRCT_F,NTOT2)
        CALL GETMEM('D1S_FULL','ALLO','REAL',LRCT_FS,NTOT2)
        If (IFINAL.eq.0) Then
*
* Use normal MOs
*
           CALL GETMEM('D1A_RCT','ALLO','REAL',LRCT,NACPAR)
           CALL GETMEM('P2MO','ALLO','REAL',ipP2MO,NACPR2)
*
* Get the total density in MOs
*
           Call DDafile(JOBIPH,2,Work(LRCT),NACPAR,jDisk)
           Call Put_dArray('D1mo',Work(LRCT),NACPAR)  ! Put it on the RUNFILE
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(Work(LRCT))
* Transform to AOs
           Call Get_D1A_RASSCF(CMO,WORK(LRCT),WORK(LRCT_F))
*
* Get the spin density in MOs
*
           IF (NACTEL.EQ.0) THEN
             CALL DCOPY_(NTOT2,[0.0D0],0,WORK(LRCT_FS),1)
           ELSE
             CALL GETMEM('D1S_RCT','ALLO','REAL',LRCT_S,NACPAR)
             Call DDafile(JOBIPH,2,Work(LRCT_S),NACPAR,jDisk)
             IF ( NASH(1).NE.NAC ) CALL DBLOCK(Work(LRCT_S))
* Transform to AOs
             Call Get_D1A_RASSCF(CMO,WORK(LRCT_S),WORK(LRCT_FS))
             CALL GETMEM('D1S_RCT','FREE','REAL',LRCT_S,NACPAR)
           END IF
*
* Get the 2-particle density in MO
*
           Call DDafile(JOBIPH,2,Work(ipP2MO),NACPR2,jDisk)
           Call Put_dArray('P2mo',Work(ipP2MO),NACPR2) ! Put it on the RUNFILE
*
           CALL SGFCIN(CMO,FMO,FI,D1I,Work(LRCT_F),Work(LRCT_FS))
*
           CALL GETMEM('P2MO','FREE','REAL',ipP2MO,NACPR2)
           CALL GETMEM('D1A_RCT','FREE','REAL',LRCT,NACPAR)
*
        Else
*
* Here the pseudo-natural orbitals are in CMO and we need to
* get the D1A of the selected state in this basis.
*
*
* Compute the density of the particular state
*

           Call mma_allocate(Dtmp,NAC**2,Label='Dtmp')
           Call mma_allocate(DStmp,NAC**2,Label='DStmp')
           Call mma_allocate(Ptmp,NACPR2,Label='Ptmp')
           If ( NAC.ge.1 ) Then
              If (NACTEL.eq.0) THEN
                 Dtmp(:)=0.0D0
                 DStmp(:)=0.0D0
                 Ptmp(:)=0.0D0
              Else
* load back 1- and 2-RDMs from previous DMRG run
                 NACT4=NAC**4
                 Call mma_allocate(PAtmp,NACPR2,Label='PAtmp')
                 CALL mma_allocate(Pscr,NACT4,Label='Pscr')
#ifdef _ENABLE_BLOCK_DMRG_
                 CALL block_densi_rasscf(IPCMRoot,Dtmp,DStmp,
     &                                   Ptmp,PAtmp,Pscr)
#elif _ENABLE_CHEMPS2_DMRG_
                 CALL chemps2_densi_rasscf(IPCMRoot,Dtmp,DStmp,
     &                                     Ptmp,PAtmp,Pscr)
#elif _ENABLE_DICE_SHCI_
                 CALL dice_densi_rasscf(IPCMRoot,Dtmp,DStmp,
     &                                  Ptmp,PAtmp,Pscr)
#endif

* NN.14 NOTE: IFCAS must be 0 for DMRG-CASSCF
c                If (IFCAS.GT.2) Call CISX(IDXSX,Dtmp,
c    &                                     DStmp,Ptmp,PAtmp,Pscr)
                 Call mma_deallocate(Pscr)
                 Call mma_deallocate(PAtmp)
              EndIf
*
           Else
              Dtmp(:)=0.0D0
              DStmp(:)=0.0D0
              Ptmp(:)=0.0D0
           End If
* Modify the symmetric 2-particle density if only partial
* "exact exchange" is included.
c          n_Det=2
c          n_unpaired_elec=(iSpin-1)
c          n_paired_elec=nActEl-n_unpaired_elec
c          If(n_unpaired_elec+n_paired_elec/2.eq.nac) n_Det=1
           If (ExFac.ne.1.0D0) Call Mod_P2(Ptmp,NACPR2,
     &                                   Dtmp,NACPAR,
     &                                   DStmp,ExFac,n_Det)
*
           Call Put_dArray('P2mo',Ptmp,NACPR2) ! Put it on the RUNFILE
*
           Call mma_deallocate(Ptmp)
*
           Call Put_dArray('D1mo',Dtmp,NACPAR) ! Put it on the RUNFILE
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(Dtmp)
           Call Get_D1A_RASSCF(CMO,Dtmp,Work(LRCT_F))
*
           IF ( NASH(1).NE.NAC ) CALL DBLOCK(DStmp)
           Call Get_D1A_RASSCF(CMO,DStmp,Work(LRCT_FS))
*
           Call mma_deallocate(Dtmp)
           Call mma_deallocate(DStmp)
*
           Call SGFCIN(CMO,FMO,FI,D1I,Work(LRCT_F),Work(LRCT_FS))
*
        End If
        CALL GETMEM('D1S_FULL','FREE','REAL',LRCT_FS,NTOT2)
        CALL GETMEM('D1A_RCT','FREE','REAL',LRCT_F,NTOT1)
*
      ELSE
*
* Normal case
*
*
*
*
        CALL GETMEM('TmpDS' ,'Allo','REAL',ipTmpDS ,NACPAR)
        CALL GETMEM('TmpD1S','Allo','REAL',ipTmpD1S,NTOT2 )
        call dcopy_(NACPAR,DS,1,Work(ipTmpDS),1)
        IF ( NASH(1).NE.NAC ) CALL DBLOCK(Work(ipTmpDS))
        Call Get_D1A_RASSCF(CMO,Work(ipTmpDS),Work(ipTmpD1S))
        CALL GETMEM('TmpDS' ,'Free','REAL',ipTmpDS ,NACPAR)
*
        CALL SGFCIN(CMO,FMO,FI,D1I,D1A,Work(ipTmpD1S))
        CALL GETMEM('TmpD1S','Free','REAL',ipTmpD1S,NTOT2 )
*
      END IF
*
      If (IfVB.eq.2) GoTo 9000
*
* SOLVE DMRG WAVEFUNCTION
*
      if(IfVB.eq.1)then
* NN.14 FIXME: I'm not sure whether this option should work?
        call cvbmn_rvb(max(ifinal,1))
      else
        If (KSDFT(1:3).ne.'SCF'.
     &      and.DFTFOCK(1:4).eq.'DIFF'.and.nac.ne.0) Then
          nTmpPUVX=nFint
          Call GetMem('TmpPUVX','Allo','Real',ipTmpPUVX,nTmpPUVX)
          Call GetMem('TmpTUVX','Allo','Real',ipTmpTUVX,NACPR2)
          Call dCopy_(NACPR2,[0.0d0],0,Work(ipTmpTUVX),1)
          Call Get_dArray('DFT_TwoEl',Work(ipTmpPUVX),nTmpPUVX)
          Call Get_TUVX(Work(ipTmpPUVX),Work(ipTmpTUVX))
          Call DaXpY_(NACPR2,1.0d0,TUVX,1,Work(ipTmpTUVX),1)
#ifdef _ENABLE_BLOCK_DMRG_
          Call BlockCtl(FMO,Work(ipTmpTUVX),IFINAL,IRst)
#elif _ENABLE_CHEMPS2_DMRG_
          Call Chemps2Ctl(FMO,Work(ipTmpTUVX),IFINAL,IRst)
#elif _ENABLE_DICE_SHCI_
          Call DiceCtl(FMO,Work(ipTmpTUVX),IFINAL,IRst)
#endif

          Call GetMem('TmpTUVX','Free','Real',ipTmpTUVX,NACPR2)
          Call GetMem('TmpPUVX','Free','Real',ipTmpPUVX,nTmpPUVX)
        Else
#ifdef _ENABLE_BLOCK_DMRG_
          Call BlockCtl(FMO,TUVX,IFINAL,IRst)
#elif _ENABLE_CHEMPS2_DMRG_
          Call Chemps2Ctl(FMO,TUVX,IFINAL,IRst)
#elif _ENABLE_DICE_SHCI_
          Call DiceCtl(FMO,TUVX,IFINAL,IRst)
#endif
        End If
      endif

*
* C
* CALCULATE DENSITY MATRICES
* SAVE DENSITY MATRICES ON FILE
* COMPUTE AVERAGE DENSITY MATRICES
* C
*
* Dtmp: ONE-BODY DENSITY
* DStmp: ONE-BODY SPIN DENSITY
* Ptmp: SYMMETRIC TWO-BODY DENSITY
* PAtmp: ANTISYMMETRIC TWO-BODY DENSITY
*
      Call Timing(Rado_1,dum1,dum2,dum3)
      Zero = 0.0d0
      Call dCopy_(NACPAR,[Zero],0,D,1)
      Call dCopy_(NACPAR,[Zero],0,DS,1)
      Call dCopy_(NACPR2,[Zero],0,P,1)
      Call dCopy_(NACPR2,[Zero],0,PA,1)
      CALL mma_allocate(Dtmp,NAC**2,Label='Dtmp')
      CALL mma_allocate(DStmp,NAC**2,Label='DStmp')
      CALL mma_allocate(Ptmp,NACPR2,Label='Ptmp')
      CALL mma_allocate(PAtmp,NACPR2,Label='PAtmp')
      jDisk = IADR15(3)

      Do jRoot = 1,lRoots
* load density matrices from DMRG run
        If ( NAC.ge.1 ) Then
          NACT4=NAC**4
          CALL mma_allocate(Pscr,NACT4,Label='Pscr')
#ifdef _ENABLE_BLOCK_DMRG_
          CALL block_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#elif _ENABLE_CHEMPS2_DMRG_
          CALL chemps2_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#elif _ENABLE_DICE_SHCI_
          CALL dice_densi_rasscf(jRoot,Dtmp,DStmp,Ptmp,PAtmp,Pscr)
#endif
          Call mma_deallocate(Pscr)
        EndIf
* Modify the symmetric 2-particle density if only partial
* "exact exchange" is included.
c        n_Det=2
c        n_unpaired_elec=(iSpin-1)
c        n_paired_elec=nActEl-n_unpaired_elec
c        If(n_unpaired_elec+n_paired_elec/2.eq.nac) n_Det=1
c
c           Write(LF,*) ' iSpin=',iSpin
c           Write(LF,*) ' n_unpaired_elec',n_unpaired_elec
c           Write(LF,*) ' n_paired_elec',  n_paired_elec
c           Write(LF,*) ' n_unpaired_elec+n_paired_elec/2',
c     &                  n_unpaired_elec+n_paired_elec/2
c           Write(LF,*) ' n_Det=',n_Det
c
c
        If (ExFac.ne.1.0D0) Call Mod_P2(Ptmp,NACPR2,
     &                                Dtmp,NACPAR,
     &                                DStmp,ExFac,n_Det)

* update average density matrices
        Scal = 0.0d0
        Do kRoot = 1,nRoots
          If ( iRoot(kRoot).eq.jRoot ) then
            Scal = Weight(kRoot)
          End If
        End Do
        call daxpy_(NACPAR,Scal,Dtmp,1,D,1)
        call daxpy_(NACPAR,Scal,DStmp,1,DS,1)
        call daxpy_(NACPR2,Scal,Ptmp,1,P,1)
        call daxpy_(NACPR2,Scal,PAtmp,1,PA,1)
* save density matrices on disk
        Call DDafile(JOBIPH,1,Dtmp,NACPAR,jDisk)
        Call DDafile(JOBIPH,1,DStmp,NACPAR,jDisk)
        Call DDafile(JOBIPH,1,Ptmp,NACPR2,jDisk)
        Call DDafile(JOBIPH,1,PAtmp,NACPR2,jDisk)
      End Do

      Call mma_deallocate(PAtmp)
      Call mma_deallocate(Ptmp)
      Call mma_deallocate(DStmp)
      Call mma_deallocate(Dtmp)
*
* C
* PREPARE DENSITY MATRICES AS USED BY THE SUPER CI SECTION
* C
*
* print matrices
      IF ( IPRLEV.GE.INSANE  ) THEN
        CALL TRIPRT('Averaged one-body density matrix, D',
     &              ' ',D,NAC)
        CALL TRIPRT('Averaged one-body spin density matrix, DS',
     &              ' ',DS,NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',P,NACPAR)
        CALL TRIPRT('Averaged antisymmetric two-body density matrix,PA',
     &              ' ',PA,NACPAR)
      END IF
      IF ( NASH(1).NE.NAC ) CALL DBLOCK(D)
      Call Timing(Rado_2,dum1,dum2,dum3)
      Rado_2 = Rado_2 - Rado_1
      Rado_3 = Rado_3 + Rado_2
*
      Call mma_deallocate(FMO)

 9000 Continue
*
*     For RF calculations make sure that the we are following the
*     correct root.
*
*     In the current implementation the overlap between the CI vectors
*     of different macro iterations is used. This criterion stricktly
*     only hold if the orbitals are not changed in between the
*     interations, to make sure that this approximately holds the
*     comparision is only to be considered to be valid if the rotmax
*     parameter is below an empirical threshold. In the future the
*     procedure for indentifying root flipping has to be made more
*     robust.
*
CSVC: if CISElect is used, roots should be traced regardless of orbital
C     rotation, so in that case ignore the automatic tracing and follow
C     the relative CISE root given in the input by the 'CIRF' keyword.
*
* ======================================================================
* NN.14 FIXME:
*     The overlap b/w the DMRG wave of different macro iterations
*     hasn't yet been implemented. Eventually this must be fixed to
*     chose the correct root for RF calculation with DMRG-CASSCF.
*     For the time, just skip following.
* ======================================================================
*
c     If (lRF.and.KeyCISE.and.KeyCIRF) Then
c       JPCMROOT=IPCMROOT
c       IPCMROOT=IROOT(ICIRFROOT)
c       Call Put_iScalar("RF CASSCF root",IPCMROOT)
c       If (JPCMROOT.ne.IPCMROOT) Then
c         Write (6,'(1X,A,I3,A,I3)') 'RF Root has flipped from ',
c    &                 JPCMROOT, ' to ',IPCMROOT
c       End If
c     Else If (lRF) Then
c       Call Qpg_iScalar('RF CASSCF root',Exist)
c       If (.NOT.Exist) Then
*
*          We are here since we are using the default values.
*
c          Call Put_iScalar("RF CASSCF root",IPCMROOT)
c          Call Put_iScalar("RF0CASSCF root",IPCMROOT)
c       End If
*
c       Call Allocate_Work(ipRF,nConf)
c       Call Qpg_dArray("RF CASSCF Vector",Exist,mConf)
*       Write (6,*) 'Exist=',Exist
c       If (Exist
c    &      .and. mConf .eq. nConf
c    &      .and. iFinal.ne.2
c    &      .and. (ABS(RotMax).lt.1.0D-3 .or. KeyCISE)
c    &     ) Then
c          Call Get_dArray("RF CASSCF Vector",Work(ipRF),nConf)
c          rNorm=Sqrt(DDot_(nConf,Work(ipRF),1,Work(ipRF),1))
*          Write (6,*) 'rNorm=',rNorm
c          JPCMROOT=IPCMROOT
c          If (rNorm.gt.1.0D-10) Then
c             Call Allocate_Work(ipTemp,nConf)
c             rMax=0.0D0
c             jDisk = IADR15(4)
c             Do i = 1, lRoots
c                Call DDafile(JOBIPH,2,Work(ipTemp),nConf,jDisk)
c                qMax=Abs(DDot_(nConf,Work(ipTemp),1,Work(ipRF),1))
*                Write (6,*) 'qMax=',qMax
c                If (qMax.gt.rMax .and.
c    &               qMax.gt.0.5D0) Then
c                   rMax = qMax
c                   JPCMROOT=i
c                End If
c             End Do
c             Call Free_Work(ipTemp)
c          End If
c       Else
c          JPCMROOT=IPCMROOT
c       End If
*
c       If (JPCMROOT.ne.IPCMROOT) Then
c          Write (6,*) ' RF Root has flipped from ',IPCMROOT, ' to ',
c    &                                              JPCMROOT
c          IPCMROOT=JPCMROOT
c          Call Put_iScalar("RF CASSCF root",IPCMROOT)
c       End If
*
c       jDisk = IADR15(4)
c       Do i=1,IPCMROOT-1
c         Call DDafile(JOBIPH,0,rdum,nConf,jDisk)
c       End Do
c       Call DDafile(JOBIPH,2,Work(ipRF),nConf,jDisk)
c       Call Put_dArray("RF CASSCF Vector",Work(ipRF),nConf)
c       Call Free_Work(ipRF)
c     End If
*
      Return
      End

* _ENABLE_BLOCK_DMRG_
#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_DMRGCtl()
      End
#endif
