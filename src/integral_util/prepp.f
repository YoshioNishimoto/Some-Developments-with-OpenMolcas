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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine PrepP()
************************************************************************
*                                                                      *
* Object: to set up the handling of the 2nd order density matrix.      *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January '92                                              *
************************************************************************
      use iSD_data
      use aces_stuff
      use pso_stuff
      use index_arrays, only: iSO2Sh
      use Basis_Info, only: nBas
      use Sizes_of_Seward, only: S
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "etwas.fh"
#include "mp2alaska.fh"
#include "nsd.fh"
#include "dmrginfo_mclr.fh"
#include "nac.fh"
#include "mspdft.fh"
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
************ columbus interface ****************************************
#include "columbus_gamma.fh"
#include "setup.fh"
************************************************************************
      Integer nFro(0:7)
      Integer Columbus
      Character*8 RlxLbl,Method, KSDFT*80
      Logical lPrint
      Logical DoCholesky
      Real*8 CoefX,CoefR
      Character Fmt*60
      Real*8, Allocatable:: D1ao(:), D1AV(:), Tmp(:,:)
*     hybrid MC-PDFT things
      Logical Do_Hybrid
      Real*8  WF_Ratio,PDFT_Ratio
*
*...  Prologue

      iRout = 250
      iPrint = nPrint(iRout)
      lPrint=iPrint.ge.6
#ifdef _CD_TIMING_
      Call CWTIME(PreppCPU1,PreppWall1)
#endif
*
      Call StatusLine(' Alaska:',' Prepare the 2-particle matrix')
*
      iD0Lbl=1
      iComp=1
*
      lsa=.False.
      Gamma_On=.False.
      Gamma_mrcisd=.FALSE.
      lPSO=.false.
      Case_2C=.False.
      Case_3C=.False.
      Case_mp2=.False.

      nDens = 0
      Do 1 iIrrep = 0, nIrrep - 1
         nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
 1    Continue
*
*
*...  Get the method label
      Call Get_cArray('Relax Method',Method,8)
      Call Get_iScalar('Columbus',columbus)
      nCMo = S%n2Tot
      mCMo = S%n2Tot
      If (Method.eq. 'KS-DFT  ' .or.
     &    Method.eq. 'MCPDFT  ' .or.
     &    Method.eq. 'MSPDFT  ' .or.
     &    Method.eq. 'CASDFT  ' ) Then
         Call Get_iScalar('Multiplicity',iSpin)
         Call Get_cArray('DFT functional',KSDFT,80)
         Call Get_dScalar('DFT exch coeff',CoefX)
         Call Get_dScalar('DFT corr coeff',CoefR)
         ExFac=Get_ExFac(KSDFT)
         CoulFac=One
      Else
         iSpin=0
         ExFac=One
         CoulFac=One
      End If
*
*...  Check the wave function type
*
*                                                                      *
************************************************************************
*                                                                      *
      If ( Method.eq.'RHF-SCF ' .or.
     &     Method.eq.'UHF-SCF ' .or.
     &     Method.eq.'IVO-SCF ' .or.
     &     Method.eq.'MBPT2   ' .or.
     &     Method.eq.'KS-DFT  ' .or.
     &     Method.eq.'ROHF    ' ) then
         If (lPrint) Then
            Write (6,*)
            Write (6,'(2A)') ' Wavefunction type: ',Method
            If (Method.eq.'KS-DFT  ') Then
               Write (6,'(2A)') ' Functional type:   ',KSDFT
               Fmt = '(1X,A26,20X,F18.6)'
               Write(6,Fmt)'Exchange scaling factor',CoefX
               Write(6,Fmt)'Correlation scaling factor',CoefR
            End If
            Write (6,*)
         End If
         If(Method.eq.'MBPT2   ') Then
            Case_mp2=.true.
            Call DecideOnCholesky(DoCholesky)
            If(.not.DoCholesky) Then
               iSeed = 10
               LuGam = IsFreeUnit(iSeed)
               Write(FnGam,'(A6)') 'LuGam'
               Call DaName_MF_WA(LuGam,FnGam)
               Gamma_on=.True.
               Call Aces_Gamma()
            End If
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'Corr. WF' ) then
         If (lPrint) Then
            Write (6,*)
            Write (6,*)
     &         ' Wavefunction type: an Aces 2 correlated wavefunction'
            Write (6,*)
         End If
         Gamma_On=.True.
         Call Aces_Gamma()
       Else if (Method(1:7).eq.'MR-CISD' .and. Columbus.eq.1) then
************ columbus interface ****************************************
*do not reconstruct the two-particle density from the one-particle
*density or partial two-particle densities but simply read them from
*file

        nShell=mSkal
        nPair=nShell*(nShell+1)/2
        nQuad=nPair*(nPair+1)/2

*---- Allocate Table Of Content for half sorted gammas.
*
      Call mma_Allocate(G_Toc,nQuad+2,Label='G_Toc')

*  find free unit number
        lgtoc=61
        lgtoc=isFreeUnit(lgtoc)
        idisk=0
*  read table of contents of half-sorted gamma file
         Call DaName(lgtoc,'gtoc')
         Call ddafile(lgtoc,2,G_Toc,nQuad+2,idisk)
         Call Daclos(lgtoc)
         n=int(G_Toc(nQuad+1))
         lbin=int(G_Toc(nQuad+2))
         if (n.ne.nQuad) then
           Call WarningMessage(2,'n.ne.nQuad')
           Write (6,*) 'n,nQuad=',n,nQuad
           Call Abend()
         endif

        Gamma_On=.True.
        Gamma_mrcisd=.TRUE.
*       open gamma file
         LuGamma=60
         LuGamma=isfreeunit(LuGamma)
*        closed in closep
         Call DaName_MF(LuGamma,'GAMMA')
*  allocate space for bins
         Call mma_Allocate(Bin,2,lBin,Label='Bin')
*  compute SO2cI array
         Call mma_Allocate(SO2cI,2,nSOs,Label='SO2cI')
         call Mk_SO2cI(SO2cI,iSO2Sh,nsos)
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'RASSCF  ' .or.
     &          Method.eq.'CASSCF  ' .or.
     &          Method.eq.'GASSCF  ' .or.
     &          Method.eq.'MCPDFT  ' .or.
     &          Method.eq.'MSPDFT  ' .or.
     &          Method.eq.'DMRGSCF ' .or.
     &          Method.eq.'CASDFT  ') then
*
         Call Get_iArray('nAsh',nAsh,nIrrep)
         nAct = 0
         Do iIrrep = 0, nIrrep-1
            nAct = nAct + nAsh(iIrrep)
         End Do
         If (nAct.gt.0) lPSO=.true.
*
         nDSO = nDens
         mIrrep=nIrrep
         Call ICopy(nIrrep,nBas,1,mBas,1)
         If (lPrint) Then
            Write (6,*)
            Write (6,'(2A)') ' Wavefunction type: ', Method
            If (Method.eq.'CASDFT  ' .or. Method.eq.'MCPDFT  ')
     &         Write (6,'(2A)') ' Functional type:   ',KSDFT
            If (Method.eq.'MSPDFT  ')
     &         Write (6,'(2A)') ' MS-PDFT Functional type:   ',KSDFT
            Write (6,*)
         End If
         If (method.eq.'MCPDFT  ') lSA=.true.
         If (method.eq.'MSPDFT  ') then
          lSA=.true.
          dogradmspd=.true.
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else if ( Method.eq.'CASSCFSA' .or.
     &          Method.eq.'DMRGSCFS' .or.
     &          Method.eq.'GASSCFSA' .or.
     &          Method.eq.'RASSCFSA' .or.
     &          Method.eq.'CASPT2  ') then
         Call Get_iArray('nAsh',nAsh,nIrrep)
         nAct = 0
         Do iIrrep = 0, nIrrep-1
            nAct = nAct + nAsh(iIrrep)
         End Do
         If (nAct.gt.0) lPSO=.true.
         nDSO = nDens
         Call Get_iScalar('SA ready',iGo)
         If (iGO.eq.1) lSA=.true.
         mIrrep=nIrrep
         Call ICopy(nIrrep,nBas,1,mBas,1)
         If (lPrint) Then
            Write (6,*)
            If (lSA) Then
               Write (6,'(2A)') ' Wavefunction type: State average ',
     &                            Method(1:6)
            Else
               Write (6,'(2A)') ' Wavefunction type: ', Method
            End If
            Write (6,*)
         End If
         If (Method.eq.'CASPT2  ') Then
            Call DecideOnCholesky(DoCholesky)
C           If(.not.DoCholesky) Then
               Gamma_On=.True.
               !! It is opened, but not used actually. I just want to
               !! use the Gamma_On flag.
               !! Just to avoid error termination.
               !! Actual working arrys are allocated in drvg1.f
               LuGamma=60
               Call DaName_MF_WA(LuGamma,'GAMMA')
               If (DoCholesky) Call mma_allocate(G_Toc,1,Label='G_Toc')
               Call mma_allocate(SO2cI,1,1,Label='SO2cI')
               Call mma_allocate(Bin,1,1,Label='Bin')
C           End If
         End If
         Method='RASSCF  '
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Call WarningMessage(2,'Alaska: Unknown wavefuntion type')
         Write (6,*) 'Wavefunction type:',Method
         Write (6,*) 'Illegal type of wave function!'
         Write (6,*) 'ALASKA cannot continue.'
         Call Quit_OnUserError()
      End If
*
*...  Read the (non) variational 1st order density matrix
*...  density matrix in AO/SO basis
         nsa=1
         If (lsa) nsa=5
         If ( Method.eq.'MCPDFT  ') nsa=5
         If ( Method.eq.'MSPDFT  ') nsa=5
!AMS modification: add a fifth density slot
         mDens=nsa+1
         Call mma_allocate(D0,nDens,mDens,Label='D0')
         D0(:,:)=Zero
         Call mma_allocate(DVar,nDens,nsa,Label='DVar')
         if (.not.gamma_mrcisd)
     &      Call Get_dArray_chk('D1ao',D0(1,1),nDens)
*

         Call Get_D1ao_Var(DVar,nDens)
*
         Call mma_Allocate(DS,nDens,Label='DS')
         Call mma_Allocate(DSVar,nDens,Label='DSVar')
         If (Method.eq.'UHF-SCF ' .or.
     &       Method.eq.'ROHF    ' .or.
     &    (Method.eq.'KS-DFT  '.and.iSpin.ne.1) .or.
     &       Method.eq.'Corr. WF' ) Then
            Call Get_dArray_chk('D1sao',DS,nDens)
            Call Get_D1sao_Var(DSVar,nDens)
         Else
            DS   (:)=Zero
            DSVar(:)=Zero
         End If
*
*   This is necessary for the ci-lag
*
*
*     Unfold density matrix
*
************ columbus interface ****************************************
*do not modify the effective density matrices and fock matrices
      if (.not.gamma_mrcisd) then
      ij = -1
      Do 10 iIrrep = 0, nIrrep-1
         Do 11 iBas = 1, nBas(iIrrep)
            Do 12 jBas = 1, iBas-1
               ij = ij + 1
               D0   (1+ij,1) = Half*D0   (1+ij,1)
               DVar (1+ij,1) = Half*DVar (1+ij,1)
               DS   (1+ij) = Half*DS   (1+ij)
               DSVar(1+ij) = Half*DSVar(1+ij)
 12         Continue
            ij = ij + 1
 11      Continue
 10   Continue
      endif

      If (iPrint.ge.99) Then
         RlxLbl='D1AO    '
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],D0)
         RlxLbl='D1AO-Var'
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DVar)
         RlxLbl='DSAO    '
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DS)
         RlxLbl='DSAO-Var'
         Call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DSVar)
      End If

*
*...  Get the MO-coefficients
************ columbus interface ****************************************
*     not that the columbus mcscf MO coefficients have been written
*     to the RUNFILE !

         If (Method.eq.'UHF-SCF ' .or.
     &       Method.eq.'ROHF    ' .or.
     &    (Method.eq.'KS-DFT  '.and.iSpin.ne.1) .or.
     &       Method.eq.'Corr. WF'      ) Then
            nsa=2
         Else
            nsa=1
            If (lsa) nsa=2
         End If
         kCMO=nsa
         Call mma_allocate(CMO,mCMO,kCMO,Label='CMO')
         Call Get_dArray_chk('Last orbitals',CMO(:,1),mCMO)
         If (iPrint.ge.99) Then
            ipTmp1 = 1
            Do iIrrep = 0, nIrrep-1
               Call RecPrt(' CMO''s',' ',
     &                     CMO(ipTmp1,1),nBas(iIrrep),
     &                     nBas(iIrrep))
               ipTmp1 = ipTmp1 + nBas(iIrrep)**2
            End Do
         End If
*
*
*...  Get additional information in the case of a RASSCF wave function
*...  Get the number of inactive, active and frozen orbitals
************ columbus interface ****************************************
*  no need for MRCI gradient
         If (.not.lpso .or. gamma_mrcisd ) Goto 1000
         Call Get_iScalar('nSym',i)
         Call Get_iArray('nIsh',nIsh,i)
         Call Get_iArray('nAsh',nAsh,i)
         Call Get_iArray('nFro',nFro,i)
         If (iPrint.ge.99) Then
            Write (6,*) ' nISh=',nISh
            Write (6,*) ' nASh=',nASh
            Write (6,*) ' nFro=',nFro
         End If
         nAct = 0
         nTst = 0
         Do iIrrep = 0, nIrrep-1
!            write(*,*)"nAsh(iIrrep)",nAsh(iIrrep)  ! yma
            nAct = nAct + nAsh(iIrrep)
            nTst = nTst + nFro(iIrrep)
         End Do
         If (nTst.ne.0) Then
            Call WarningMessage(2,
     &                  '; No frozen orbitals are allowed!'//
     &                  '; ALASKA cannot continue;')
            Call Quit_OnUserError()
         End If
*
*...  Get the one body density for the active orbitals
*     (not needed for SA-CASSCF)
         nG1 = nAct*(nAct+1)/2
         nsa=1
         If (lsa) nsa=0
         mG1=nsa
         Call mma_allocate(G1,nG1,mG1,Label='G1')
         If (nsa.gt.0) Then
            Call Get_dArray_chk('D1mo',G1(:,1),nG1)
            If (iPrint.ge.99) Call TriPrt(' G1',' ',G1(:,1),nAct)
         End If
*
*...  Get the two body density for the active orbitals
         nG2 = nG1*(nG1+1)/2
         nsa=1
         if (lsa) nsa=2
         mG2=nsa
         Call mma_allocate(G2,nG2,mG2,Label='G2')
!       write(*,*) 'got the 2rdm, Ithink.'
         if(Method.eq.'MCPDFT  '.or.Method.eq.'MSPDFT  ') then
           Call Get_dArray_chk('P2MOt',G2,nG2)!PDFT-modified 2-RDM
         else
           Call Get_dArray_chk('P2mo',G2,nG2)
         end if
         If (iPrint.ge.99) Call TriPrt(' G2',' ',G2(1,1),nG1)
         If (lsa) Then

*  CMO1 Ordinary CMOs
*
*  CMO2 CMO*Kappa
*
           Call Get_dArray_chk('LCMO',CMO(:,2),mCMO)
           If (iPrint.ge.99) Then
            ipTmp1 = 1
            Do iIrrep = 0, nIrrep-1
               Call RecPrt('LCMO''s',' ',
     &                     CMO(ipTmp1,2),nBas(iIrrep),
     &                     nBas(iIrrep))
               ipTmp1 = ipTmp1 + nBas(iIrrep)**2
            End Do
           End If
*
* P are stored as
*                            _                     _
*   P1=<i|e_pqrs|i> + sum_i <i|e_pqrs|i>+<i|e_pqrs|i>
*   P2=sum_i <i|e_pqrs|i>
*
           Call Get_dArray_chk('PLMO',G2(:,2),nG2)
           ndim1=0
           if(doDMRG)then
             ndim0=0  !yma
             do i=1,8
               ndim0=ndim0+LRras2(i)
             end do
             ndim1=(ndim0+1)*ndim0/2
             ndim2=(ndim1+1)*ndim1/2
             do i=1,ng2
               if(i.gt.ndim2) G2(i,2)=0.0D0
             end do
           end if
           Call Daxpy_(ng2,One,G2(:,2),1,G2(:,1),1)
           If(iPrint.ge.99)Call TriPrt(' G2L',' ',G2(:,2),nG1)
           If(iPrint.ge.99)Call TriPrt(' G2T',' ',G2(:,1),nG1)
*
           Call Get_dArray_chk('D2av',G2(:,2),nG2)
           If (iPrint.ge.99) Call TriPrt('G2A',' ',G2(:,2),nG1)
*
*
*  Densities are stored as:
*
*       ipd0 AO:
*
*       D1 = inactive diagonal density matrix
*                                _                 _
*       D2 = <i|E_pq|i> + sum_i <i|E_pq|i>+<i|E_pq|i> + sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> - 1/2 D2
*
*       D3 = sum_i <i|E_pq|i> (active)
*
*       D4 = sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> (inactive)
*
*       G1 = <i|e_ab|i>
*       G2 = sum i <i|e_ab|i>
*
!************************
         RlxLbl='D1AO    '
!         Call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0)

           Call mma_allocate(Tmp,nDens,2,Label='Tmp')
           Call Get_D1I(CMO(1,1),D0(1,1),Tmp,nIsh,nBas,nIrrep)
           Call mma_deallocate(Tmp)

!************************
         RlxLbl='D1AO    '
!         Call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0)
*
           Call dcopy_(ndens,DVar,1,D0(1,2),1)
           If (.not.isNAC) call daxpy_(ndens,-Half,D0(1,1),1,D0(1,2),1)
!          RlxLbl='D1COMBO  '
!          Call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0(1,2))
*
*   This is necessary for the kap-lag
*
           nG1 = nAct*(nAct+1)/2
           Call mma_allocate(D1AV,nG1,Label='D1AV')
           Call Get_dArray_chk('D1av',D1AV,nG1)
           Call Get_D1A(CMO(1,1),D1AV,D0(1,3),
     &                 nIrrep,nbas,nish,nash,ndens)
           Call mma_deallocate(D1AV)
!************************
!          RlxLbl='D1AOA   '
!          Call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0(1,3))
*
           Call Get_dArray_chk('DLAO',D0(:,4),nDens)

*        Getting conditions for hybrid MC-PDFT
         Do_Hybrid=.false.
         CALL qpg_DScalar('R_WF_HMC',Do_Hybrid)
         If(Do_Hybrid) Then
          CALL Get_DScalar('R_WF_HMC',WF_Ratio)
          PDFT_Ratio=1.0d0-WF_Ratio
         End If
!ANDREW - modify D2: should contain only the correction pieces
         If ( Method.eq.'MCPDFT  ') then
!Get the D_theta piece
            Call mma_allocate(D1ao,nDens)
            Call Get_dArray_chk('D1ao',D1ao,ndens)
            ij = 0
            Do  iIrrep = 0, nIrrep-1
               Do iBas = 1, nBas(iIrrep)
                  Do jBas = 1, iBas-1
                     ij = ij + 1
                     D1ao(ij) = Half*D1ao(ij)
                  end do
                  ij = ij + 1
                end do
            end do
            call daxpy_(ndens,-1d0,D0(1,1),1,D1ao,1)
!           write(*,*) 'do they match?'
!           do i=1,ndens
!             write(*,*) d1ao(i),DO(i,3)
!           end do

            call daxpy_(ndens,-Half,D0(1,1),1,D0(1,2),1)
            call daxpy_(ndens,-1.0d0,D1ao,1,D0(1,2),1)
!ANDREW -   Generate new D5 piece:
            D0(:,5)=Zero
            call daxpy_(ndens,0.5d0,D0(1,1),1,D0(1,5),1)
            call daxpy_(ndens,1.0d0,D1ao,1,D0(1,5),1)

          if(do_hybrid) then
*           add back the wave function parts that are subtracted
*           this might be inefficient, but should have a clear logic
            call daxpy_(ndens,Half*WF_Ratio,D0(1,1),1,D0(1,2),1)
            call daxpy_(ndens,WF_Ratio,D1ao,1,D0(1,2),1)
*           scale the pdft part
            call dscal_(ndens,PDFT_Ratio,D0(1,5),1)
          end if
            Call mma_deallocate(D1ao)
          else If (Method.eq.'MSPDFT  ') Then
            Call Get_DArray('MSPDFTD5        ',D0(1,5),nDens)
            Call Get_DArray('MSPDFTD6        ',D0(1,6),nDens)
            call daxpy_(ndens,-1.0d0,D0(1,5),1,D0(1,2),1)
          end if

!          call dcopy_(ndens*5,0.0d0,0,D0,1)
!          call dcopy_(nG2,0.0d0,0,G2,1)


!************************
           !Call dscal_(Length,0.5d0,D0(1,4),1)
           !Call dscal_(Length,0.0d0,D0(1,4),1)

         RlxLbl='DLAO    '
!         Call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0(1,4))
! DMRG with the reduced AS
           !if(doDMRG)then
           !  length=ndim1  !yma
           !end if
         End If
         If (iPrint.ge.99) Call TriPrt(' G2',' ',G2(1,1),nG1)
*
*...  Close 'RELAX' file
1000     Continue
*
*...  Epilogue, end
#ifdef _CD_TIMING_
      Call CWTIME(PreppCPU2,PreppWall2)
      Prepp_CPU  = PreppCPU2 - PreppCPU1
      Prepp_Wall = PreppWall2 - PreppWall1
#endif

      Return
      End

      Subroutine Get_D1I(CMO,D1It,D1I,nish,nbas,nsym)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Dimension CMO(*), D1I(*),D1It(*)
      Integer nbas(nsym),nish(nsym)

      iOff1 = 0
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nIsh(iSym)
        If ( iBas.ne.0 ) then
          iOff2 = iOff1
          Do i = 1,iBas
            Do j = 1,iBas
              Sum = Zero
              Do k = 0,iOrb-1
                Sum = Sum + Two * CMO(iOff1+k*iBas+i)
     &                          * CMO(iOff1+k*iBas+j)
              End Do
              D1I(iOff2 + j) = Sum
            End Do
            iOff2 = iOff2 + iBas
          End Do
          iOff1 = iOff1 + iBas*iBas
        End If
      End Do
      Call Fold2(nsym,nBas,D1I,D1It)
      Return
      End

      Subroutine Get_D1A(CMO,D1A_MO,D1A_AO,
     &                    nsym,nbas,nish,nash,ndens)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
      Dimension CMO(*) , D1A_MO(*) , D1A_AO(*)
      Integer nbas(nsym),nish(nsym),nash(nsym)
      Real*8, Allocatable:: Scr1(:), Tmp1(:,:), Tmp2(:,:)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)


      iOff2 = 1
      iOff3 = 1
      ii=0
      Call mma_allocate(Scr1,2*nDens,Label='Scr1')
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iAsh = nAsh(iSym)
        iIsh = nIsh(iSym)
        Call dCopy_(iBas*iBas,[Zero],0,Scr1(iOff3),1)
        If ( iAsh.ne.0 ) then
          Call mma_allocate(Tmp1,iAsh,iAsh,Label='Tmp1')
          Call mma_allocate(Tmp2,iBas,iAsh,Label='Tmp2')
          Do i=1,iAsh
           Do j=1,iAsh
            Tmp1(j,i)=D1A_MO(itri(i+ii,j+ii))
           End do
          End do
          ii=ii+iash
          Call DGEMM_('N','T',
     &                iBas,iAsh,iAsh,
     &                One,CMO(iOff2+iIsh*iBas),iBas,
     &                    Tmp1,iAsh,
     &               Zero,Tmp2,iBas)
          Call DGEMM_('N','T',
     &                iBas,iBas,iAsh,
     &                One,Tmp2,iBas,
     &                    CMO(iOff2+iIsh*iBas),iBas,
     &               Zero,Scr1(iOff3),iBas)
          Call mma_deallocate(Tmp2)
          Call mma_deallocate(Tmp1)
        End If
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + iBas*iBas
      End Do
      Call Fold2(nSym,nBas,Scr1,D1A_AO)
      Call mma_deallocate(Scr1)
      Return
      End


      Subroutine Fold2(nSym,nBas,A,B)

      Implicit Real*8 (A-H,O-Z)

      Dimension nBas(*) , A(*) , B(*)

      iOff1 = 0
      iOff2 = 0
      Do iSym = 1, nSym
        mBas = nBas(iSym)
        Do iBas= 1, mBas
          Do jBas = 1 , iBas-1
            B(iOff2+jBas) =   A(iOff1+jBas)
          End Do
          B(iOff2+iBas) =  A(iOff1+iBas)
          iOff1 = iOff1 + mBas
          iOff2 = iOff2 + iBas
        End Do
      End Do

      Return
      end
************ columbus interface ****************************************
*read table of contents for gamma file

        subroutine read_lgtoc(lgtoc,gtoc,n)
        integer n,lgtoc
        real*8 gtoc(n)
          read(lgtoc) gtoc
        return
         end
