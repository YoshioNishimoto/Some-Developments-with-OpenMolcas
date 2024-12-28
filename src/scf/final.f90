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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!               2016,2017, Roland Lindh                                *
!***********************************************************************
      SubRoutine Final()
      use SCF_Arrays, only: Dens, OneHam, Ovrlp, TwoHam, CMO, EOrb, FockAO, OccNo, KntE, MssVlc, Darwin
      use InfSCF, only: nD
      Implicit None
!
!
!---- Read remaining one-electron integrals
      Call R1IntB()

      Call Final_Internal(nD)
!
      Return

      Contains
      SubRoutine Final_Internal(nD)
!***********************************************************************
!                                                                      *
!     purpose: perform final calculations                              *
!                                                                      *
!***********************************************************************
#ifdef _EFP_
      use EFP_Module, only: EFP_Instance
      use EFP, only: EFP_Shutdown
#endif
#ifdef _HDF5_
      Use mh5, Only: mh5_put_dset
      use SCFWFn, only: wfn_mocoef, wfn_mocoef_a,  wfn_mocoef_b, wfn_occnum, wfn_occnum_a, wfn_occnum_b,                 &
                        wfn_orbene, wfn_orbene_a, wfn_orbene_b, wfn_tpidx, wfn_tpidx_a, wfn_tpidx_b
#endif
      Use Interfaces_SCF, Only: dOne_SCF
      use OFembed, only: Do_OFemb, FMaux, NDSD
#ifdef _FDE_
      use Embedding_Global, only: embPot, embWriteEsp
#endif
      use SpinAV, only: DSc
      use InfSCF, only: nBT, nDens, DMOMax, FMOMax, kIVO, MaxBas, nSym, KSDFT, EneV, Falcon, iPrint, NoProp, DSCF, &
                        TotCPU, nFld, iStatPrn, E1V, E2V, FThr, iPrForm, MaxBXO, Name, NamFld, nBas, nBB, nBO,           &
                        nDel, nFro, nIter, nIterP, nnB, nnO, nOcc, nOrb, TimFld
#ifdef _FDE_
      use InfSCF, only: nAtoms
#endif
      use Constants, only: Zero, One, Two
      use stdalloc, only: mma_allocate, mma_deallocate
      use Files, only: LuOut
      use AddCorr, only: Do_Addc, Do_Tw
      Implicit None
      Integer nD
!
#ifdef _EFP_
      Logical, External:: EFP_On
#endif
!
!---- Define local variable
      Integer iD, iRC, iOpt, iSymLb, iFock, jFock, iCMO, iVirt, jVirt, ij, iBas, jBas, iSym, kl, lk, iRef, jRef, iiOrb, iOrb,   &
              nOccMax, nOccMin, iWFType,  kBas, iFld

      Real*8 TCPU1, TCPU2, Dummy, TWall1, TWall2
      Logical FstItr
      Character(LEN=8) RlxLbl,Method
      Character(LEN=60) Fmt
      Character(LEN=128) OrbName
      Logical RF_On,Langevin_On,PCM_On
      Character(LEN=80) Note
      Character(LEN=8) What
      Integer IndType(7,8)
      Real*8, Dimension(:), Allocatable:: Temp, CMOn, Etan, Epsn
      Real*8, Dimension(:,:), Allocatable:: GVFck, Scrt1, Scrt2, DMat, EOr
#ifdef _HDF5_
#include "Molcas.fh"
      character(Len=1), allocatable :: typestring(:)
      Integer nSSh(mxSym), nZero(mxSym)
      Integer i
#endif
      Integer nFldP
      Dimension Dummy(1)
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
      Call CWTime(TCpu1,TWall1)
!
      What='COEI'

      Call SorbCMOs(CMO,nBB,nD,EOrb,OccNo,nnB,nBas,nOrb,nSym)
!
      Call Put_darray('SCF orbitals',CMO(1,1),nBB)
      If (nD.eq.2) Then
         Call Put_darray('SCF orbitals_ab',CMO(1,2),nBB)
      End If
!
      If (nIter(nIterP).le.0) Then
!
         FstItr=.TRUE.
         Call SCF_Energy(FstItr,E1V,E2V,EneV)
         call dcopy_(nBT*nD,Dens(1,1,1),1,Dens(1,1,nDens),1)
         call dcopy_(nBT*nD,TwoHam(1,1,1),1,TwoHam(1,1,nDens),1)
!
         DMOMax = Zero
         FMOMax = Zero
!
      Else
!
!------- Compute improved virtuals (if needed)
         If (kIvo.eq.1) then
            Do iD = 1, nD
               Call IvoGen(OneHam,nBT,CMO(1,iD),nBO,EOrb(1,iD),nnO,nOcc(1,iD))
            End Do
         End If
!
!--------------------------------------------------
!------- Write stuff for gradient calculation -----
!--------------------------------------------------
!
!...  Add Fock operator matrix
!
!     Note that we need an array with four additional
!     elements due to the use of WrOne, which adds
!     four elements for some auxiliary information!
!
      Call mma_allocate(Temp,nBT+4,Label='Temp')
      Call FZero(Temp,nBT+4)
!
      Do iD = 1, nD
         Call DCopy_(nBT,FockAO(1,iD),1,Temp,1)
         iRc=-1
         iOpt=0
         RlxLbl='Fock Op '
         If (iD.eq.2) RlxLbl='Fock Op2'
         iSymLb=1
         Call WrOne(iRc,iOpt,RlxLbl,1,Temp,iSymLb)
         If ( iRc.ne.0 ) Then
            Write (6,*) 'Final: Error writing on ONEINT'
            Write (6,'(A,A)') 'RlxLbl=',RlxLbl
            Call Abend()
         End If
      End Do
      Call mma_deallocate(Temp)
!
!------- Construct the generalized variational Fock matrix
!
         Call mma_allocate(GVFck,nBT,nD,Label='GVFck')
         Call mma_allocate(Scrt1,MaxBas**2,nD,Label='Scrt1')
         Call mma_allocate(Scrt2,MaxBxO,nD,Label='Scrt2')
!
         iFock = 1
         jFock = 1
         iCMo  = 1
         Do iSym = 1, nSym
            If (nOrb(iSym).le.0) Cycle
!
            Do iD = 1, nD
               Call Square(FockAO(jFock,iD),Scrt1(1,iD),1,nBas(iSym),nBas(iSym))
!----------    Transform to MO basis
               Call DGEMM_('T','N',nOrb(iSym),nBas(iSym),nBas(iSym),  &
                           One,CMO(iCMo,iD),nBas(iSym),             &
                                 Scrt1(1,iD),nBas(iSym),              &
                           Zero,Scrt2(1,iD),nOrb(iSym))
               Call DGEMM_('N','N',nOrb(iSym),nOrb(iSym),nBas(iSym),  &
                           One,Scrt2(1,iD),nOrb(iSym),              &
                                 CMO(iCMo,iD),nBas(iSym),             &
                           Zero,Scrt1(1,iD),nOrb(iSym))
!
!----------    Set elements with both indices virtual to zero
               Do iVirt = 1, nOrb(iSym)-nOcc(iSym,iD)
                  jVirt = 1 + nOcc(iSym,iD)*nOrb(iSym) + nOcc(iSym,iD) + (iVirt-1)*nOrb(iSym)
                  call dcopy_(nOrb(iSym)-nOcc(iSym,iD),[Zero],0,Scrt1(jVirt,iD),1)
               End Do
!----------    Now project back to the SO basis
               Call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nOrb(iSym), &
                           One,CMO(iCMo,iD),nBas(iSym),            &
                                 Scrt1(1,iD),nOrb(iSym),             &
                           Zero,Scrt2(1,iD),nBas(iSym))
               Call DGEMM_('N','T',nBas(iSym),nBas(iSym),nOrb(iSym), &
                           One,Scrt2(1,iD),nBas(iSym),             &
                                 CMO(iCMo,iD),nBas(iSym),            &
                           Zero,Scrt1(1,iD),nBas(iSym))

               ij = iFock
               Do iBas = 1, nBas(iSym)
                  Do jBas = 1, iBas-1
                     kl =  nBas(iSym)*(jBas-1) + iBas
                     lk =  nBas(iSym)*(iBas-1) + jBas
                     GVFck(ij,iD) = Scrt1(kl,iD) + Scrt1(lk,iD)
                     ij = ij + 1
                  End Do
                  kl =  nBas(iSym)*(iBas-1) + iBas
                  GVFck(ij,iD) = Scrt1(kl,iD)
                  ij = ij + 1
               End Do
!
            End Do ! iD
!
            iFock = iFock + nBas(iSym)*(nBas(iSym)+1)/2
            jFock = jFock + nBas(iSym)*(nBas(iSym)+1)/2
            iCMo  = iCMo  + nBas(iSym)*nOrb(iSym)

         End Do  ! iSym
         Call mma_deallocate(Scrt2)
         Call mma_deallocate(Scrt1)
!
!...  Add elementary info
         Method='RHF-SCF '
         If (nD==2) Method='UHF-SCF '

         If (kIvo.ne.0) Method='IVO-SCF '
         If (KSDFT.ne.'SCF') Method='KS-DFT  '
         Call Put_cArray('Relax Method',Method,8)
!        Call Put_Energy(EneV)
         Call Store_Energies(1,[EneV],1)
         Call Put_dScalar('SCF energy',EneV)
!        If (nD==2) Call Put_dScalar('Ener_ab',EneV_ab)
         Call Put_iArray('nIsh',nOcc(1,1),nSym)
         If (nD==2) Call Put_iArray('nIsh_ab',nOcc(1,2),nSym)
         Call Put_iArray('nOrb',nOrb,nSym)
         Call Put_iArray('nDel',nDel,nSym)
         Call Put_iArray('nFroPT',nFro,nSym)    ! for Cholesky-CC
         Call Put_iArray('nDelPT',nDel,nSym)    !
!
!...  Add MO-coefficients
         Call Put_dArray('Last orbitals',CMO(1,1),nBB)
         If (nD==2) Call Put_dArray('CMO_ab',CMO(1,2),nBB)
!
!...  Add one body density matrix in AO/SO basis
!
!     If the density matrix is explicitly calculated in the beginning we will
!     use that one instead.
         If (nIter(nIterP).le.0) Then
            Call DensAB(nBT,nDens,nD,Dens)
         Else
            Call mma_allocate(DMat,nBT,nD,Label='DMat')
!
            if(nD==1) then
               Call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(1,1),nBB,OccNo(1,1),DMat(1,1),.true.)
            else
               Call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(1,1),nBB,OccNo(1,1),DMat(1,1),.true.)
               Call dOne_SCF(nSym,nBas,nOrb,nFro,CMO(1,2),nBB,OccNo(1,2),DMat(1,2),.false.)
               Call Put_dArray('D1aoVar',Dummy,0) ! Undefined the field.
            End If

            Call DensAB(nBT,1,nD,DMat)
            Call mma_deallocate(DMat)
         End If

         If(nD.eq.2) Then
            Call mma_allocate(CMOn,nBB,Label='CMOn')
            Call mma_allocate(Etan,nnB,Label='Etan')
            Call mma_allocate(Epsn,nnB,Label='Epsn')
            Call NatoUHF(Dens(1,1,1),Dens(1,2,1),FockAO(1,1),FockAO(1,2),nBT,CMO(1,1),nBB,Ovrlp,   &
                         CMOn,Etan,Epsn,nnB,nSym,nBas,nOrb)
            Call PadCMO(CMOn,CMOn,nSym,nBas,nOrb)
            Call PadEor(Etan,Etan,nSym,nBas,nOrb)
            Call PadEor(Epsn,Epsn,nSym,nBas,nOrb)
         End If
!
!...  Add generalized Fock matrix
         If (nD==2) Then
            Call DaXpY_(nBT,One,GVFck(1,2),1,GVFck(1,1),1)
         Else
            Call DScal_(nBT,Two,GvFck(1,1),1)
         End If
         Call Put_dArray('FockOcc',GVFck(1,1),nBT)
         Call mma_deallocate(GVFck)
!
!... Add SCF orbital energies
         Call Put_dArray('OrbE',EOrb(1,1),nnB)
         if(nD==2) Call Put_dArray('OrbE_ab',EOrb(1,2),nnB)
!
         Call Put_dScalar('Thrs',Fthr)
      EndIf

#ifdef _FDE_
      ! Embedding
      if (embPot.and.(embWriteEsp)) then
         Call embPotOutput(nAtoms,Dens)
      end if
#endif
!
! t.t.;
!    store Fock matrix in Runfile in case of fragment calculation
      If (Falcon) Then
!                 write(6,*)'Fock matrix is written in RunFile.'
!                 write(6,*)'fck:'
!                 write(6,*)'nbt=',nbt
!                 write(6,*)'ndens=',ndens
!                 write(6,*) (FockAO(itt),itt=1,nbt)
         Call Put_dArray('Fragment_Fock',FockAO(1,1),nBT)
      End if

! t.t.; end
!---- Final printout
! original PrFin was splitted to 3 parts to run in UHF mode
!
      Do iSym=1,nSym
         nDel(iSym)=nBas(iSym)-nOrb(iSym)
      End Do
      Call PrFin0(Dens(1,1,nDens),Dens(1,2,nDens),nBT,EOrb(1,1),nnB,CMO(1,1),nBO,KntE)
!
      Do iD = 1, nD
         Call PrFin(OneHam,Ovrlp,Dens(1,iD,nDens),TwoHam(1,iD,nDens),nBT,EOrb(1,iD),OccNo(1,iD),nnB,CMO(1,iD),nBO,Note,   &
                    iD-1,MssVlc,Darwin)
!
         Call PrFin2(Ovrlp,nBT,OccNo(1,iD),nnB,CMO(1,iD),nBO,Note)
      End Do
!
      If (iPrint.ge.2.and.nD==2) Then
         Call PriMO('Natural orbitals',.true.,.true.,Zero,2.0d9,nSym,nBas,nOrb,Name,Epsn,Etan,CMOn,iPrForm)
      End If
!------- Calculate expectation values
      If (.not.NoProp) Then
         If (iPrint.ge.3) Then
            Write(6,'(/6X,A)') 'Expectation values of various operators'
         End If
         Call Prpt()
      End If
!
! make a fix for energies for deleted orbitals
      Call mma_allocate(EOr,nnB,nD,Label='EOr')
!
      iRef=1
      jRef=1
      Do iSym=1,nSym
         iiOrb=nOrb(iSym)-nDel(iSym)
         Do iOrb=1,iiOrb
            Do iD = 1, nD
               EOr(iRef,iD)=EOrb(jRef,iD)
            End Do
            iRef=iRef+1
            jRef=jRef+1
         End do
         Do iOrb=1,nDel(iSym)
            Do iD = 1, nD
               EOr(iRef,iD)=1000
            End Do
            iRef=iRef+1
         End do
      End Do

      If(nD==1) Then
         Do iSym=1,nSym
            IndType(1,iSym)=nFro(iSym)
            IndType(2,iSym)=nOcc(iSym,1)
            IndType(3,iSym)=0
            IndType(4,iSym)=0
            IndType(5,iSym)=0
            IndType(6,iSym)=nOrb(iSym)-nFro(iSym)-nOcc(iSym,1)-nDel(iSym)
            IndType(7,iSym)=nDel(iSym)
         End Do
      Else
         Do iSym=1,nSym
            nOccMax=Max(nOcc(iSym,1),nOcc(iSym,2))
            nOccMin=Min(nOcc(iSym,1),nOcc(iSym,2))
            IndType(1,iSym)=nFro(iSym)
            IndType(2,iSym)=nOccMin
            IndType(3,iSym)=0
            IndType(4,iSym)=nOccMax-nOccMin
            IndType(5,iSym)=0
            IndType(6,iSym)=nOrb(iSym)-nFro(iSym)-nOccMax-nDel(iSym)
            IndType(7,iSym)=nDel(iSym)
         End Do
      End If
      If(nD==1) then
         OrbName='SCFORB'
         If(KSDFT.EQ.'SCF') Then
            iWFtype=2
         Else
            Note=Trim(Note)//' / '//Trim(KSDFT)
            iWFtype=3
         End If
         Call WrVec_(OrbName,LuOut,What,nD-1,nSym,nBas,nBas,CMO(1,1),Dummy, OccNo(1,1),Dummy,EOr(1,1),  &
           Dummy,IndType,Note,iWFtype)
#ifdef _HDF5_
         nZero = 0
         call mma_allocate(typestring, nnB)
         Do i = 1, nSym
         nSSh(i) = nBas(i) - nFro(i) - nOcc(i,1) - nDel(i)
         End Do
         call orb2tpstr(NSYM,NBAS,NFRO,NOCC(1,1),NZERO,NZERO,NZERO,NSSH,NDEL,typestring)
         call mh5_put_dset(wfn_tpidx, typestring)
         call mma_deallocate(typestring)
         call mh5_put_dset(wfn_mocoef, CMO(1,1))
         call mh5_put_dset(wfn_occnum, OccNo(1,1))
         call mh5_put_dset(wfn_orbene, EOrb(1,1))
#endif
      Else
         OrbName='UHFORB'
         If(KSDFT.EQ.'SCF') Then
            iWFtype=4
         Else
            Note=Trim(Note)//' / '//Trim(KSDFT)
            iWFtype=5
         End If
         Call WrVec_(OrbName,LuOut,What,nD-1,nSym,nBas,nBas,CMO(1,1),CMO(1,2), OccNo(1,1),OccNo(1,2),      &
                     EOr(1,1),EOr(1,2), IndType, Note,iWFtype)
#ifdef _HDF5_
         nZero = 0
         call mma_allocate(typestring, nnB)
         Do i = 1, nSym
         nSSh(i) = nBas(i) - nFro(i) - nOcc(i,1) - nDel(i)
         End Do
         call orb2tpstr(NSYM,NBAS,NFRO,NOCC(1,1),NZERO,NZERO,NZERO,NSSH,NDEL,typestring)
         call mh5_put_dset(wfn_tpidx_a, typestring)
         Do i = 1, nSym
         nSSh(i) = nBas(i) - nFro(i) - nOcc(i,2) - nDel(i)
         End Do
         call orb2tpstr(NSYM,NBAS,NFRO,NOCC(1,2),NZERO,NZERO,NZERO,NSSH,NDEL,typestring)
         call mh5_put_dset(wfn_tpidx_b, typestring)
         call mh5_put_dset(wfn_mocoef_a, CMO(1,1))
         call mh5_put_dset(wfn_occnum_a, OccNo(1,1))
         call mh5_put_dset(wfn_orbene_a, EOrb(1,1))
         call mh5_put_dset(wfn_mocoef_b, CMO(1,2))
         call mh5_put_dset(wfn_occnum_b, OccNo(1,2))
         call mh5_put_dset(wfn_orbene_b, EOrb(1,2))
#endif
         iBas=0
         Do iSym=1,nSym
            IndType(1,iSym)=nFro(iSym)
            IndType(2,iSym)=0
            IndType(3,iSym)=0
            IndType(4,iSym)=0
            IndType(5,iSym)=0
            IndType(6,iSym)=nOrb(iSym)-nFro(iSym)-nDel(iSym)
            IndType(7,iSym)=nDel(iSym)
            Do kBas=1,nBas(iSym)
               iBas=iBas+1
               If(Etan(iBas).gt.1.99d0) Then
                  IndType(2,iSym)=IndType(2,iSym)+1
                  IndType(6,iSym)=IndType(6,iSym)-1
               Else If(Etan(iBas).gt.0.01d0) Then
                  IndType(4,iSym)=IndType(4,iSym)+1
                  IndType(6,iSym)=IndType(6,iSym)-1
               End If
            End Do
         End Do
         OrbName='UNAORB'
         Note='UHF natural orbitals'
         If(KSDFT.EQ.'SCF') Then
            iWFtype=6
         Else
            iWFtype=7
         End If
         Call WrVec_(OrbName,LuOut,What,0,nSym,nBas,nBas,CMOn,Dummy,Etan,Dummy,Epsn,Dummy,IndType, Note,iWFtype)
#ifdef _HDF5_
         call orb2tpstr(NSYM,NBAS,NFRO,IndType(2,:),NZERO,IndType(4,:),NZERO,IndType(6,:),NDEL,typestring)
         call mh5_put_dset(wfn_tpidx, typestring)
         call mma_deallocate(typestring)
         call mh5_put_dset(wfn_mocoef, CMOn)
         call mh5_put_dset(wfn_occnum, Etan)
         call mh5_put_dset(wfn_orbene, Epsn)
#endif
         Call mma_deallocate(Epsn)
         Call mma_deallocate(Etan)
         Call mma_deallocate(CMOn)
      End If

      Call mma_deallocate(EOr)
!
!---- release Buffers for semi-direct SCF
      If (DSCF) Call ClsBuf()
!---- release SEWARD
      If ( DSCF             .or.   &
           RF_On()          .or.   &
           Langevin_On()    .or.   &
           PCM_On()         .or.   &
           Do_OFemb         .or.   &
           Do_Tw            .or.   &
           Do_Addc          .or.   &
#ifdef _EFP_
           EFP_On()         .or.   &
#endif
           KSDFT.ne.'SCF'        ) Call ClsSew()
!
      If (Allocated(FMaux)) Call mma_deallocate(FMaux)
      If (Allocated(NDSD)) Call mma_deallocate(NDSD)
      If (Allocated(DSc)) Call mma_deallocate(DSc)
#ifdef _EFP_
      If (EFP_On()) Call EFP_ShutDown(EFP_Instance)
#endif
!
      Call CWTime(TCpu2,TWall2)
      TimFld(15) = TimFld(15) + (TCpu2 - TCpu1)
      TotCpu=Max(TCpu2,0.1D0)
!
!---- Write out timing informations
      If(.not.NoProp) Then
         nFldP=nFld-1
      Else
         nFldP=nFld-2
      End If
      Fmt='(2x,A)'
      if(iStatPRN.gt.0) then
      Write(6,*)
      Call CollapseOutput(1,'Statistics and timing')
      Write(6,'(3X,A)')     '---------------------'
      Write(6,*)
      Write(6,Fmt)'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,Fmt)'   Part of the program                              CPU    fraction'
      Write(6,Fmt)'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
      Do iFld = 1, nFldP
         If ((iFld.eq.11).or.(iFld.eq.12)) Cycle
         Write(6,'(2x,A45,2f10.2)')NamFld(iFld),TimFld(iFld),TimFld(iFld)/TotCpu
      End Do
      Write(6,*)
      Write(6,'(2x,A45,2F10.2)')NamFld(nFld),TotCpu
      Write(6,Fmt)'- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
      Call CollapseOutput(0,'Statistics and timing')
      Write(6,*)
      endif

      End SubRoutine Final_Internal

      End SubRoutine Final
