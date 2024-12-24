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
! Copyright (C) 2024, Yoshio Nishimoto                                 *
!***********************************************************************
! Some memo
! "RASSCF energy for state X" is D^SS*V(D^SS)
! State averaged energies, which are used for the optimization in RASSCF, are D^SS*V(D^SS)
!
! electron-electron contributions are symmetric: D1*V(D2) = D2*V(D1),
! whereas electron-nuclear contributions are not always: D1*V(N) /= D2*V(N),
! so D*V/2 + D*V/2 = D*V is valid only if the explicit density used for defining the (SS) energy
! is the same to the ASC density
! For D^SS*V(D^SA), the correct nuclear contribution is (D^SS + D^SA)/2 * V^N
!
! The definition of the PCM density during SCF:
! def_solv = 1     : state-specific
! def_solv = 3,4,5 : state-average
! State-averaged PCM is invariant wrt rotations between internal states,
! whereas state-specific PCM is not invariant, so the state rotation contributions have to be
! included during Z-vector.
! The state rotation contribution may be computed at the first cycle unless def_solv = 3 (only?).
Module PCM_grad

  use definitions, only: iwp,wp,u6
  use stdalloc, only: mma_allocate, mma_deallocate
  use DWSol, only: DWSol_init,DWSol_final,DWSol_wgt,DWSolv,W_SOLV

  Implicit None

  ! RF_On() = .true. is equivalent to lRF = .true.
  ! PCM is for PCM (must be combined with cond: C-PCM),
  ! so PCM can be false even if lRF = .true. (not considered, though)
  ! do_RF below is different: do_RF is whether PCM equation is solved in MCLR or not.
  ! If RFPERT is true, we use the fixed reaction field and do not solve the PCM equation.

  ! iStpPCM = 0 : before initization
  ! iStpPCM = 1 : after initialization, before Z-vector
  ! iStpPCM = 2 : during Z-vector
  ! iStpPCM = 3 : after Z-vector
  integer(kind=iwp) :: iStpPCM = 0
  ! flag for PCM
  logical(kind=iwp) :: do_RF = .false.
! logical(kind=iwp), pointer :: do_PCM => do_RF
  ! definition of the solvation energy (WIP)
  ! def_solv = 1 : \Delta F = D^SS*V(D^SS)/2
  ! def_solv = 3 : \Delta F = D^SS * (V^N + V(D^SA)) - D^SA*V(D^SA)/2
  ! def_solv = 4 : \Delta F = D^SS * V(D^SA)/2 = D^SS * (V^N + V(D^SA)) - D^SS*V(D^SA)/2
  ! def_solv = 5 : \Delta F = D^SA * V(D^SA)/2
  integer(kind=iwp) :: def_solv = 0_iwp
  ! whether RFpert has been used in CASPT2
  logical(kind=iwp) :: PT2_solv = .false.
  ! RFPERT
  logical(kind=iwp) :: RFPERT = .false.

  !!!  One-electron integral contirbutions (in AO)
  ! from the SCF density
  ! PCMSCFAO(:,1) :: nuclear + electron contributions
  ! PCMSCFAO(:,2) :: nuclear contributions (1st index of DrvXV)
  ! PCMSCFAO(:,3) :: electron contributions (2nd index of DrvXV)
  real(kind=wp), Allocatable :: PCMSCFAO(:,:)
  ! from the state-specific (RLXROOT) density
  ! The meaning of the second index is the same above
  real(kind=wp), Allocatable :: PCMSSAO(:,:)
  ! PCMSSAO and PCMSSOori can be different if internal state rotations have to be considered,
  ! that is not single-state, def_solv = 3, or CASPT2 (other than XMS)
  ! PCMSSAO will contain rotated quantities, whereas PCMSSAOori will contain state-specific quantities
  real(kind=wp), Allocatable :: PCMSSAOori(:,:)
  ! from the inactive density
  real(kind=wp), Allocatable :: PCMIAO(:,:)

  !!!  One-electron integral contirbutions (in MO)
  ! from the SCF density
  ! The meaning of the second index is the same to AO
  real(kind=wp), Allocatable :: PCMSCFMO(:,:)
  ! from the state-specific (RLXROOT) density
  real(kind=wp), Allocatable :: PCMSSMO(:,:)
  ! from the state-specific (RLXROOT) density
  real(kind=wp), Allocatable :: PCMSSMOori(:,:)
  ! from the unrelaxed PT2 density
  real(kind=wp), Allocatable :: PCMPT2MO(:,:)

  !!! Square active density (in MO)
  ! used for PCM in SCF
  real(kind=wp), Allocatable, target :: DSCFMO(:,:)
  ! state-specific (RLXROOT) density
  real(kind=wp), Allocatable, target :: DSSMO(:,:)
  ! state-specific (RLXROOT) density
  real(kind=wp), Allocatable, target :: DSSMOori(:,:)

  !!! Square active density (in AO)
  ! used for PCM in SCF
  real(kind=wp), Allocatable :: DSCFAO(:)
  ! state-specific (RLXROOT) density
  real(kind=wp), Allocatable :: DSSAO(:)
  ! state-specific (RLXROOT) density
  real(kind=wp), Allocatable :: DSSAOori(:)
  ! state-specific (RLXROOT) density
! real(kind=wp), Allocatable :: DIAO(:)

  !!! Response densities and integrals
  real(kind=wp), Allocatable :: DZMO(:)      ! square density in MO (all orbitals)
  real(kind=wp), Allocatable :: DZACTMO(:,:) ! active response density in MO
  real(kind=wp), Allocatable :: DZAO(:)      ! square density in AO
  real(kind=wp), Allocatable :: PCMZMO(:,:)  ! square integral in MO
  real(kind=wp), Allocatable :: PCMZAO(:,:)  ! square integral in AO

  ! ISRot(:,:,1) :: internal rotations during Z-vector (Ap = ipS1)
  ! ISRot(:,:,2) :: trial vector (p = ipCId)
  ! ISRot(:,:,3) :: residue (ipST)
  ! ISRot(:,:,4) :: can be preconditioned ipST
  ! ISRot(:,:,5) :: solution (= ipCIT)
  ! ISRot(:,:,6) :: Pvec
  ! ISRot(:,:,7) :: Qvec
  ! ISRot(:,:,8) :: Uvec
  ! ISRot(:,:,9) :: R0
! real(kind=wp), Allocatable :: ISRot(:,:,:)  ! internal state rotations
  real(kind=wp), Allocatable :: ISRot(:,:)  ! internal state rotations, used in CIDens_sa

  !! weight of solvent
! real(kind=wp), Allocatable :: W_SOLV(:)

  !! E_NN^PCM: nuclear-nuclear interaction energies (< 0)
  !! To be precise, interaction energies between nuclei (positive) and
  !! ASCs induced by nuclei (negative)
  real(kind=wp) :: potnuc_pcm = 0.0d+00

  integer(kind=iwp) :: iCharge_PCM
  integer(kind=iwp) :: IPCMROOT

  !!! Apparent surface charges (ASCs)
  !Real*8, Allocatable :: ASCSCF(:)

  logical(kind=iwp), pointer :: isNAC

  integer(kind=iwp), pointer :: iCharge_Ref
  integer(kind=iwp), pointer :: ipCI
  integer(kind=iwp), pointer :: iRlxRoot
  integer(kind=iwp), pointer :: iSpin
  integer(kind=iwp), pointer :: NSSA(:)
  integer(kind=iwp), pointer :: nBas(:)
  integer(kind=iwp), pointer :: nConf
  integer(kind=iwp), pointer :: ncsf(:)
  integer(kind=iwp), pointer :: nRoots
  integer(kind=iwp), pointer :: nSym
  integer(kind=iwp), pointer :: ntAsh
  integer(kind=iwp), pointer :: ntBsqr
  integer(kind=iwp), pointer :: ntBtri
  integer(kind=iwp), pointer :: State_Sym
  integer(kind=iwp), pointer :: LUJOB

  integer(kind=iwp), pointer :: nFro(:)
  integer(kind=iwp), pointer :: nIsh(:)
  integer(kind=iwp), pointer :: nAsh(:)
  integer(kind=iwp), pointer :: nA(:)
  integer(kind=iwp), pointer :: nOrb(:)
  integer(kind=iwp), pointer :: ipCM(:)
  integer(kind=iwp), pointer :: ipMat(:,:)
  integer(kind=iwp), pointer :: iToc(:)

  real(kind=wp), pointer :: xispsm(:,:)
  real(kind=wp), pointer :: ERASSCF(:)
  real(kind=wp), pointer :: weight(:)


contains
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_grad_init(isNAC_, &
                           iCharge_Ref_,ipCI_,iRlxRoot_,iSpin_,NSSA_,nConf_,ncsf_,nRoots_, nSym_, &
                           ntAsh_,ntBsqr_,ntBtri_,State_Sym_, LUJOB_, &
                           nFro_,nIsh_,nAsh_,nA_,nBas_,nOrb_,ipCM_,ipMat_,iToc_, &
                           xispsm_,ERASSCF_,weight_)

  use ISRotation, only: InvEne,InvSCF
  use rctfld_module, only: RSlPar

  implicit none

#include "detdim.fh"

  logical(kind=iwp), intent(in), target :: isNAC_
  integer(kind=iwp), intent(in), target :: iCharge_Ref_,ipCI_, iRlxRoot_, iSpin_, NSSA_(*), nConf_, ncsf_(*), &
    nRoots_, nSym_, ntAsh_, ntBsqr_, ntBtri_, State_Sym_, LUJOB_, &
    nFro_(*), nIsh_(*), nAsh_(*), nA_(*), nBas_(*), nOrb_(*), ipCM_(*), ipMat_(*), iToc_(*)
  real(kind=wp), target :: xispsm_(*),ERASSCF_(*),weight_(*)

  integer(kind=iwp) :: i

  isNAC          => isNAC_

  iCharge_Ref    => iCharge_Ref_
  ipCI           => ipCI_
  iRlxRoot       => iRlxRoot_
  iSpin          => iSpin_
  NSSA(1:2)      => NSSA_(1:2)
  nConf          => nConf_
  ncsf(1:8)      => ncsf_(1:8)
  nRoots         => nRoots_
  nSym           => nSym_
  ntAsh          => ntAsh_
  ntBsqr         => ntBsqr_
  ntBtri         => ntBtri_
  State_Sym      => State_Sym_
  LUJOB          => LUJOB_

  nFro(1:8)      => nFro_(1:8)
  nIsh(1:8)      => nIsh_(1:8)
  nAsh(1:8)      => nAsh_(1:8)
  nA(1:8)        => nA_(1:8)
  nBas(1:8)      => nBas_(1:8)
  nOrb(1:8)      => nOrb_(1:8)
  ipCM(1:8)      => ipCM_(1:8)
  ipMat(1:8,1:8) => ipMat_(1:64)
! write (*,*) "ipCM"
! write (*,'(8i4)') ipcm(1:8)
! write (*,*) "ipMat"
! do i = 1, 8
! write (*,'(8i4)') ipmat(i,1:8)
! end do
  iToc(1:3)      => iToc_(1:3) !! 3 is sufficient for the moment

  xispsm(1:MXPCSM,1:MXPICI) => xispsm_(1:MXPCSM*MXPICI)
  ERASSCF(1:nRoots) => ERASSCF_(1:nRoots)
  weight(1:nRoots) => weight_(1:nRoots)

  call mma_allocate(PCMSCFAO,ntBsqr,3,Label='PCMSCFAO')
  call mma_allocate(PCMSSAO,ntBsqr,3,Label='PCMSSAO')
  call mma_allocate(PCMSSAOori,ntBsqr,3,Label='PCMSSAOori')

  call mma_allocate(PCMSCFMO,ntBsqr,3,Label='PCMSCFMO')
  call mma_allocate(PCMSSMO,ntBsqr,3,Label='PCMSSMO')
  call mma_allocate(PCMSSMOori,ntBsqr,3,Label='PCMSSMOori')
  call mma_allocate(PCMPT2MO,ntBsqr,3,Label='PCMPT2MO')

  call mma_allocate(DSCFMO,ntAsh,ntAsh,Label='DSCFMO')
  call mma_allocate(DSSMO,ntAsh,ntAsh,Label='DSSMO')
  call mma_allocate(DSSMOori,ntAsh,ntAsh,Label='DSSMOori')

  call mma_allocate(DSCFAO,ntBsqr,Label='DSCFAO')
  call mma_allocate(DSSAO,ntBsqr,Label='DSSAO')
  call mma_allocate(DSSAOori,ntBsqr,Label='DSSAOori')

  call mma_allocate(DZMO,ntBsqr,Label='DZMO')
  call mma_allocate(DZACTMO,ntAsh,ntAsh,Label='DZACTMO')
  call mma_allocate(DZAO,ntBsqr,Label='DZAO')
  call mma_allocate(PCMZMO,ntBsqr,3,Label='PCMZMO')
  call mma_allocate(PCMZAO,ntBsqr,3,Label='PCMZAO')

! call mma_allocate(ISRot,nRoots,nRoots,9,Label='ISRot')
  call mma_allocate(ISRot,nRoots,nRoots,Label='ISRot')

! call mma_allocate(W_SOLV,nRoots,Label='W_SOLV')

  iStpPCM = 1 !! after initialization, before Z-vector
  do_RF = .true.
  if (RFPERT) do_RF = .false.

  if (def_solv==0_iwp) then
!   call Get_iScalar('PCMRoot',iPCMRoot)
    call Get_iScalar('RF CASSCF root',iPCMRoot)
    if (iPCMRoot>  0) def_solv = 1_iwp
    if (iPCMRoot== 0) def_solv = 3_iwp
    if (iPCMRoot==-1) def_solv = 4_iwp
    if (iPCMRoot==-2) def_solv = 5_iwp
    if (iPCMRoot==-3) def_solv = 6_iwp

    if (def_solv==0_iwp) then
      write (u6,'(" iPCMRoot = ",i3," is not considered yet")')
      call abend()
    end if
  end if
! write (u6,*) "PCMRoot  = ", iPCMRoot
! write (u6,*) "def_solv = ", def_solv

  !! For weighted solvation
  !! assume that noneq is irrelevant in MCLR
  call DWSol_init(IPCMROOT,nRoots,.false.)
  call DWSol_wgt(2,ERASSCF)
  write (u6,*) "weight of the solvation density"
  do i = 1, nroots
    write (u6,*) "weight:",i,W_solv(i)
  end do

! call dcopy_(9*nRoots*nRoots,[0.0d+00],0,ISRot,1)
  call dcopy_(nRoots*nRoots,[0.0d+00],0,ISRot,1)

  if (def_solv==1) then
    InvEne    = .false.
    InvSCF    = .false.
  end if
  if (def_solv==3) then
    InvEne    = .true.
    InvSCF    = .true.
!   InvSCF    = .false.
  end if
  if (def_solv==4) then
    InvEne    = .false.
    InvSCF    = .true.
  end if

  if (DWSolv%DWZeta /= 0.0d+00) InvSCF = .false.

  End Subroutine PCM_grad_init
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_grad_final()

  implicit none

  !! nullify the pointer

  if (allocated(PCMSCFAO)) call mma_deallocate(PCMSCFAO)
  if (allocated(PCMSSAO))  call mma_deallocate(PCMSSAO)
  if (allocated(PCMSSAOori)) call mma_deallocate(PCMSSAOori)

  if (allocated(PCMSCFMO)) call mma_deallocate(PCMSCFMO)
  if (allocated(PCMSSMO))  call mma_deallocate(PCMSSMO)
  if (allocated(PCMSSMOori)) call mma_deallocate(PCMSSMOori)
  if (allocated(PCMPT2MO)) call mma_deallocate(PCMPT2MO)

  if (allocated(DSCFMO))   call mma_deallocate(DSCFMO)
  if (allocated(DSSMO))    call mma_deallocate(DSSMO)
  if (allocated(DSSMOori)) call mma_deallocate(DSSMOori)

  if (allocated(DSCFAO))   call mma_deallocate(DSCFAO)
  if (allocated(DSSAO))    call mma_deallocate(DSSAO)
  if (allocated(DSSAOori)) call mma_deallocate(DSSAOori)

  if (allocated(DZMO))     call mma_deallocate(DZMO)
  if (allocated(DZACTMO))  call mma_deallocate(DZACTMO)
  if (allocated(DZAO))     call mma_deallocate(DZAO)
  if (allocated(PCMZMO))   call mma_deallocate(PCMZMO)
  if (allocated(PCMZAO))   call mma_deallocate(PCMZAO)

  if (allocated(ISRot))    call mma_deallocate(ISRot)

! if (allocated(W_SOLV))   call mma_deallocate(W_SOLV)
  call DWSol_final()

  End Subroutine PCM_grad_final
!
!-----------------------------------------------------------------------
!
  Subroutine PrepPCM()

  use Arrays, only: CMO
!
! Prepare PCM-related integrals etc
!
  implicit none

  real(kind=wp) :: PotNuc,rdum(1)
  real(kind=wp), allocatable :: Htmp(:),Gtmp(:),D1ao(:)
  integer(kind=iwp) :: leng,iSym,ip1
  logical(kind=iwp) :: NonEq,First,Dff,Do_DFT,lRF

  leng = ntBtri
  Call mma_allocate(Htmp,leng,Label='Htmp')
  Call mma_allocate(Gtmp,leng,Label='Gtmp')
  Call mma_allocate(D1ao,ntBsqr,Label='D1ao')

  !! Compute the density used for solvation during SCF
  call PCM_grad_dens(1) ! SCF
  call PCM_grad_dens2(1,DSCFMO,DSCFAO)
  call fold(nSym,nBas,DSCFAO,D1ao)

! write (6,*) "DSCFAO"
! do isym = 1, leng
!   write (6,'(i3,f20.10)') isym,d1ao(isym)
! end do

  NonEq  = .False.
  First  = .True.
  Dff    = .False.
  Do_DFT = .True.
  lRF    = .True.
  PotNuc = 0.0d+00
  Call Get_dScalar('PotNuc',PotNuc)

  !! Htmp: nuclear-electron contributions
  !! Gtmp: electron-electron contributions
  Htmp(:) = 0.0d+00
  Gtmp(:) = 0.0d+00
!     write (*,*) "icharge = ", icharge_PCM
! write (6,*) "doing SA density"
  Call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF, &
             'SCF',0.0d+00,iCharge_PCM,iSpin,rdum,rdum,0,'1234',Do_DFT)
!           do i = 1, leng
!           write (*,'(i3,3f20.10)')i,htmp(i),gtmp(i),d1ao(i)
!           end do

  ip1 = 1
  do iSym = 1, nSym
    call square(Htmp(ip1),PCMSCFAO(ipCM(iSym),2),1,nBas(iSym),nBas(iSym))
    call square(Gtmp(ip1),PCMSCFAO(ipCM(iSym),3),1,nBas(iSym),nBas(iSym))
    ip1 = ip1 + nBas(iSym)*(nBas(iSym)+1)/2
  end do
  call dcopy_(ntBsqr,PCMSCFAO(1,2),1,PCMSCFAO(1,1),1)
  call daxpy_(ntBsqr,1.0d+00,PCMSCFAO(1,3),1,PCMSCFAO(1,1),1)

  PCMSCFMO = PCMSCFAO
! write (*,*) "inside PrepPCM PCMSCFAO"
! call sqprt(pcmscfao(1,1),nbas(1))
! call sqprt(pcmscfao(1,2),nbas(1))
! call sqprt(pcmscfao(1,3),nbas(1))
  call tcmo(PCMSCFMO(1,1),1,1)
  call tcmo(PCMSCFMO(1,2),1,1)
  call tcmo(PCMSCFMO(1,3),1,1)
! write (*,*) "inside PrepPCM PCMSCFMO"
! call sqprt(pcmscfmo(1,1),nbas(1))
! call sqprt(pcmscfmo(1,2),nbas(1))
! call sqprt(pcmscfmo(1,3),nbas(1))

  !! Compute SS density
  call PCM_grad_dens(2) ! SS
! call sqprt(dssmo,ntash)
  if (isNAC) then
    call PCM_grad_dens2(2,DSSMO,DSSAO)
  else
    call PCM_grad_dens2(1,DSSMO,DSSAO)
  end if
  call fold(nSym,nBas,DSSAO,D1ao)
! write (6,*) "DSSAO"
! do isym = 1, leng
!   write (6,'(i3,f20.10)') isym,d1ao(isym)
! end do

  Htmp(:) = 0.0d+00
  Gtmp(:) = 0.0d+00
! write (6,*) "doing SS density"
! do isym = 1, 10
! write (6,*) isym,d1ao(isym)
! end do
  Call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF, &
             'SCF',0.0d+00,iCharge_PCM,iSpin,rdum,rdum,0,'1234',Do_DFT)

  ip1 = 1
  do iSym = 1, nSym
    call square(Htmp(ip1),PCMSSAO(ipCM(iSym),2),1,nBas(iSym),nBas(iSym))
    call square(Gtmp(ip1),PCMSSAO(ipCM(iSym),3),1,nBas(iSym),nBas(iSym))
    ip1 = ip1 + nBas(iSym)*(nBas(iSym)+1)/2
  end do
  call dcopy_(ntBsqr,PCMSSAO(1,2),1,PCMSSAO(1,1),1)
  call daxpy_(ntBsqr,1.0d+00,PCMSSAO(1,3),1,PCMSSAO(1,1),1)

! write (*,*) "inside PrepPCM PCMSSAO"
! call sqprt(pcmssao(1,1),nbas(1))
! call sqprt(pcmssao(1,2),nbas(1))
! call sqprt(pcmssao(1,3),nbas(1))
  PCMSSMO = PCMSSAO
  call tcmo(PCMSSMO(1,1),1,1)
  call tcmo(PCMSSMO(1,2),1,1)
  call tcmo(PCMSSMO(1,3),1,1)

  if (def_solv /= 3) then
    call dcopy_(ntBsqr*3,PCMSSAO,1,PCMSSAOori,1)
    call dcopy_(ntBsqr*3,PCMSSMO,1,PCMSSMOori,1)
    call dcopy_(ntBsqr,DSSAO,1,DSSAOori,1)
    call dcopy_(ntAsh**2,DSSMO,1,DSSMOori,1)
  end if
! write (*,*) "inside PrepPCM PCMSSMO"
! call sqprt(pcmssmo(1,1),nbas(1))
! call sqprt(pcmssmo(1,2),nbas(1))
! call sqprt(pcmssmo(1,3),nbas(1))

  if (allocated(Htmp)) Call mma_deallocate(Htmp)
  if (allocated(Gtmp)) Call mma_deallocate(Gtmp)
  if (allocated(D1ao)) Call mma_deallocate(D1ao)

  End Subroutine PrepPCM
!
!-----------------------------------------------------------------------
!
  Subroutine PrepPCM2(mode,DMO,DAO,PCMAO,PCMMO)
!
! Prepare PCM-related integrals etc
! mode = 1: gradient
! mode = 2: non-adiabatic coupling
!
  use Arrays, only: CMO

  implicit none

  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: DMO(ntAsh**2)
  real(kind=wp), intent(inout) :: DAO(ntBsqr),PCMAO(ntBsqr,3),PCMMO(ntBsqr,3)

  real(kind=wp) :: PotNuc,rdum(1)
  real(kind=wp), allocatable :: Htmp(:),Gtmp(:),D1ao(:)
  integer(kind=iwp) :: leng,iSym,ip1!,i
  logical(kind=iwp) :: NonEq,First,Dff,Do_DFT,lRF

  leng = ntBtri
  Call mma_allocate(Htmp,leng,Label='Htmp')
  Call mma_allocate(Gtmp,leng,Label='Gtmp')
  Call mma_allocate(D1ao,leng,Label='D1ao')

  !! Compute SCF density
  call PCM_grad_dens2(mode,DMO,DAO)
  call fold(nSym,nBas,DAO,D1ao)

  NonEq  = .False.
  First  = .True.
  Dff    = .False.
  Do_DFT = .True.
  lRF    = .True.
  PotNuc = 0.0d+00
  Call Get_dScalar('PotNuc',PotNuc)

  !! Htmp: nuclear-electron contributions
  !! Gtmp: electron-electron contributions
  Htmp(:) = 0.0d+00
  Gtmp(:) = 0.0d+00
  Call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF, &
             'SCF',0.0d+00,iCharge_PCM,iSpin,rdum,rdum,0,'1234',Do_DFT)

  ip1 = 1
  do iSym = 1, nSym
    call square(Htmp(ip1),PCMAO(1,2),1,nBas(iSym),nBas(iSym))
    call square(Gtmp(ip1),PCMAO(1,3),1,nBas(iSym),nBas(iSym))
    ip1 = ip1 + nBas(iSym)*(nBas(iSym)+1)/2
  end do
  call dcopy_(ntBsqr,PCMAO(1,2),1,PCMAO(1,1),1)
  call daxpy_(ntBsqr,1.0d+00,PCMAO(1,3),1,PCMAO(1,1),1)

  PCMMO = PCMAO
  call tcmo(PCMMO(1,1),1,1)
  call tcmo(PCMMO(1,2),1,1)
  call tcmo(PCMMO(1,3),1,1)

  if (allocated(Htmp)) Call mma_deallocate(Htmp)
  if (allocated(Gtmp)) Call mma_deallocate(Gtmp)
  if (allocated(D1ao)) Call mma_deallocate(D1ao)

  return

  End Subroutine PrepPCM2
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_grad_dens(mode)
! Compute active density in MO
! mode=1: construct density used for PCM in SCF
! mode=2: mostly state-specific density

  use ipPage, only: W

  implicit none

  integer(kind=iwp), intent(in) :: mode

  real(kind=wp), allocatable:: CIL(:), CIR(:), G1r(:), G2r(:)
  real(kind=wp), pointer :: G1q(:,:) => null()
  integer(kind=iwp) :: nconfl,nconfr,nconf1,i,ij,j,jR,kR,ng1,ng2

  nconf1 = ncsf(state_sym)
  nConfL =Max(nconf1,int(xispsm(state_sym,1)))
  nConfR =nConfL

  ng1 = ntash*ntash
  ng2 = ng1*(ng1+1)/2

  call mma_allocate(CIL,nConfL,Label='CIL')
  call mma_allocate(CIR,nConfR,Label='CIR')
  Call mma_allocate(G1r,ng1,Label='G1r')
  Call mma_allocate(G2r,ng2,Label='G2r') ! not used?

  if (mode==1) G1q(1:ntAsh,1:ntAsh) => DSCFMO(:,1:ntAsh)
  if (mode==2) G1q(1:ntAsh,1:ntAsh) => DSSMO(:,1:ntAsh)
  G1r = 0.0d+00
  G1q = 0.0d+00

! write (*,*) "PCM_grad_dens with mode = ", mode
  if (mode==1) then
    Do jR = 1, nRoots
    ! if (mode==2 .and. jR /= iRlxRoot) cycle
      if (W_SOLV(jR).le.1.0d-10) cycle
      Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,state_sym)
      Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIR,state_sym)
      G1r = 0.0d+00
      Call Densi2(1,G1r,G2r,CIL,CIR,0,0,0,ng1,ng2)
!     write (*,*) "jR = ", jR
!     call sqprt(g1r,nash(1))
      !! For RDM1
      call daxpy_(ntAsh*ntAsh,W_SOLV(jR),G1r,1,G1q,1)
!     ij=0
!     Do i=0,ntAsh-1
!       Do j=0,i-1
!         ij=ij+1
!         G1q(ij)=G1q(ij)+ (G1r(1+i*ntAsh+j)+G1r(1+j*ntAsh+i))*0.5d+00
!       End Do
!       ij=ij+1
!       G1q(ij)=G1q(ij) + G1r(1+i*ntAsh+i)
!     End Do
    End Do
  else if (mode==2) then
    if (isNAC) then
      jR = NSSA(2)
      kR = NSSA(1)
    else
      jR = iRlxRoot
      kR = iRlxRoot
    end if
    Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,state_sym)
    Call CSF2SD(W(ipCI)%Vec(1+(kR-1)*nconf1),CIR,state_sym)
    Call Densi2(1,G1r,G2r,CIL,CIR,0,0,0,ng1,ng2)
    !! For RDM1
    call daxpy_(ntAsh*ntAsh,1.0D+00,G1r,1,G1q,1)
    do i = 1, ntAsh
      do j = 1, i-1
        G1q(i,j) = (G1q(i,j)+G1q(j,i))*0.5d+00
        G1q(j,i) = G1q(i,j)
      end do
    end do
  end if
  call mma_deallocate(CIL)
  call mma_deallocate(CIR)
  call mma_deallocate(G1r)
  call mma_deallocate(G2r)

  nConf=ncsf(state_sym)
  !! use weights
! if (mode==1) call dscal_(ntAsh**2,1.0d+00/dble(nRoots),G1q,1)

  Return

  End Subroutine PCM_grad_dens
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_grad_dens2(mode,DMO,DAO)
! Compute AO density using MO density (DMO)

  use ipPage, only: W
       use Arrays, only: Int1

  implicit none

  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in)  :: DMO(ntAsh,ntAsh)
  real(kind=wp), intent(out) :: DAO(:)
  real(kind=wp), allocatable:: WRK(:)

  integer(kind=iwp) :: i,iB,iiB,ijB,iOrb,iS,j,jB
  real(kind=wp) :: rd,rd2

! write (6,*) "PCM_grad_dens2, mode = ", mode
  call mma_allocate(WRK,ntBsqr,Label='WRK')

  !! transform MO to AO
  ! write (6,*) "ipmat(1,1) = ", ipmat(1,1) ! = 1
  WRK(:) = 0.0d+00
  Do iS=1,nSym
    iOrb = nOrb(iS)
    if (mode.eq.1) then
      !! For NAC, do not include the inactive density
      Do iB=1,nFro(is)+nIsh(is)
        WRK(ipmat(iS,iS)+iB-1+iOrb*(iB-1)) = 2.0D+00
      End Do
    end if

    Do iB=1,nAsh(iS)
      Do jB=1,nAsh(iS)
        iiB=nA(iS)+ib
        ijB=nA(iS)+jb
!       rd=DMO(iiB,ijB)
        rd=(DMO(iiB,ijB)+DMO(ijB,iiB))*0.5d+00
!       iij=iTri(iib,ijb)
        iiB=nFro(iS)+nIsh(iS)+ib
        ijB=nFro(iS)+nIsh(iS)+jb
        WRK(ipmat(iS,iS)+iiB-1+iOrb*(ijB-1)) = rd
      End Do
    End Do
  End Do
! write (6,*) "D in MO"
! call sqprt(WRK,norb(1))
  call tcmo(WRK,1,-2)
! write (6,*) "D in AO"
! call sqprt(WRK,norb(1))
  call dcopy_(ntBsqr,WRK,1,DAO,1)

! if (mode.eq.1) then
!   rd = 0.0d+00
!   rd2= 0.0d+00
!   write (6,*) "pcmscfao"
!   call sqprt(pcmscfao,12)
!   write (6,*) "dssao"
!   call sqprt(dssao,12)
!   do i = 1, 12*12
!   do j = 1, 12
!   ib = max(i,j)*(max(i,j)-1)/2
!   !  rd = rd + pcmscfao(i)*dssao(i)
!   rd = rd + pcmscfao(i,1)*dao(i)
!   rd2= rd2+ pcmssao(i,1)*dao(i)
!   rd = rd + pcmscfao(i+12*(j-1))*dssao(i+12*(j-1))
!   end do
!   end do
!   write (6,'("reference (PCM) en  = ",f20.10)') -112.14594735d+00
!   write (6,'("reference (PCM) en2 = ",f20.10)') -112.14600563d+00
!   write (6,'("PCM free energy SCF = ",f20.10)') rd*0.5d+00
!   write (6,'("PCM free energy SCF = ",f20.10)') rd*0.5d+00 + potnuc_pcm
!   write (6,'("PCM free energy SS  = ",f20.10)') rd*0.5d+00
!   write (6,'("PCM free energy SS  = ",f20.10)') rd2*0.5d+00 + potnuc_pcm
!   write (6,'("expected energy corr= ",f20.10)') -0.0002208582
!   write (6,'("expected energy corr= ",f20.10)') -0.0002791382
!     Call Put_dScalar('RF Self Energy',rd)
!   write (6,'("RF Self Energy      = ",f20.10)') rd
! end if

  call mma_deallocate(WRK)

  Return

  End Subroutine PCM_grad_dens2
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_grad_D2V(dens,vintMO,vintAO,first_,Dff_,NonEq_,Do_DFT_)

  use ipPage, only: W
  use rctfld_module, only: nTs

  implicit none

  real(kind=wp), intent(in)  :: dens(ntBsqr)
  real(kind=wp), intent(out) :: vintMO(ntBsqr,3),vintAO(ntBsqr,3)
  logical(kind=iwp), intent(in) :: first_, Dff_, NonEq_, Do_DFT_

  integer(kind=iwp) :: leng,iSym,ip1 ! ,i,iB,iiB,ijB,iOrb,iS,j,jB
  real(kind=wp) :: PotNuc,rdum(1)
  real(kind=wp), allocatable :: Htmp(:),Gtmp(:),D1ao(:)
  logical(kind=iwp) :: first, Dff, NonEq, Do_DFT, lRF

  real(kind=wp), allocatable :: PCM_Charges(:,:)
!
! density matrix (D; AO) --> electrostatic potential integral (V; AO)
!
  vintMO = 0.0d+00
  vintAO = 0.0d+00

  leng    = ntBtri
  First   = First_
  Dff     = Dff_
  NonEq   = NonEq_
  Do_DFT  = Do_DFT_
  lRF     = .True.
  PotNuc  = 0.0d+00 ! not used actually

  Call mma_allocate(Htmp,leng,Label='Htmp')
  Call mma_allocate(Gtmp,leng,Label='Gtmp')
  Call mma_allocate(D1ao,leng,Label='D1ao')

  call fold(nSym,nBas,dens,D1ao)
  Htmp(:) = 0.0d+00
  Gtmp(:) = 0.0d+00

  ! h1, vint, and dens is triangular in the AO basis
! Call DrvRF(Htmp,Gtmp,D1ao,RepNuc,nh1,First,Dff,NonEq,iCharge)
  Call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF, &
             'SCF',0.0d+00,iCharge_PCM,iSpin,rdum,rdum,0,'1234',Do_DFT)

  ip1 = 1
  do iSym = 1, nSym
    call square(Htmp(ip1),vintAO(ipCM(iSym),2),1,nBas(iSym),nBas(iSym))
    call square(Gtmp(ip1),vintAO(ipCM(iSym),3),1,nBas(iSym),nBas(iSym))
    ip1 = ip1 + nBas(iSym)*(nBas(iSym)+1)/2
  end do
  call dcopy_(ntBsqr,vintAO(1,2),1,vintAO(1,1),1)
  call daxpy_(ntBsqr,1.0d+00,vintAO(1,3),1,vintAO(1,1),1)

  vintMO = vintAO
! PCMSCFMO = PCMSCFAO
  call tcmo(vintMO(1,1),1,1)
  call tcmo(vintMO(1,2),1,1)
  call tcmo(vintMO(1,3),1,1)

  if (allocated(Htmp)) Call mma_deallocate(Htmp)
  if (allocated(Gtmp)) Call mma_deallocate(Gtmp)
  if (allocated(D1ao)) Call mma_deallocate(D1ao)

! if (iStpPCM==3) then
!   Call mma_allocate(PCM_Charges,2,nTS,Label='PCM_Charges')
!   Call Get_dArray('PCM Charges',PCM_Charges,2*nTs)
!   write (6,'("PCM_Charges induced by the relaxed density matrix")')
!   do iSym = 1, nTS
!     write (6,'(i4,f20.10)') iSym,PCM_Charges(1,iSym)+PCM_Charges(2,iSym)
!   end do
!   if (allocated(PCM_Charges)) Call mma_deallocate(PCM_Charges)
! end if

  Return

  End Subroutine PCM_grad_D2V
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_grad_TimesE2(idSym,rKappa,FockOut,ipCIOut)
  use ISRotation, only: ISR

  implicit none

  integer(kind=iwp), intent(in) :: idSym,ipCIOut
  real(kind=wp), intent(in) :: rKappa(*)
  real(kind=wp), intent(inout) :: FockOut(*)

  integer(kind=iwp) :: iS,jS,n1dens,n2dens,iA,jA,ip1,ip2
  real(kind=wp) :: rd,rdum(1)

  real(kind=wp),allocatable :: WRK(:)
  integer(kind=iwp) :: i,j

  Interface

  SubRoutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)
    Integer iispin, iCsym, iSSym
    Integer nInt1, nInt2s, nInt2a
    Real*8, Target:: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
    Integer ipCI1, ipCI2
    Logical Have_2_el
  End SubRoutine CISigma_sa

  End Interface
!
! Compute implicit contributions
! Explicit contributions are computed in existing subroutines
!
! The PCM contribution has a quartic dependence on the CI coefficient.
! The energy contribution to the S state is like (with state-averaged PCM):
!   E^S = sum_pq Dpq^S K_{pq,rs} D_rs^SA
!       = sum_pq sum_IJ <CI^S|Epq^IJ|CJ^S> K_{pq,rs} sum_T sum_rs w_T <CK^T|E_rs|CL^T>
! Second derivatives of the same density matrix can be evaluated by one-center-type
! integrals, whereas those of the different density matrix is evaluated by using
! response density
! The former is like: CId(I) = <CI|Epq^IJ|ZJ>*Hpq(D)
! The latter is like: CId(I) = <CI|Epq^IJ|CJ>*Hpq(Z)
! The former is evaluated through FIMO (or FAMO), but the latter is newly evaluated
!
! Different from the CI Lagrangian, TimesE2 is not dependent on the choice of the solvation energy,
! because Z-vector uses the diagonality of the Hamiltonian, but not the definition of the state-specific energy.
!

  !! Construct the orbital response density using SCF (state-averaged) density
  !! The Lagrangian (for generalized Brillouin theorem) is defined with the SCF density
  Call mma_allocate(WRK,ntBsqr,Label='WRK')
  Call OITD(rKappa,1,DZMO,WRK,.True.)
  Call mma_deallocate(WRK)

  !! Then, CI response density (DZACTMO), which has been computed in CIDens_sa
  !! through ci_kap (iStpPCM=2) or out_pt2 (iStpPCM=3)
! write (*,*) "dzactmo"
! call sqprt(dzactmo,nash(1))
  Do iS=1,nSym
    If (nAsh(iS).gt.0) Then
      jS=iEOr(is-1,iDSym-1)+1
      Do jA=1,nAsh(js)
!       rd=DZACTMO(iA+nA(iS),jA+nA(js))
!       ip1=nOrb(iS)*(nIsh(is)+iA-1)+ipCM(is)
        ip2=nIsh(iS) + nOrb(iS)*(nIsh(js)+jA-1)+ipmat(is,js)
        Call DaXpY_(nAsh(iS),1.0d+00,DZACTMO(1+nA(iS),jA+nA(jS)),1,DZMO(ip2),1)
      End Do
    End If
  End Do

  !! Transform to AO
  DZAO = DZMO
  Call TCMO(DZAO,1,-2)

  !! Compute V(D^Z) integrals
  !! integrals used for implicit contributions are just V(D^SA), so we already have them
  !! we only need PCMZMO(:,3) actually
  call PCM_grad_D2V(DZAO,PCMZMO,PCMZAO,.false.,.false.,.false.,.true.)

  !! For kap*Z (= RInt_generic + CI_KAP) contributions
  !! Similar to fockgen, but only the implicit contribution in the eigenenergy has to be
  !! evaluated here and the evaluation of the erfx term is not needed
  DZMO = 0.0D+00
  Do iS=1,nSym
    If (nBas(iS).gt.0) Then
      jS=iEOr(is-1,iDSym-1)+1
      Do iA=1,nAsh(is)
        Do jA=1,nAsh(js)
          rd=DSCFMO(nA(is)+iA,nA(js)+jA) ! solvent density during SCF
          ip1=nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)
          ip2=nBas(iS)*(nIsh(js)+jA-1) +ipMat(is,js)
          Call DaXpY_(nBas(iS),Rd,PCMZMO(ip1,3),1,DZMO(ip2),1)
        End Do
      End Do
    End If
  End Do

  If (iDsym.eq.1) Then
! If (state_sym.eq.1) Then !?
    Do iS=1,nSym
      If (nBas(iS)*nIsh(iS).gt.0) &
        Call DaXpY_(nBas(iS)*nIsh(is),2.0d+00,PCMZMO(ipMat(is,is),3),1,DZMO(ipMat(is,is)),1)
    End Do
  End If

  if (iStpPCM==2) then
    !! during Z-vector, compute kappa (DZAO)
    DZAO = 0.0D+00
    Do iS=1,nSym
      jS=iEOR(iS-1,idSym-1)+1
      If (nBas(is)*nBas(jS).ne.0) then
         Call DGeSub(DZMO(ipMat(iS,jS)),nBas(iS),'N',DZMO(ipMat(jS,iS)),nBas(jS),'T', &
                     DZAO(ipMat(iS,jS)),nBas(iS),nBas(iS),nBas(jS))
      end if
    End Do
    call daxpy_(ntBsqr,2.0d+00,DZAO,1,FockOut,1)
  else if (iStpPCM==3) then
    !! only for the connection term
    call daxpy_(ntBsqr,2.0d+00,DZMO,1,FockOut,1)
    !! No CI response
    return
  end if


  !! For ci*Z (= Kap_CI + Ci_Ci) contributions
  ! First, construct MO
! write (*,*) "pcmzmo"
! call sqprt(pcmzmo(1,3),norb(1))
  DZMO = 0.0d+00
  do iS = 1, nSym
    If (nBas(iS).gt.0) Then
      jS=iEOr(is-1,iDSym-1)+1
      Do jA=1,nAsh(js)
!       ip1=nIsh(iS)+1+nBas(jS)*(nIsh(js)+jA-1)!+ipCM(is)
        ip1=nIsh(iS)+1 + nBas(iS)*(nIsh(js)+jA-1)+ipMat(iS,iS)-1
        ip2=nIsh(iS)+1 + nBas(iS)*(nIsh(js)+jA-1)+ipMat(iS,jS)-1
        call daxpy_(nAsh(iS),+1.0d+00,PCMZMO(ip1,3),1,DZMO(ip2),1)
      end do
    end if
  end do
! write (*,*) "DZMO"
! call sqprt(dzmo,norb(1))

  !! Evaluate the CI derivative
  !! Note that ipCIOUT will be weighted with (2)W_SOLV
  call dswap_(nRoots,weight,1,W_SOLV,1)
  call dscal_(nRoots,2.0d+00,weight,1)
  Call CISigma_sa(0,state_sym,state_sym,DZMO,SIZE(DZMO),rdum,1,rdum,1,ipCI,ipCIOUT,.False.)
  call dscal_(nRoots,0.5d+00,weight,1)
  call dswap_(nRoots,weight,1,W_SOLV,1)

  !! consider the weight of the derivative
  if (DWSolv%DWZeta /= 0.0d+00) then
    call DWder(2,idsym,DZMO,SIZE(DZMO),DZMO,SIZE(DZMO),ISR%Ap, &
               ntash,nRoots,ERASSCF,itoc,LUJOB, &
               nSym,nAsh,nIsh,nBas,ipCM)
  end if

  End Subroutine PCM_grad_TimesE2
!
!-----------------------------------------------------------------------
!

  Subroutine PCM_grad_CLag(mode,ipCI,ipCID,SLag)

  use Arrays, only: FIMO, INT2
  use ipPage, only: W
  use ISRotation, only: InvSCF,ISR,ISR_RHS,InvEne,ISR_Projection

  implicit none

  integer(kind=iwp), intent(in) :: mode,ipCI,ipCID
  real(kind=wp), intent(inout) :: SLag(*)

  integer(kind=iwp) :: i,j,nconf1,idsym,is,js,ia,ja,jR,kR,ip1
  real(kind=wp) :: rdum,rtmp(1)
  real(kind=wp), allocatable :: SCFcont(:),SScont(:),SScori(:),vectmp(:),vectmpSS(:,:),vectmpSA(:,:)

  Interface

  SubRoutine CISigma_sa(iispin,iCsym,iSSym,Int1,nInt1,Int2s,nInt2s,Int2a,nInt2a,ipCI1,ipCI2,Have_2_el)
    Integer iispin, iCsym, iSSym
    Integer nInt1, nInt2s, nInt2a
    Real*8, Target:: Int1(nInt1), Int2s(nInt2s), Int2a(nInt2a)
    Integer ipCI1, ipCI2
    Logical Have_2_el
  End SubRoutine CISigma_sa

  End Interface

  nconf1 = ncsf(State_Sym)
  rtmp(1) = 0.0d+00
  if (isNAC) then
    jR = NSSA(2)
    kR = NSSA(1)
  else
    jR = iRlxRoot
    kR = iRlxRoot
  end if

  call mma_allocate(SCFcont,SIZE(FIMO),Label='SCFcont')
  call mma_allocate(SScont ,SIZE(FIMO),Label='SScont')
  call mma_allocate(SScori ,SIZE(FIMO),Label='SScori')
  call mma_allocate(vectmp,nconf1,Label='vectmp')
  call mma_allocate(vectmpSS,nconf1,nroots,Label='vectmpSS')
  call mma_allocate(vectmpSA,nconf1,nroots,Label='vectmpSA')
!
! CI derivative of the Lagrangian (dependent on the density that polarizes ASCs)
!
  if (PT2_solv .and. mode==1) then
    call daxpy_(ntBsqr,1.0d+00,PCMPT2MO(1,3),1,PCMSSMO(1,3),1)
  end if

  if ( mode == 2) then
    PCMSSMO(:,3) = PCMSSMO(:,3) - PCMSSMOori(:,3)
!   call dcopy_(nconf1*nRoots,W(ipCID)%Vec(1),1,vectmpSS,1)
!   call daxpy_(ntBsqr,-1.0d+00,PCMSSMOori(1,3),1,PCMSSMO(1,3),1)

!   call daxpy_(ntash**2,-1.0d+00,dssmoori,1,dssmo,1)
!   write (6,*) "active density due to slag rotation"
!   call sqprt(dssmo,ntash)
!   call daxpy_(ntash**2,+1.0d+00,dssmoori,1,dssmo,1)
!   write (6,*) "PCMSSMO due to slag rotation"
!   call sqprt(pcmssmo(1,3),nbas(1))
  end if

!
! CI derivative of the state-specific energy that is dependent on the definition of the solvation
!
  SCFcont = 0.0d+00 !! contributions to the density that polarizes ASCs during SCF
  SScont  = 0.0d+00 !! to rotated (due to non-iterative internal rotations) state-specific density
  SScori  = 0.0d+00 !! to original (unrotated) state-specific density

  do iS = 1, nSym
    If (nAsh(iS).gt.0) Then
      jS=iEOr(is-1,state_sym-1)+1
      Do jA=1,nAsh(js)
!       ip1=nIsh(iS)+1+nBas(jS)*(nIsh(js)+jA-1)!+ipCM(is)
        ip1=nIsh(iS)+1 + nBas(iS)*(nIsh(js)+jA-1) + ipMat(iS,jS)-1 !+ipCM(is)
        if (def_solv == 1) then
          !! subtract implicit D^SS*V(e,SCF)/2
          call daxpy_(nAsh(iS),+1.0d+00,PCMSSMO(ip1,3),1,SCFcont(ip1),1)
          !! explicit and implicit D^SS*V(e,SCF)/2
          if (.not.isNAC .and. mode==1) then
            call daxpy_(nAsh(iS),-1.0d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
          end if
          !! sum of these contributions are usually (RlxRoot = PCMRoot) zero,
          !! because D^SS is used for polarizing ASCs
        else if (def_solv == 2) then
          !! subtract explicit D^SS*(V(e,SA)/2
          call daxpy_(nAsh(iS),-0.5d+00,PCMSCFMO(ip1,1),1,SScont(ip1),1)
          !! implicit D^SS*V(e,SA)/2
          call daxpy_(nAsh(iS),+0.5d+00,PCMSSMO(ip1,1),1,SCFcont(ip1),1)
        else if (def_solv == 3) then
          !! implicit V(e) in FIMO
          call daxpy_(nAsh(iS),+1.0d+00,PCMSSMO(ip1,3),1,SCFcont(ip1),1)
          !! explicit + implicit erfx
          if (.not.isNAC .and. mode==1) then
            call daxpy_(nAsh(iS),-1.0d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
          end if
        else if (def_solv == 4) then
          !! implicit D^SS*V(e,SA) in FIMO
          call daxpy_(nAsh(iS),+1.0d+00,PCMSSMO(ip1,3),1,SCFcont(ip1),1)
          if (.not.isNAC .and. mode == 1) then
            !! explicit + implicit erfx
            call daxpy_(nAsh(iS),-1.0d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
            !! explicit correction term (1)
            call daxpy_(nAsh(iS),+0.5d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
            !! implicit correction term
            call daxpy_(nAsh(iS),+0.5d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
            call daxpy_(nAsh(iS),-0.5d+00,PCMSSMOori(ip1,3),1,SCFcont(ip1),1)
            !! explicit correction term (1)
            call daxpy_(nAsh(iS),-0.5d+00,PCMSCFMO(ip1,3),1,SScori(ip1),1)
          end if
        else if (def_solv == 5) then
          !! implicit D^SS*V(e,SA) in FIMO
!         call daxpy_(nAsh(iS),+1.0d+00,PCMSSMO(ip1,3),1,SCFcont(ip1),1)
!         !! explicit + implicit erfx
!         call daxpy_(nAsh(iS),-1.0d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
!         !! explicit correction term (1)
!         call daxpy_(nAsh(iS),+1.0d+00,PCMSCFMO(ip1,1),1,SCFcont(ip1),1)
!         !! implicit correction term
!         call daxpy_(nAsh(iS),+1.0d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
!         call daxpy_(nAsh(iS),-1.0d+00,PCMSSMO(ip1,3),1,SCFcont(ip1),1)
!         !! explicit correction term (1)
!         call daxpy_(nAsh(iS),-1.0d+00,PCMSCFMO(ip1,1),1,SScont(ip1),1)
          !! subtract explicit FIMO
          call daxpy_(nAsh(iS),-1.0d+00,PCMSCFMO(ip1,1),1,SScont(ip1),1)
          !! explicit correction term (V^N)
          call daxpy_(nAsh(iS),+1.0d+00,PCMSCFMO(ip1,2),1,SCFcont(ip1),1)
          !! explicit + implicit correction term (V^e)
          call daxpy_(nAsh(iS),+1.0d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
        else if (def_solv == 6) then
          !! subtract explicit V^N
          call daxpy_(nAsh(iS),-1.0d+00,PCMSCFMO(ip1,2),1,SScont(ip1),1)
          !! add explicit V^N
          call daxpy_(nAsh(iS),+1.0d+00,PCMSCFMO(ip1,2),1,SCFcont(ip1),1)
          !! add implicit FIMO
          call daxpy_(nAsh(iS),+1.0d+00,PCMSSMO(ip1,3),1,SCFcont(ip1),1)
          !! explicit + implicit erfx
          call daxpy_(nAsh(iS),-1.0d+00,PCMSCFMO(ip1,3),1,SCFcont(ip1),1)
        end if
      end do
    end if
  end do

  if (PT2_solv .and. mode==1) then
    call daxpy_(ntBsqr,-1.0d+00,PCMPT2MO(1,3),1,PCMSSMO(1,3),1)
  end if

  if (mode == 2) then
!   call daxpy_(ntBsqr,+1.0d+00,PCMSSMOori(1,3),1,PCMSSMO(1,3),1)
    PCMSSMO(:,3) = PCMSSMO(:,3) + PCMSSMOori(:,3)
  end if

  !!! check the sign

  !call dscal_(SIZE(SCFcont),-1.0d+00,SCFcont,1)
  !call dscal_(SIZE(SScont),-1.0d+00,SScont,1)

! call daxpy_(SIZE(FIMO),1.0d+00,FIMO,1,SScont,1)
  if (mode==1) call daxpy_(SIZE(FIMO),1.0d+00,SScori,1,SScont,1)
  if (mode==2) call daxpy_(SIZE(FIMO),1.0d+00,SScori,1,SScont,1)

! write (6,*) "SS cont"
! call sqprt(sscont,nbas(1))

  !! CISigma rather than CISigma_sa?
  !! effort can be halved (but handling the CI vector is not straightforward)

  !! state-specific quantities (this should not be weighted, but this does not matter)
  !! with def_solv = 3, this is zero (after projection)
  if (mode==1) then
    !! Need to evaluate internal state rotations
    !! However, FIMO and Int2 will result in eigenstates, so we just need to evaluate one-electron state-specific contributions
  ! Call daxpy_(SIZE(FIMO),+1.0D+00,FIMO,1,SScont,1)
  ! Call CISigma_sa(0,state_sym,state_sym,SScont,SIZE(FIMO),Int2,SIZE(Int2),rtmp,1,ipCI,ipCID,.True.)
  ! Call daxpy_(SIZE(FIMO),-1.0D+00,FIMO,1,SScont,1)
    Call CISigma_sa(0,state_sym,state_sym,SScont,SIZE(FIMO),rtmp,1,rtmp,1,ipCI,ipCID,.False.)
  else if (mode==2) then
  ! Call daxpy_(SIZE(FIMO),+1.0D+00,FIMO,1,SScori,1)
  ! Call CISigma_sa(0,state_sym,state_sym,SScori,SIZE(FIMO),Int2,SIZE(Int2),rtmp,1,ipCI,ipCID,.True.)
  ! Call daxpy_(SIZE(FIMO),-1.0D+00,FIMO,1,SScori,1)
    Call CISigma_sa(0,state_sym,state_sym,SScori,SIZE(FIMO),rtmp,1,rtmp,1,ipCI,ipCID,.False.)
  end if
  call dcopy_(nconf1,W(ipCID)%Vec(1+nconf1*(jR-1)),1,vectmp,1)

  !! state averaged quantities (ipCID is overwritten)
  !! it is not actually state-averaged; it is solvent density, so use W_SOLV
  call dswap_(nRoots,weight,1,W_SOLV,1)
  Call CISigma_sa(0,state_sym,state_sym,SCFcont,SIZE(FIMO),rtmp,1,rtmp,1,ipCI,ipCID,.False.)
  call dswap_(nRoots,weight,1,W_SOLV,1)

  !! consider the weight of the derivative
  idsym = state_sym !?
  if (DWSolv%DWZeta /= 0.0d+00) then
    call DWder(2,idsym,SCFcont,SIZE(FIMO),SCFcont,SIZE(FIMO),ISR%RVec, &
               ntash,nRoots,ERASSCF,itoc,LUJOB, &
               nSym,nAsh,nIsh,nBas,ipCM)
  end if

  !! Add the (unweighted) state-specific contributions
  !! Putting real(nRoots) directly into daxpy_ does not work, somehow
  !! Maybe kind=wp is needed?
  if (mode == 1) then
    !! rather, weight should be modified when calling CISigma_sa for SS
    rdum = real(nRoots)
    call daxpy_(nconf1,rdum,vectmp,1,W(ipCID)%Vec(1+nconf1*(kR-1)),1)
!   if (def_solv /= 3) call dcopy_(nconf1*nroots,W(ipCID)%Vec(1),1,vectmpSS,1)
  else if (mode == 2 .and. def_solv /= 3) then
!   call daxpy_(nconf1*nroots,1.0d+00,vectmpSS,1,W(ipCID)%Vec(1),1)
!   write (6,*) "vectmpSS"
!   do i = 1, nRoots
!     write (6,*) "iroot = ", i
!     do j = 1, nconf1
!      write (6,'(i4,f20.10)') j,vectmpSS(j,i)*2.0d+00
!     end do
!   end do
  end if

! write (6,*) "CLag before projection"
! do i = 1, nRoots
!   write (6,*) "iroot = ", i
!   do j = 1, nconf1
!    write (6,'(i4,f20.10)') j,w(ipCID)%Vec(j+nconf1*(i-1))*2.0d+00
!   end do
! end do

  !! State rotations are not needed, because CI vectors are obtained by diagonalization (?)
! call PCM_grad_projection(1,W(ipCI)%Vec(1),W(ipCID)%Vec(1),SLag,vectmpSS,vectmpSA)

  call dscal_(nconf1*nRoots,2.0d+00,W(ipCID)%Vec(1),1)

! if (.not.InvEne .and. InvSCF) then
  if (mode==1 .and. InvSCF) then
!   ISR%Rvec = 0.0d+00
    call ISR_RHS(W(ipCI)%Vec,W(ipCID)%Vec)
    Call ISR_Projection(W(ipCI)%Vec,W(ipCID)%Vec)
    !! check the sign
!   call daxpy_(nRoots**2,+1.0d+00,ISR%Rvec,1,SLag,1)
!   ISR%Rvec = 0.0d+00
!   call dscal_(nRoots**2,2.0d+00,SLag,1)
  end if

! write (6,*) "CLag after projection"
! do i = 1, nRoots
!   write (6,*) "iroot = ", i
!   do j = 1, nconf1
!    write (6,'(i4,f20.10)') j,w(ipCID)%Vec(j+nconf1*(i-1))
!   end do
! end do

! if (mode == 1 .and. def_solv /= 3) call dcopy_(nconf1*nroots,vectmpSS,1,W(ipCID)%Vec(1),1)
! if (mode == 1 .and. def_solv /= 3 .and. InvSCF) call dscal_(nconf1*nroots,0.5d+00,W(ipCID)%Vec(1),1)

  call mma_deallocate(SCFcont)
  call mma_deallocate(SScont)
  call mma_deallocate(SScori)
  call mma_deallocate(vectmp)
  call mma_deallocate(vectmpSS)
  call mma_deallocate(vectmpSA)

  return

  End Subroutine PCM_grad_CLag
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_grad_projection(mode,CI,CIDER,SLAG,vectmpSS,vectmpSA)

  use ipPage, only: W
  use ISRotation, only: InvSCF

  implicit none

  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: CI(ncsf(State_Sym),nRoots),vectmpSS(ncsf(State_Sym),nRoots),vectmpSA(ncsf(State_Sym),nRoots)
  real(kind=wp), intent(inout) :: CIDER(ncsf(State_Sym),nRoots),SLAG(*)

  integer(kind=iwp) :: iRoots,jRoots,i
  real(kind=wp) :: scal
  real(kind=wp), external :: ddot_
  real(kind=wp), allocatable :: DERtmp(:)

  call mma_allocate(DERtmp,ncsf(State_Sym),Label='DERtmp')

  !! Numerical state lagrangian
! if (mode==1) call dcopy_(nRoots**2,[0.0d+00],0,SLAG,1)
! SLAG(1:nRoots**2) = 0.0D+00

! if (mode==1 .and. def_solv/=3) then
  if (mode==1) then

    SLag(1:nRoots**2) = 0.0d+00
!   write (6,'("State Lagrangian")')
    Do iRoots = 1, nRoots
!   write (6,*) "ERASSCF = ", ERASSCF(iRoots),iROots
      Do jRoots = 1, iRoots-1
        scal = DDot_(ncsf(State_Sym),CI(1,iRoots),1,CIDER(1,jRoots),1) &
             - DDot_(ncsf(State_Sym),CI(1,jRoots),1,CIDER(1,iRoots),1)
        if (iStpPCM==1) then
!         write (6,*) "iStpPCM==1"
          if (InvSCF) then
            ! non-iterative internal state rotations, if invariant
            SLag(iRoots+nRoots*(jRoots-1)) = scal/(ERASSCF(jRoots)-ERASSCF(iRoots))
          else
            ! initial residue (for iterative solution)
            scal = DDot_(ncsf(State_Sym),CI(1,iRoots),1,CIDER(1,jRoots),1) &
                 - DDot_(ncsf(State_Sym),CI(1,jRoots),1,CIDER(1,iRoots),1)
            SLag(iRoots+nRoots*(jRoots-1)) = +scal
          end if
        else if (iStpPCM==2.and..not.InvSCF) then
!         write (6,*) "iStpPCM==2"
          SLag(iRoots+nRoots*(jRoots-1)) = +scal!*(ERASSCF(jRoots)-ERASSCF(iRoots))
          if (def_solv==1) &
        scal = DDot_(ncsf(State_Sym),CI(1,iRoots),1,CIDER(1,jRoots),1) &
             + DDot_(ncsf(State_Sym),CI(1,jRoots),1,CIDER(1,iRoots),1)
          SLag(iRoots+nRoots*(jRoots-1)) = scal/2!*(ERASSCF(jRoots)-ERASSCF(iRoots))
        else if (iStpPCM==3) then
          write (6,*) "should not be called..."
          call abend
        end if
!       write (6,'(2i3,2f20.10)') iRoots,jRoots,scal,slag(iroots+nroots*(jroots-1))
      End Do
    End Do

    call dscal_(nRoots**2,-2.0d+00,SLag,1)

!   if (.not.InvSCF .and. iStpPCM==1) call dscal_(nRoots**2,-1.0d+00,ISRot(:,:,1),1)
    if (InvSCF .and. iStpPCM==1) call dscal_(nRoots**2,-1.0d+00,SLag,1)
  end if

  !! Projection
  Do iRoots = 1, nRoots
    Call DCopy_(ncsf(State_Sym),CIDER(1,iRoots),1,DERtmp,1)
    Do jRoots = 1, nRoots
      Scal = DDot_(ncsf(State_Sym),DERtmp,1,CI(1,jRoots),1)
      Call DaXpY_(ncsf(State_Sym),-Scal,CI(1,jRoots),1,CIDER(1,iRoots),1)
    End Do
  End Do

  call mma_deallocate(DERtmp)

  return

  end subroutine PCM_grad_projection
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_mod_ERASSCF(ERASSCF_)

  implicit none

  real(kind=wp), intent(inout) :: ERASSCF_(nRoots)

  real(kind=wp) :: ecorr
  integer(kind=iwp) :: iRoot,iA,jA,ip1,ip2,iS,jS,iRlxRoot_sav,idSym
  logical(kind=iwp) :: isNAC_sav
!
! Subtract (add) the surplus solvation energies from the ERASSCF.
! ERASSCF at present is not eigenvalues of the Hamiltonian, because of some additional
! contributions from the solvent effect, so the purpose of this subroutine is
! to subtract the surplus energies so that ERASSCF is eigenvalues.
!
  idSym = 1

  iRlxRoot_sav = iRlxRoot
  isNAC_sav    = isNAC

  do iRoot = 1, nRoots
  write (6,*) "original erasscf = ", erasscf(iroot),erasscf_(iroot)
    ecorr = 0.0d+00

    do iS = 1, nSym
!     jS=iEOr(is-1,state_sym-1)+1
      !! inactive
      do iA = 1, nIsh(iS)
        ip2 = iA-1 + nOrb(iS)*(iA-1) + ipMat(iS,iS)
        if (def_solv.eq.1) then
          !! - D^SS * V(e,SS) / 2
          ecorr = ecorr - 2.0d+00*pcmscfmo(ip2,3)
        else if (def_solv.eq.2) then
          !! - D^SS * (V(N,SA)+V(e,SA)) / 2
          ecorr = ecorr - 2.0d+00*pcmscfmo(ip2,1) &
          !! + D^SA * V(N,SS) / 2
                        + 2.0d+00*pcmssmo(ip2,2)
        else if (def_solv.eq.3) then
          !! - D^SS * V(e,SA) / 2
          ecorr = ecorr - 2.0d+00*pcmscfmo(ip2,3)
        else if (def_solv.eq.4) then
          !! - D^SS * V(e,SA) / 2
          ecorr = ecorr - 2.0d+00*pcmscfmo(ip2,3)
        else if (def_solv.eq.5) then
          !! - D^SS * V(e,SA) / 2
          ecorr = ecorr - 2.0d+00*pcmscfmo(ip2,3)
        else if (def_solv.eq.6) then
          !! - D^SS * V(e,SA) / 2
          ecorr = ecorr - 2.0d+00*pcmscfmo(ip2,3)
        end if
      end do
    end do
    ecorr = ecorr * 0.5d+00

    if (def_solv==4 .or. def_solv==5 .or. def_solv==6) then
      !! hack for the moment to get state-specific density
      iRlxRoot = iRoot
      isNAC    = .false.
      call PCM_grad_dens(2)
    end if

    !! active
    do iS = 1, nSym
!     jS=iEOr(is-1,state_sym-1)+1
      jS=iEOr(is-1,idSym-1)+1
      do iA = 1, nash(iS)
        do jA = 1, nash(jS)
!         ip1 = iA+nA(iS) + ntAsh*(jA+nA(jS)-1)
          ip2 = nIsh(iS)+iA-1 + nOrb(iS)*(nIsh(jS)+jA-1) + ipMat(iS,jS) ! ipCM(iS)
          if (def_solv.eq.1) then
            !! - D^SS * V(e,SA) / 2
            ecorr = ecorr - 0.5d+00*DSCFMO(iA+nA(iS),jA+nA(jS))*pcmscfmo(ip2,3)
          else if (def_solv.eq.2) then
            !! - D^SS * (V(N,SA)+V(e,SA)) / 2
            ecorr = ecorr - 0.5d+00*DSCFMO(iA+nA(iS),jA+nA(jS))*pcmscfmo(ip2,1) &
            !! + D^SA * V(N,SS) / 2
                            + 0.5d+00*DSCFMO(iA+nA(iS),jA+nA(jS))*pcmssmo(ip2,2)
          else if (def_solv.eq.3) then
            ecorr = ecorr - 0.5d+00*DSCFMO(iA+nA(iS),jA+nA(jS))*pcmscfmo(ip2,3)
          else if (def_solv.eq.4) then
            ecorr = ecorr - 0.5d+00*DSSMO(iA+nA(iS),jA+nA(jS))*pcmscfmo(ip2,3)
          else if (def_solv.eq.5) then
            ecorr = ecorr + DSCFMO(iA+nA(iS),jA+nA(jS))*(pcmscfmo(ip2,2)+0.5d+00*pcmscfmo(ip2,3))
            ecorr = ecorr - DSSMO(iA+nA(iS),jA+nA(jS))*pcmscfmo(ip2,1)
          else if (def_solv.eq.6) then
            ecorr = ecorr + DSCFMO(iA+nA(iS),jA+nA(jS))*(pcmscfmo(ip2,2)+0.5d+00*pcmscfmo(ip2,3))
            ecorr = ecorr - DSSMO(iA+nA(iS),jA+nA(jS))*pcmscfmo(ip2,2)
          end if
        end do
      end do
    end do

    ERASSCF_(iRoot) = ERASSCF_(iRoot) - ecorr
  write (6,*) "modified erasscf = ", erasscf_(iroot)
  end do

  if (def_solv==4 .or. def_solv==5 .or. def_solv==6) then
    iRlxRoot = iRlxRoot_sav
    isNAC    = isNAC_sav
    call PCM_grad_dens(2)
  end if

  End Subroutine PCM_mod_ERASSCF
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_grad_PT2()

  implicit none

  logical(kind=iwp) :: Found
  integer(kind=iwp) :: iAO1,iAO2,iSym,mData
  real(kind=wp), external :: ddot_
!
! Add implicit contributions due to the unrelaxed CASPT2 density
!
  !! check D1aoVar
  !! if it is non-zero, RFpert is on in CASPT2,
  !! so the response contributions should be evaluated.
  !! Otherwise, exit this subroutine
  call qpg_dArray('D1aoVar',Found,mData)
  if (.not.Found .or. mData==0) then
    call warningMessage(2,'Gradients without RFpert in &CASPT2 is incorrect.')
    return
  end if

  !! PCMPT2MO here is used as a temporary array
  call Get_dArray('D1aoVar',PCMPT2MO,ntBtri)
  if (ddot_(ntBtri,PCMPT2MO,1,PCMPT2MO,1).le.1.0d-10) return

  PT2_solv = .true.

  call unfold(PCMPT2MO,ntBtri,DZAO,ntBsqr,nSym,nBas)
  call PCM_grad_D2V(DZAO,PCMPT2MO,PCMZAO,.false.,.false.,.false.,.true.)

  return

  End Subroutine PCM_grad_PT2

End Module PCM_grad
