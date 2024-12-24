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
Module PCM_alaska

  use definitions, only: iwp,wp
  use pso_stuff, only: nDens
  use rctfld_module, only: nTS, PCM
  use stdalloc, only: mma_allocate, mma_deallocate

  Implicit None

  !! second ASCs used in many ways
  !!! PCM_SQ_ind(1,:) :: nuclear charges
  !!! PCM_SQ_ind(2,:) :: electron charges
  real(kind=wp), Allocatable :: PCM_SQ_ind(:,:)
  !! unrotated state-specific ASC
  real(kind=wp), Allocatable :: PCM_SQ_ori(:,:)
  !! CMO obtained in RASSCF (PCM + CMO)
  real(kind=wp), Allocatable :: PCMO(:)
  !! state-averaged density matrix (in AO)
  real(kind=wp), Allocatable :: DSA_AO(:)
  !! unrotated state-specific density matrix (in AO)
  real(kind=wp), Allocatable :: DSS_ori(:)

  logical(kind=iwp) :: lSA = .false.
  integer(kind=iwp) :: def_solv = 0_iwp

  contains
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_alaska_lSA()

  implicit none

  integer(kind=iwp) :: iGo
  Character(Len=8) Method

  Call Get_cArray('Relax Method',Method,8)
  lSA=.false.
  if (Method=='CASSCFSA' .or. Method=='DMRGSCFS' .or. Method=='GASSCFSA' .or. &
      Method=='RASSCFSA' .or. Method=='CASPT2  ') then
    Call Get_iScalar('SA ready',iGo)
    if (iGO==1) lSA=.true.
  end if

  if (.not.PCM) return
  if (.not.lSA) return

  call PCM_alaska_init()

  return

  End Subroutine PCM_alaska_lSA
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_alaska_init()

  use Sizes_of_Seward, only: S

  implicit none

  integer(kind=iwp) :: mCMO,iPCMRoot

  if (def_solv==0_iwp) then
!   call Get_iScalar('PCMRoot',iPCMRoot)
    call Get_iScalar('RF CASSCF root',iPCMRoot)
    if (iPCMRoot>  0) def_solv = 1_iwp
    if (iPCMRoot== 0) def_solv = 3_iwp
    if (iPCMRoot==-1) def_solv = 4_iwp
    if (iPCMRoot==-2) def_solv = 5_iwp
    if (iPCMRoot==-3) def_solv = 6_iwp
  end if

  call mma_allocate(PCM_SQ_ind,2,nTS,Label='PCM_SQ_ind')
  PCM_SQ_ind = 0.0d+00

! if (def_solv /= 3) then
    call mma_allocate(PCM_SQ_ori,2,nTS,Label='PCM_SQ_ori')
    PCM_SQ_ori = 0.0d+00
! end if

  mCMO = S%n2Tot
  Call mma_allocate(PCMO,mCMO,Label='CMO')
  Call Get_dArray_chk('Last orbitals',PCMO,mCMO)

  return

  End Subroutine PCM_alaska_init
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_alaska_final()

  implicit none

  if (.not.PCM) return
  if (.not.lSA) return

  if (allocated(PCM_SQ_ind)) call mma_deallocate(PCM_SQ_ind)
  if (allocated(PCM_SQ_ori)) call mma_deallocate(PCM_SQ_ori)
  if (allocated(PCMO))       call mma_deallocate(PCMO)
  if (allocated(DSA_AO))     call mma_deallocate(DSA_AO)
  if (allocated(DSS_ori))    call mma_deallocate(DSS_ori)

  return

  End Subroutine PCM_alaska_final
!
!-----------------------------------------------------------------------
!
  Subroutine PCM_alaska_prep()

  use Basis_Info, only: nBas
  use etwas, only: nAsh,nIsh
! use NAC, only: isNAC
  use PCM_arrays, only: PCM_SQ
  use pso_stuff, only: nDens
  use Symmetry_Info, only: nIrrep

  implicit none

  real(kind=wp), allocatable :: Dtmp(:),Tmp(:),D1AV(:),h1(:),TwoHam(:)

  integer(kind=iwp) :: iIrrep,nAct,nG1,ij,iBas,jBas,iPCMRoot
  logical(kind=iwp) :: First,Dff,NonEq
  real(kind=wp) :: RepNuc

  call Set_Basis_Mode('Valence')
  call Setup_iSD()
  call Get_iArray('nAsh',nAsh,nIrrep)
  Call Get_iArray('nIsh',nIsh,nIrrep)

! First, construct state-averaged density matrix (in AO)

  nAct = 0
  nDens = 0
  Do iIrrep = 0, nIrrep-1
    nAct = nAct + nAsh(iIrrep)
    nDens = nDens + nBas(iIrrep)*(nBas(iIrrep)+1)/2
  End Do
  nG1 = nAct*(nAct+1)/2

  call mma_allocate(DSA_AO,nDens,Label='DSA_AO')
  call mma_allocate(Dtmp,nDens,Label='Dtmp') ! ntBtri length

  if (def_solv==1) then
    !! state-specific density, which polarizes ASC during SCF
    !! So, it is not necessarily identical to D1ao, determined by RlxRoot
    !! This density is already dynamically or fixed weighted
    DSA_AO = 0.0D+00
    Call Get_dArray_chk('D1ao_PCM',DSA_AO,nDens)
!   Call Get_dArray_chk('D1ao',DSA_AO,nDens)
    !! "D1ao_PCM" is from &RASSCF
    !! "Reaction field" is from DrvRF
  else
    !! inactive
    call mma_allocate(Tmp,nDens*2,Label='Tmp') ! ntBtri*2
    Call Get_D1I(PCMO,Dtmp,Tmp,nIsh,nBas,nIrrep)
    if (allocated(Tmp)) call mma_deallocate(Tmp)
    call dcopy_(nDens,Dtmp,1,DSA_AO,1)

    !! active
    Call mma_allocate(D1AV,nG1,Label='D1AV')
    Call Get_dArray_chk('D1av',D1AV,nG1)
    Call Get_D1A(PCMO,D1AV,Dtmp,nIrrep,nbas,nish,nash,ndens)
    call daxpy_(nDens,1.0d+00,Dtmp,1,DSA_AO,1)
    ij = -1
    Do iIrrep = 0, nIrrep-1
      Do iBas = 1, nBas(iIrrep)
        Do jBas = 1, iBas-1
          ij = ij + 1
          DSA_AO(1+ij) = 2.0d+00*DSA_AO(1+ij)
        end do
        ij = ij + 1
      end do
    end do

    if (allocated(D1AV)) call mma_deallocate(D1AV)
  end if
!         write (6,*) "DSA_AO"
!         call sqprt(DSA_AO,nbas(0))

! Second, obtain ASCs induced by the effective density (after Z-vector)
! The induced charges are put in PCM_SQ_ind

  ! some setup are required
  call Get_D1ao_Var(Dtmp,nDens) ! Read from D1aoVar (if not found, D1ao then)
!         write (6,*) "D1ao_Var"
!         call sqprt(Dtmp,nbas(0))
  call mma_allocate(h1,nDens,Label='h1')
  call mma_allocate(TwoHam,nDens,Label='TwoHam')
  RepNuc = 0.0d+00
  First = .true.
  Dff = .false.
  NonEq = .false.
  Call DrvPCM(h1,TwoHam,Dtmp,RepNuc,nDens,First,Dff,NonEq)
  Call Get_dArray('PCM Charges',PCM_SQ_ind,2*nTs)
  !! in order to remove nuclear contributions
! if (isNAC) PCM_SQ_ind(1,:) = 0.0d+00

! Third, obtain ASCs induced by the state-averaged density
! The standard charges are put in PCM_SQ.
! They have been read in Init_PCM, but they are wrong because they are overwritten during MCLR

  Call DrvPCM(h1,TwoHam,DSA_AO,RepNuc,nDens,First,Dff,NonEq)
  Call Get_dArray('PCM Charges',PCM_SQ,2*nTs)

  !! PCM_SQ_ind could have been rotated in MCLR
  !! we need unrotated state-specific density and ASC (w/o response) for some definitions
! if (def_solv /= 3) then
    call mma_allocate(DSS_ori,nDens,Label='DSS_ori')
    DSS_ori = 0.0d+00

    Call Get_dArray_chk('D1ao',DSS_ori,nDens)
!         write (6,*) "DSS_ori"
!         call sqprt(DSS_ori,nbas(0))
    Call DrvPCM(h1,TwoHam,DSS_ori,RepNuc,nDens,First,Dff,NonEq)
    Call Get_dArray('PCM Charges',PCM_SQ_ori,2*nTs)
! end if

! do ij = 1, nts
! write (6,'(i4,6f20.10)') ij,pcm_sq(1,ij),pcm_sq(2,ij),pcm_sq_ind(1,ij),pcm_sq_ind(2,ij),pcm_sq_ori(1,ij),pcm_sq_ori(2,ij)
! end do

! if (.true.) then
!   Call Get_dArray_chk('D1ao_PCM',DSA_AO,nDens)
!   Call DrvPCM(h1,TwoHam,DSA_AO,RepNuc,nDens,First,Dff,NonEq)
! end if

  call Free_iSD()

  !! tentative patch for MECI calculations with XMS-type CASPT2
  !! "Reaction Field" has been overwritten, but the correct reaction field is needed
  !! for computing H0 in XMS-type calculations
  !! For PTED, ALASKA should not be called.
  h1 = 0.0d+00
  TwoHam = 0.0d+00
  Call DrvRF(h1,TwoHam,DSA_AO,RepNuc,nDens,First,Dff,NonEq)

  if (allocated(Dtmp))   call mma_deallocate(Dtmp)
  if (allocated(h1))     call mma_deallocate(h1)
  if (allocated(TwoHam)) call mma_deallocate(TwoHam)

  End Subroutine PCM_alaska_prep

End Module PCM_alaska
