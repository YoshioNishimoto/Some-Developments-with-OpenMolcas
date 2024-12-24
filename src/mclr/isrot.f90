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
!
! This module is used when we need to consider rotations between
! internal states, which are usually state-averaged RASSCF states. The
! standard SA-RASSCF (CASSCF) is invariant with respect to rotations
! within internal states, but not for the following example:
!   1. unequally weighted SA-RASSCF
!   2. SA-RASSCF/PCM (unless state-averaging and density that polarizes
!      ASCs are the same)
!   3. CASPT2/RASPT2
! In some cases, the Lagrangian multipliers relevant to the rotations
! can be obtained without iteration. For instance, for 3, if the
! reference SA-CASSCF is invariant, the multipliers can be immediately
! obtained. Or, for 2, if the density that polarizes ASCs are equally
! state-averaged, the same applies. For other cases, they shoould be
! obtained iteratively.
!
! (InvSCF, InvEne) = (.true., .true.)
!   -> conventional SA-CASSCF/RASSCF and SA-CASSCF/PCM (with def_solv = 3)
!      There is no need to consider the internal state rotations
! (InvSCF, InvEne) = (.true., .false.)
!   -> CASPT2 (w/o PCM), SA-CASSCF/PCM (with def_solv = 4)
!      Internal state rotations are non-iteratively evaluated by rotating
!      the 1e and 2e density matrix in rhs_sa and rhs_nac through the
!      SLag matrix.
! (InvSCF, InvEne) = (.false., .false.)
!   -> anything based on SA-CASSCF/PCM (with def_solv /= 3,4)
!      Some internal rotations may be evaluated non-iteratively as in
!      (InvSCF, InvEne) = (.true., .false.) . This comes from the
!      rotation in MS-CASPT2 and should not be considered as rotation
!      parameters to be optimized, since it is not iteratively
!      determined during energy.
!      The other internal rotations come from the overlap between the
!      original CI and CI derivative. This is optimized during Z-vector.
! So, if InvEne = .true., no special treatment is required.
!
Module ISRotation

  use definitions, only: iwp,wp,u6
  use stdalloc, only: mma_allocate, mma_deallocate

  Implicit None

  ! Whether the target energy is invariant or non-invariant wrt internal
  ! state rotations. The conventional SA-CASSCF/RASSCF is invariant, so
  ! InvEne = .true. For the above cases, the energy is non-invariant,
  ! so InvEne = .false. (for instance, CASPT2 or SA-CASSCF with PCM etc.
  logical(kind=iwp) :: InvEne = .true.
  ! Whether the reference SCF is invariant or not. That means rotation
  ! between internal states are considered during Z-vector.
  ! See the above note for the details
  ! In short, this flag determines whether internal state rotations are
  ! iterative (InvSCF = .false.) or non-iterative (InvSCF = .true.) determined
  logical(kind=iwp) :: InvSCF = .true.
  ! Whether we consider the coupling of O-S and C-S coupling
  ! Purely development purpose
  logical(kind=iwp) :: IntRotOff = .true.
  ! Activate conjugate gradient squared method
  ! This option may be useful for non-PCM runs
  logical(kind=iwp) :: CGS = .false.

  type ISR_param
    real(kind=wp), Allocatable :: Ap(:,:)   ! results of A*p
    real(kind=wp), Allocatable :: p(:,:)    ! trial vector during CG
    real(kind=wp), Allocatable :: Rvec(:,:) ! initial residue (RHS)
    real(kind=wp), Allocatable :: prec(:,:) ! preconditioned something
    real(kind=wp), Allocatable :: Xvec(:,:) ! solution
    real(kind=wp), Allocatable :: Pvec(:,:) ! P vector (for CGS)
    real(kind=wp), Allocatable :: Qvec(:,:) ! Q vector (for CGS)
    real(kind=wp), Allocatable :: Uvec(:,:) ! U vector (for CGS)
    real(kind=wp), Allocatable :: R0(:,:)   ! initial residue (for CGS)
  end type ISR_param

  type(ISR_param) :: ISR

  integer(kind=iwp), pointer :: nRoots
  integer(kind=iwp), pointer :: ncsf(:)
  integer(kind=iwp), pointer :: State_Sym

  real(kind=wp),     pointer :: ERASSCF(:)

contains
!
!-----------------------------------------------------------------------
!
  Subroutine ISR_init(iPL,nRoots_,ncsf_,State_Sym_, &
                      ERASSCF_)

  implicit none

  integer(kind=iwp), intent(in) :: iPL
  integer(kind=iwp), intent(in), target :: nRoots_,ncsf_(*),State_Sym_

  real(kind=wp), intent(in), target :: ERASSCF_(*)

  if (iPL>=2) then
    if (InvSCF) then
      write (6,*) "internal rotations are non-iteratively determined"
    else if (.not.InvSCF) then
      write (6,*) "internal rotations are iteratively determined"
    end if
  end if

  nRoots            => nRoots_
  ncsf(1:8)         => ncsf_(1:8)
  State_Sym         => State_Sym_

  ERASSCF(1:nRoots) => ERASSCF_(1:nRoots)

! if (.not.CGS) return
  if (.not.CGS .and. InvEne) return

  call mma_allocate(ISR%Ap  ,nRoots,nRoots,Label='ISR%Ap')
  call mma_allocate(ISR%p   ,nRoots,nRoots,Label='ISR%p')
  call mma_allocate(ISR%Rvec,nRoots,nRoots,Label='ISR%Rvec')
  call mma_allocate(ISR%prec,nRoots,nRoots,Label='ISR%prec')
  call mma_allocate(ISR%Xvec,nRoots,nRoots,Label='ISR%Xvec')
  call mma_allocate(ISR%Pvec,nRoots,nRoots,Label='ISR%Pvec')
  call mma_allocate(ISR%Qvec,nRoots,nRoots,Label='ISR%Qvec')
  call mma_allocate(ISR%Uvec,nRoots,nRoots,Label='ISR%Uvec')
  call mma_allocate(ISR%R0  ,nRoots,nRoots,Label='ISR%R0')

  ISR%Ap   = 0.0d+00
  ISR%p    = 0.0d+00
  ISR%Rvec = 0.0d+00
  ISR%prec = 0.0d+00
  ISR%Xvec = 0.0d+00
  ISR%Pvec = 0.0d+00
  ISR%Qvec = 0.0d+00
  ISR%Uvec = 0.0d+00
  ISR%R0   = 0.0d+00

  !! actually, if InvEne = .false., it is usually better to use CGS
  CGS = .true.

  if (CGS .and. iPL>=2) then
    Write(6,'(x,A)') 'Using the preconditioned conjugate '// &
                     'gradient squared (CGS) method'
  end if

  End Subroutine ISR_init
!
!-----------------------------------------------------------------------
!
  Subroutine ISR_final()

  implicit none

! if (.not.CGS) return
  if (.not.CGS .and. InvEne) return

  if (allocated(ISR%Ap))    call mma_deallocate(ISR%Ap)
! if (allocated(ISR%p))     call mma_deallocate(ISR%p)
  if (allocated(ISR%Rvec))  call mma_deallocate(ISR%Rvec)
  if (allocated(ISR%prec))  call mma_deallocate(ISR%prec)
  if (allocated(ISR%Xvec))  call mma_deallocate(ISR%Xvec)
  if (allocated(ISR%Pvec))  call mma_deallocate(ISR%Pvec)
  if (allocated(ISR%Qvec))  call mma_deallocate(ISR%Qvec)
  if (allocated(ISR%Uvec))  call mma_deallocate(ISR%Uvec)
  if (allocated(ISR%R0))    call mma_deallocate(ISR%R0)

  End Subroutine ISR_final
!
!-----------------------------------------------------------------------
!
  Subroutine ISR_final2()

  implicit none

! if (.not.CGS) return
  if (.not.CGS .and. InvEne) return

  if (allocated(ISR%p))     call mma_deallocate(ISR%p)

  End Subroutine ISR_final2
!
!-----------------------------------------------------------------------
!
  Subroutine ISR_RHS(CI,CIDER)
  implicit none

  real(kind=wp), intent(in) :: CI(ncsf(State_Sym),nRoots),CIDER(ncsf(State_Sym),nRoots)

  integer(kind=iwp) :: i,j
  real(kind=wp) :: scal
  real(kind=wp), external :: ddot_
!
! Construct the RHS (with minus) of state rotation
! Note that CIDER (W(ipST)%Vec) already has the minus sign
!
  do i = 1, nRoots
    do j = 1, i-1
      scal = DDot_(ncsf(State_Sym),CI(1,i),1,CIDER(1,j),1) &
           - DDot_(ncsf(State_Sym),CI(1,j),1,CIDER(1,i),1)
      if (InvSCF) then
        ! non-iterative internal state rotations, if invariant
!       ISR%Rvec(i,j) = ISR%Rvec(i,j) + 2.0d+00*scal/(ERASSCF(j)-ERASSCF(i))
        ISR%Rvec(i,j) = ISR%Rvec(i,j) + 1.0d+00*scal/(ERASSCF(j)-ERASSCF(i))
      else if (.not.InvSCF) then
        ! initial residue (for iterative solution)
!       ISR%Rvec(i,j) = ISR%Rvec(i,j) + 2.0d+00*scal
        ISR%Rvec(i,j) = ISR%Rvec(i,j) + 1.0d+00*scal
      end if
    end do
  end do

! ISR%Rvec =  ISR%Rvec*0.5d+00 !! because CIDER has been multiplied by two

  write (6,*) "initial state rotation"
  call sqprt(isr%rvec,nroots)

  Return

  End Subroutine ISR_RHS
!
!-----------------------------------------------------------------------
!
  Subroutine ISR_projection(CI,CIDER)

  implicit none

  real(kind=wp), intent(in) :: CI(ncsf(State_Sym),nRoots)
  real(kind=wp), intent(inout) :: CIDER(ncsf(State_Sym),nRoots)

  integer(kind=iwp) :: i,j
  real(kind=wp) :: scal
  real(kind=wp), external :: ddot_
!
! Project out the internal rotation contribution
!
  do i = 1, nRoots
    do j = 1, nRoots
      scal = DDot_(ncsf(State_Sym),CI(1,i),1,CIDER(1,j),1)
      call daxpy_(ncsf(State_Sym),-scal,CI(1,i),1,CIDER(1,j),1)
    end do
  end do

  Return

  End Subroutine ISR_projection
!
!-----------------------------------------------------------------------
!
  Subroutine ISR_TimesE2(CI,CIDER,W_SOLV)

  implicit none

  real(kind=wp), intent(in) :: CI(ncsf(State_Sym),nRoots),W_SOLV(nRoots)
  real(kind=wp), intent(inout) :: CIDER(ncsf(State_Sym),nRoots)

  integer(kind=iwp) :: i,j
  real(kind=wp) :: scal
  real(kind=wp), external :: ddot_
!
! Compute the S-S block of the Ap operation
!
  if (.not.InvSCF) then
    do i = 1, nRoots
      do j = 1, i-1
        scal = 0.0d+00
        if (IntRotOff) then
          scal = DDot_(ncsf(State_Sym),CI(1,j),1,CIDER(1,i),1)*W_SOLV(j) &
               - DDot_(ncsf(State_Sym),CI(1,i),1,CIDER(1,j),1)*W_SOLV(i)
          scal = scal*dble(nRoots) !! because CIDER has been multiplied by 1/nRoots
        end if
        ISR%Ap(i,j) = ISR%Ap(i,j) + scal - (ERASSCF(j)-ERASSCF(i))*ISR%p(i,j)*2.0d+00
      end do
      ISR%Ap(i,i) = ISR%Ap(i,i) + ISR%p(i,i)
    end do

    !! CIDER has been multiplied by 1/nRoots only, so compensate this
    do i = 1, nRoots
      call dscal_(ncsf(State_sym),W_SOLV(i)*dble(nRoots),CIDER(1,i),1)
    end do
  end if

  Return

  End Subroutine ISR_TimesE2
!
!-----------------------------------------------------------------------
!
  Subroutine DMInvISR(ISRotIn,ISRotOut)

  implicit none

  real(kind=wp), intent(in)    :: ISRotIn(nRoots,nRoots)
  real(kind=wp), intent(inout) :: ISRotOut(nRoots,nRoots)

  integer(kind=iwp) :: i,j
!
! diagonal preconditioning for the internal state rotation Hessian
!
  do i = 1, nRoots
    do j = 1, i-1
      ISRotOut(i,j) = -0.5d+00*ISRotIn(i,j)/(ERASSCF(j)-ERASSCF(i))
    end do
    ISRotOut(i,i) = ISRotIn(i,i)
  end do

  Return

  End Subroutine DMInvISR

End Module ISRotation
