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
! (InvSCF, InvEne) = (.false., .true.)
!   -> unequally (including dynamically) weighted SA-MCSCF
!      The initial residue of the internal rotation is zero, though.
!
Module ISRotation

  use definitions, only: iwp,wp,u6
  use stdalloc, only: mma_allocate, mma_deallocate

  Implicit None

  ! Whether the target energy, specified by RLXROOT, is invariant or
  ! non-invariant wrt internal state rotations. The conventional
  ! SA-MCSCF is invariant, so InvEne = .true. For the above cases, the
  ! energy is non-invariant, so InvEne = .false. (for instance, CASPT2
  ! or SA-MCSCF with PCM etc.). This flag only controls the evaluation
  ! of the initial internal rotation.
  logical(kind=iwp) :: InvEne = .true.
  ! Whether the reference SCF (SA-MCSCF energy) is invariant or not.
  ! That means rotation between internal states are considered during
  ! Z-vector.
  ! See the above note for details
  ! In short, this flag determines whether internal state rotations are
  ! iteratively (InvSCF = .false.) or non-iteratively (InvSCF = .true.)
  ! determined
  logical(kind=iwp) :: InvSCF = .true.
  ! Whether we consider the coupling of O-S and C-S coupling
  ! Purely development purpose
  logical(kind=iwp) :: IntRotOff = .true.
  ! unequal state-averag?
  logical(kind=iwp) :: unequal_SA = .false.

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
  real(kind=wp),     pointer :: weight(:)   ! weight of the state-averaging

  logical(kind=iwp), pointer :: do_RF

contains
!
!-----------------------------------------------------------------------
!
  Subroutine ISR_init(iPL,do_RF_,PT2,nRoots_,ncsf_,State_Sym_,def_solv, &
                      ERASSCF_,weight_)
  use DWSol, only: DWSCF
  use cgs_mod, only: CGS

  implicit none

  integer(kind=iwp), intent(in) :: iPL,def_solv
  logical(kind=iwp), intent(in) :: PT2
  logical(kind=iwp), intent(in), target :: do_RF_
  integer(kind=iwp), intent(in), target :: nRoots_,ncsf_(*),State_Sym_

  real(kind=wp), intent(in), target :: ERASSCF_(*),weight_(*)

  integer(kind=iwp) :: iR

  !! Check whether we need to consider internal state rotations
  !! iteratively or non-interatively
  InvEne = .true.
  if (.not.InvSCF) InvEne = .false.
  if (PT2 .and. nRoots_ > 0) InvEne = .false.
  if (do_RF_ .and. def_solv /= 3) InvEne = .false.
  do iR = 2, nRoots_
    if (weight_(1).ne.weight_(iR)) then
      InvSCF = .false.
      InvEne = .false.
    end if
  end do

  if (iPL>=2) then
    if (InvSCF) then
      write (6,*) "internal rotations are non-iteratively determined"
    else if (.not.InvSCF) then
      write (6,*) "internal rotations are iteratively determined"
    end if
  end if

  do_RF             => do_RF_

  nRoots            => nRoots_
  ncsf(1:8)         => ncsf_(1:8)
  State_Sym         => State_Sym_

  ERASSCF(1:nRoots) => ERASSCF_(1:nRoots)
  weight(1:nRoots)  => weight_(1:nRoots)

  unequal_SA = .false.
  do iR = 2, nRoots
    if (weight(1).ne.weight(iR)) unequal_SA = .true.
  end do

  !! If the target energy is invariant, nothing is needed
! if (.not.CGS) return
  if (.not.CGS .and. InvEne) return

  !! Otherwise, we need state rotations, so allocate the minimum vectors
  !! for state-state rotations
  call mma_allocate(ISR%Ap  ,nRoots,nRoots,Label='ISR%Ap')
  call mma_allocate(ISR%p   ,nRoots,nRoots,Label='ISR%p')
  call mma_allocate(ISR%Rvec,nRoots,nRoots,Label='ISR%Rvec')
  call mma_allocate(ISR%prec,nRoots,nRoots,Label='ISR%prec')
  call mma_allocate(ISR%Xvec,nRoots,nRoots,Label='ISR%Xvec')

  ISR%Ap   = 0.0d+00
  ISR%p    = 0.0d+00
  ISR%Rvec = 0.0d+00
  ISR%prec = 0.0d+00
  ISR%Xvec = 0.0d+00

  if (.not.CGS .and. (.not.do_RF .and. .not.DWSCF%do_DW)) return

  !! actually, if InvEne = .false., it is usually better to use CGS
  CGS = .true.

  if (CGS .and. iPL>=2) then
    Write(6,'(x,A)') 'Using the preconditioned conjugate '// &
                     'gradient squared (CGS) method'
  end if

  !! The rest is used for CGS
  call mma_allocate(ISR%Pvec,nRoots,nRoots,Label='ISR%Pvec')
  call mma_allocate(ISR%Qvec,nRoots,nRoots,Label='ISR%Qvec')
  call mma_allocate(ISR%Uvec,nRoots,nRoots,Label='ISR%Uvec')
  call mma_allocate(ISR%R0  ,nRoots,nRoots,Label='ISR%R0')

  ISR%Pvec = 0.0d+00
  ISR%Qvec = 0.0d+00
  ISR%Uvec = 0.0d+00
  ISR%R0   = 0.0d+00

  End Subroutine ISR_init
!
!-----------------------------------------------------------------------
!
  Subroutine ISR_final()

  use cgs_mod, only: CGS

  implicit none

! if (.not.CGS) return
  if (.not.CGS .and. InvEne) return

  !! Not quite safe for unequally weighted SA-RASSCF...
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

  use cgs_mod, only: CGS

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
!
  do i = 1, nRoots
    do j = 1, i-1
      scal = DDot_(ncsf(State_Sym),CI(1,i),1,CIDER(1,j),1) &
           - DDot_(ncsf(State_Sym),CI(1,j),1,CIDER(1,i),1)
      scal = -scal ! CIDER (W(ipST)%Vec) has been multiplied by -one in wfctl_sa
      if (InvSCF) then
        ! non-iterative internal state rotations, if invariant
!       ISR%Rvec(i,j) = ISR%Rvec(i,j) + 2.0d+00*scal/(ERASSCF(j)-ERASSCF(i))
        ISR%Rvec(i,j) = ISR%Rvec(i,j) + 1.0d+00*scal/(ERASSCF(i)-ERASSCF(j))
      else if (.not.InvSCF) then
        ! initial residue (for iterative solution)
!       ISR%Rvec(i,j) = ISR%Rvec(i,j) + 2.0d+00*scal
        ISR%Rvec(i,j) = ISR%Rvec(i,j) + 1.0d+00*scal
      end if
    end do
    write (6,*) "diagonal = ", DDot_(ncsf(State_Sym),CI(1,i),1,CIDER(1,i),1)
  end do

  ISR%Rvec = -ISR%Rvec !! the initial residue is the negative of the partial derivative

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
  Subroutine ISR_TimesE2(MODE,CI,CIDER)

  use DWSol, only: DWSCF,DWSolv,DWSol_Der

  implicit none

  integer(kind=iwp), intent(in) :: MODE
  real(kind=wp), intent(in) :: CI(ncsf(State_Sym),nRoots)
  real(kind=wp), intent(inout) :: CIDER(ncsf(State_Sym),nRoots)

  integer(kind=iwp) :: i,j
  real(kind=wp) :: scal,fact
  real(kind=wp), external :: ddot_
  real(kind=wp), allocatable :: DERHII(:),DEROMG(:)
!
! Compute the C-S and S-S blocks of the Ap operation
!
  !! if both are true, this subroutine is called twice, so the S-S block has to be halved
  fact = 1.0d+00
  if (do_RF .and. unequal_SA) fact = 0.5d+00

  if (.not.InvSCF) then
    do i = 1, nRoots
      do j = 1, i-1
        scal = 0.0d+00
        if (IntRotOff) then
          scal = DDot_(ncsf(State_Sym),CI(1,i),1,CIDER(1,j),1) &
               - DDot_(ncsf(State_Sym),CI(1,j),1,CIDER(1,i),1)
        end if
        ISR%Ap(i,j) = ISR%Ap(i,j) + scal + (ERASSCF(i)-ERASSCF(j))*ISR%p(i,j)*2.0d+00*fact
      end do
      !! The diagonal contribution is for dynamically weighted methods
      ISR%Ap(i,i) = ISR%Ap(i,i) + ISR%p(i,i)*fact
    end do
  end if
! call sqprt(isr%ap,nroots)

  !! Derivative of H for DW-MCSCF is evaluated with CI derivatives
  !! DW solvation is evaluated in DWder
  if (mode==1.and.DWSCF%do_DW) then
    call mma_allocate(DEROMG,nRoots,Label='DEROMG')
    DEROMG = 0.0D+00
    if (IntRotOff) then
      do i = 1, nRoots
        ! CIDER has been scaled with the weight in cisigma_sa
        if (weight(i).ge.1.0d-08) then
          DEROMG(i) = DDot_(ncsf(State_Sym),CI(1,i),1,CIDER(1,i),1)/weight(i)
        else
          DEROMG(i) = 0.0d+00 !! under consideration
        end if
      end do
    end if

    call mma_allocate(DERHII,nRoots,Label='DERHII')
    DERHII = 0.0D+00
    call DWSol_Der(mode,DEROMG,DERHII,ERASSCF,weight)

    do i = 1, nRoots
      !! not sure why 1/4
      !! One reason is CIDER has been scaled by two, and the other is ?
      ISR%Ap(i,i) = ISR%Ap(i,i) + DERHII(i)*0.25d+00
    end do
    call mma_deallocate(DERHII)
    call mma_deallocate(DEROMG)
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
!       for degeneracy?
        ISRotOut(i,j) = +0.5d+00*ISRotIn(i,j)/(ERASSCF(i)-ERASSCF(j))
    end do
    ISRotOut(i,i) = ISRotIn(i,i)
  end do

  Return

  End Subroutine DMInvISR

End Module ISRotation
