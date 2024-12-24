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
! (Preconditioned) Conjugate Gradient Squared (CGS) method
! Ref: Sonneveld, P. SIAM J. Sci. Stat. Comput. 1989, 10, 36-52.
!      Itoh, S.; Sugihara, M. Trans. Jpn. Soc. Ind. Appl. Math. 2013, 23, 253-286.
! The preconditioning comes from Algorithm 4 in the second ref (Japanese literature).
!
! CGS does not assume the Hessian is symmetric positive definite.
! If PCM or dynamically weighted things, CGS is preferred.
! The actual calculation is cgs_x.F90.
!
Module cgs_mod

  use definitions, only: iwp,wp,u6
  use stdalloc, only: mma_allocate, mma_deallocate

  Implicit None

  ! Activate conjugate gradient squared method
  ! This option may be useful for non-PCM runs
  logical(kind=iwp) :: CGS = .false.

  type CGS_type
    !! CI rotations
    integer(kind=iwp) :: ipPvec
    integer(kind=iwp) :: ipQvec
    integer(kind=iwp) :: ipUvec
    integer(kind=iwp) :: ipR0
    !! orbital rotations
    real(kind=wp), Allocatable :: Pvec(:)
    real(kind=wp), Allocatable :: Qvec(:)
    real(kind=wp), Allocatable :: Uvec(:)
    real(kind=wp), Allocatable :: R0(:)
  end type CGS_type

  type(CGS_type) :: CGSvec

contains
!
!-----------------------------------------------------------------------
!
  Subroutine CGS_init(nconf1,nRoots,nDens2)

  implicit none

  integer(kind=iwp), intent(in) :: nconf1,nRoots,nDens2
  integer(kind=iwp), external :: ipGet

  Call mma_allocate(CGSvec%Pvec,nDens2+6,Label='Pvec')
  Call mma_allocate(CGSvec%Qvec,nDens2+6,Label='Qvec')
  Call mma_allocate(CGSvec%Uvec,nDens2+6,Label='Uvec')
  Call mma_allocate(CGSvec%R0,nDens2+6,Label='R0')

  CGSvec%ipPvec=ipGet(nconf1*nroots)
  CGSvec%ipQvec=ipGet(nconf1*nroots)
  CGSvec%ipUvec=ipGet(nconf1*nroots)
  CGSvec%ipR0  =ipGet(nconf1*nroots)

  End Subroutine CGS_init
!
!-----------------------------------------------------------------------
!
  Subroutine CGS_final()

  implicit none

  Call mma_deallocate(CGSvec%Pvec)
  Call mma_deallocate(CGSvec%Qvec)
  Call mma_deallocate(CGSvec%Uvec)
  Call mma_deallocate(CGSvec%R0)

  End Subroutine CGS_final

End Module cgs_mod
