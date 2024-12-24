!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

! Variables related to CASPT2 gradients
! TODO: move here everything that is in caspt2_grad.fh
module caspt2_gradient

  use definitions, only: iwp,wp

  ! some gradient stuff
  integer(kind=iwp) :: iVecL           = 7_iwp ! Solution of the Lambda equation
  ! iVecG (G is probably gradient stuff) is used in
  ! caspt2_res.f to temporarily store residual vectors in solving the lambda equation
  ! sigder.f and clagx.f to temporarily store derivatives of overlap
  integer(kind=iwp) :: iVecG           = 8_iwp
  integer(kind=iwp) :: idSDMat(8,13)   = 0_iwp  ! offset of overlap derivative; can be defined with 11
  logical(kind=iwp) :: if_ssdm         = .false. ! State-dependent DM is used in Fock or not

  ! unit numbers
  integer(kind=iwp) :: LuPT2           = 0_iwp
  integer(kind=iwp) :: LuGAMMA         = 0_iwp
  integer(kind=iwp) :: LuCMOPT2        = 0_iwp
  integer(kind=iwp) :: LuSTD           = 0_iwp
  integer(kind=iwp) :: LuAPT2          = 0_iwp
  integer(kind=iwp) :: LuPT2GRD        = 0_iwp

  ! gradients and NAC switches
  logical(kind=iwp) :: do_grad         = .false.
  logical(kind=iwp) :: do_nac          = .false.
  logical(kind=iwp) :: do_csf          = .false. ! CSF term in deriv. coup.
  integer(kind=iwp) :: iRoot1          = 0_iwp
  integer(kind=iwp) :: iRoot2          = 0_iwp
  integer(kind=iwp) :: nStpGrd         = 1_iwp

  ! for removing the weired loop
  integer(kind=iwp) :: iStpGrd         = 1_iwp
  integer(kind=iwp) :: LUGRAD          = 0_iwp

  ! for IPEA
  logical(kind=iwp) :: do_lindep       = .false.
  logical(kind=iwp) :: if_invar        = .true. ! active invariance
  integer(kind=iwp) :: IDSAVGRD        = 0_iwp
  integer(kind=iwp) :: idBoriMat(8,13) = 0_iwp
  real(kind=wp)     :: ConvInvar       = 0.0_wp

  ! whether PT2 energy is invariant wrt rotations among inactive
  ! and secondary orbitals
  logical(kind=iwp) :: if_invaria      = .true.

  ! natural <-> quasi-canonical transformation of frozen orbitals
  real(kind=wp),allocatable :: TraFro(:)

  ! for mkfg3.f and derfg3.f compatibility
  ! number of CI vectors per batch
  integer(kind=iwp) :: nbuf1_grad      = 0_iwp
  ! number of tasks executed on the node (for parallel)
  integer(kind=iwp) :: nTasks_grad
  ! index of the executed tasks (for parallel)
  integer(kind=iwp),allocatable :: iTasks_grad(:)

end module caspt2_gradient
