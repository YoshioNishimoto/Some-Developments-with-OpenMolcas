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
! Copyright (C) 2023, Matthew R. Hennefarth                             *
!***********************************************************************

subroutine ao2mo_1e(CMO, Mat)
  ! Convert a 1-electron matrix from AO to MO basis. Writes over the
  ! input matrix.
  !
  ! Args:
  !   CMO: ndarray (ntot2)
  !     MO coefficients
  !
  !   Mat: ndarray(ntot2)
  !     1-electron matrix in AO basis
  !
  ! Returns:
  !   Mat: ndarray(ntot2)
  !     1-electron Matrix in MO basis
  use definitions, only: wp
  use stdalloc, only: mma_allocate, mma_deallocate
  implicit none

  ! For ntot2, nsym, nbas, nOrb, nfro
#include "rasdim.fh"
#include "general.fh"
  real(kind=wp), dimension(ntot2), intent(in) :: CMO
  real(kind=wp), dimension(ntot2), intent(inout) :: Mat

  real(kind=wp), parameter :: one=1.0D0, zero=0.0D0

  integer :: iSym, iBasis, iOrbitals, iFrozen
  integer, dimension(3) :: offset
  real(kind=wp), dimension(:), allocatable :: scr1, scr2

  offset = 1
  do iSym = 1, nSym
    iBasis = nBas(iSym)
    iOrbitals = nOrb(iSym)
    if (iBasis == 0 .or. iOrbitals == 0) then
      cycle
    end if
    iFrozen = nFro(iSym)

    call mma_allocate(scr1,iBasis*iBasis, "Scr1")
    call mma_allocate(scr2,iOrbitals*iBasis, "Scr2")

    call square(Mat(offset(1)), scr1, 1, iBasis, iBasis)

    call dgemm_('N','N', iBasis, iOrbitals, iBasis, one, scr1, iBasis, &
            CMO(offset(2)+(iFrozen*iBasis)), iBasis, zero, scr2, iBasis)
    call dgemm_tri('T','N', iOrbitals, iOrbitals, iBasis, one, scr2, &
            iBasis, CMO(offset(2)+(iFrozen*iBasis)), iBasis, zero, &
            Mat(offset(3)), iOrbitals)

    call mma_deallocate(scr1)
    call mma_deallocate(scr2)

    offset(1) = offset(1) + (iBasis*iBasis+iBasis)/2
    offset(2) = offset(2) + iBasis*iBasis
    offset(3) = offset(3) + (iOrbitals*iOrbitals+iOrbitals)/2
  end do
end subroutine ao2mo_1e