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
! Copyright (C) 2014, Giovanni Li Manni                                *
!               2019, Oskar Weser                                      *
!***********************************************************************
module fcidump
  use fcidump_tables, only : OrbitalTable, FockTable, TwoElIntTable,&
    mma_allocate, mma_deallocate, fill_orbitals, fill_fock, fill_2ElInt
  use fcidump_transformations, only : get_orbital_E, fold_Fock
  use fcidump_reorder, only : reorder
  use fcidump_dump, only : dump_ascii, dump_hdf5
  implicit none
  private
  public :: make_fcidumps, transform, DumpOnly
  logical :: DumpOnly = .false.
  save
contains

  subroutine make_fcidumps(orbital_energies, folded_Fock, TUVX, core_energy, permutation)
    use rasscf_data, only : nacpar
    use general_data, only : nAsh
    implicit none
    real(8), intent(in) :: orbital_energies(:), folded_Fock(:), TUVX(:), core_energy
    integer, intent(in), optional :: permutation(:)
    type(OrbitalTable) :: orbital_table
    type(FockTable) :: fock_table
    type(TwoElIntTable) :: two_el_table

    call mma_allocate(fock_table, nacpar)
    call mma_allocate(two_el_table, size(TUVX))
    call mma_allocate(orbital_table, sum(nAsh))

    call fill_orbitals(orbital_table, orbital_energies)
    call fill_fock(fock_table, folded_Fock)
    call fill_2ElInt(two_el_table, TUVX)

    if (present(permutation)) then
      call reorder(orbital_table, fock_table, two_el_table, permutation)
    end if

    call dump_ascii(core_energy, orbital_table, fock_table, two_el_table)
    call dump_hdf5(core_energy, orbital_table, fock_table, two_el_table)

    call mma_deallocate(fock_table)
    call mma_deallocate(two_el_table)
    call mma_deallocate(orbital_table)
  end subroutine make_fcidumps



  subroutine transform(iter, CMO, DIAF, D1I_MO, F_IN, orbital_E, folded_Fock)
    implicit none
    integer, intent(in) :: iter
    real(8), intent(in) :: DIAF(:), CMO(:), D1I_MO(:)
    real(8), intent(inout) :: F_IN(:)
    real(8), intent(out) :: orbital_E(:), folded_Fock(:)
!
    call get_orbital_E(iter, DIAF, orbital_E)
    call fold_Fock(CMO, F_In, D1I_MO, folded_Fock)
  end subroutine transform
end module fcidump
