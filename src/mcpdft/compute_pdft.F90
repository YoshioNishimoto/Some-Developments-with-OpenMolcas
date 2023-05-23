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
! Copyright (C) 2022, Matthew R Hennefarth                             *
!***********************************************************************

module compute_pdft
  use definitions, only: wp
  use mcpdft_output, only: lf, iPrGlb, debug
  implicit none
  private

  public :: load_1e_integrals, get_total_charge

  contains
  subroutine get_total_charge(total_charge)
    ! Computes the total charge of the system. This should be charge of
    ! nuclei - total number of electrons. Obviously, it is 0 for neutral
    ! systems. Always an integer.
    use OneDat, only: sNoOri

    ! For nTot1, nSym
#include "rasdim.fh"
#include "general.fh"

    integer, intent(out) :: total_charge

    real(kind=wp), dimension(ntot1+4) :: overlap
    real(kind=wp) :: tot_nuc_charge, tot_el_charge = 0.0
    integer :: iSym

    integer :: irc = -1
    integer :: iOpt = ibset(0, sNoOri)
    integer :: iComp = 1
    integer :: iSyLbl = 1
    character*8 :: label = 'Mltpl  0'

    call rdone(irc, iopt, label, icomp, overlap, isylbl)
    if (irc /= 0) then
      call rdone_error(label, irc)
    end if

    tot_nuc_charge = overlap(ntot1+4)
    do iSym = 1, nSym
      tot_el_charge = tot_el_charge - 2.0D0*DBLE(nFro(iSym)+nIsh(iSym))
    end do
    tot_el_charge = tot_el_charge-DBLE(nActEl)

    total_charge = int(tot_nuc_charge + tot_el_charge)
  end subroutine

  subroutine load_1e_integrals(h1, kincore, nucelcore)
    ! Loads the one-electron integrals in the AO basis. That is the
    ! terms of h_pq but actually h_{mu, nu}
    !
    ! Returns:
    !   h1: ndarray of length nTot1 (1/2 num of AO Basis Functions)
    !     Elements of h_{mu, nu} in a packed lower triangular array
    !
    !   kincore: ndarray of length nTot1
    !     Electron kinetic energy integrals in AO basis stored in a
    !     lower triangular array
    !
    !   nucelcore: ndarray of length nTot1
    !     Nuclear-electron attraction integrals in AO basis stored in a
    !     lower triangular array

    use OneDat, only: sNoNuc, sNoOri

    ! For nTot1
#include "rasdim.fh"
#include "general.fh"

    real(kind=wp), dimension(nTot1), intent(out) :: h1, kincore, nucelcore

    integer :: iComp = 1
    integer :: iSyLbl = 1
    integer :: iRc = -1
    integer :: iOpt = ibset(ibset(0, sNoOri), sNoNuc)
    character*8 :: Label = 'OneHam  '

    call RdOne(iRc, iOpt, Label, iComp, H1, iSyLbl)
    if (iRc /= 0) then
      call rdone_error(label, irc)
    end if

    irc = -1
    label = 'Kinetic '
    call rdone(irc, iOpt, label, iComp, kincore, iSyLbl)
    if (iRc /= 0) then
      call rdone_error(label, irc)
    end if

    irc = -1
    label = 'Attract '
    call rdone(irc, iopt, label, iComp, nucelcore, iSylbl)
    if (iRc /= 0) then
      call rdone_error(label, irc)
    end if
  end subroutine

  subroutine rdone_error(label, irc)
    character*8, intent(in) :: label
    integer , intent(in) :: irc

    write(lf,*) 'MC-PDFT: iRc from Call RdOne not 0'
    write(lf,*) 'Label = ', Label
    write(lf,*) 'iRc = ', iRc
    call abend
  end subroutine rdone_error

end module compute_pdft