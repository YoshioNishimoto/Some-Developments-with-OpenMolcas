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
! Copyright (C) 2021, Stefan Knecht                                    *
!***********************************************************************

subroutine sompsoo_(iReturn)



!> soMPSoo stuff
use orbopt_header, only: print_orbopt_header    
!> OpenMOLCAS stuff    
use Definitions, only: iwp
use Para_Info, only: King
use Definitions, only: iwp, u6, r4

implicit none

integer(kind=iwp), intent(inout) :: iReturn
character(len=256) :: refwfnfile
logical(kind=iwp), parameter :: rel_ham = .false., from_molcas = .true.
! ----------------------------------------------------------------------

!> make sure that in an MPI parallel setting, soMPSoo is called exclusively by the KING
if (KING()) then
    !> initialize
    refwfnfile = ''
    !> print soMPSoo header
    call print_orbopt_header(u6)
    call xflush(u6)
    !> read and process input
    call rdinput(refwfnfile)
    call xflush(u6)
    write(u6,*) ' I was in soMPSoo ...'
endif
iReturn = 0

!!> set DMRG driver as active space solver
!call set_as_solver()

!!> read DMRG settings (driver-specific input)
!call set_dmrg_settings()

!> call wave function optimizer

!call rasscf(iReturn)

!#ifdef _DMRG_
!!> reset in case we call RASSCF (or RASSI or CASPT2) afterwards requesting a CI driver
!if (doDMRG) doDMRG = .false.
!#endif

end subroutine sompsoo_
! ----------------------------------------------------------------------
