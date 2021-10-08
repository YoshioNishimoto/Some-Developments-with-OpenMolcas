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
!
! OpenMolcas wrapper for the soMPSoo program
! rdinput reads input from the OpenMolcas input file
!
!*******************

subroutine rdinput(refwfnfile)

! set global variables directly from the soMPSoo program

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(out) :: refwfnfile
character(len=180) :: Line, key
character(len=9001) :: dline
character(len=9001) :: frozen_str
integer(kind=iwp) :: LuSpool, iError, i, isplit
integer(kind=iwp), external :: isFreeUnit
logical(kind=iwp), external :: next_non_comment
character(len=180), external :: Get_Ln

! Initial values

refwfnfile = 'JOBIPH'

LuSpool = isFreeUnit(18)
call SpoolInp(LuSpool)

rewind(LuSpool)
call RdNLst(LuSpool,'soMPSoo')

!> position is now at &soMPSoo

do
  key = Get_Ln(LuSpool)
  call LeftAd(key)
  Line = key
  if (Line(1:1) == '*') cycle
  if (Line == ' ') cycle
  call UpCase(Line)
  select case (Line(1:4))
    case ('NOPC')
      !========= NOPC =============
      !no_pc = .true.
      write(u6,*) 'YEAH here we go'

    case ('END ')
      exit

    case default
      write(u6,*) 'Unidentified key word  : '
      call FindErrorLine()
      call Quit_OnUserError()

  end select

end do
! END of Input

contains

subroutine error(code)

  integer(kind=iwp), intent(in) :: code

  if (code == 1) call WarningMessage(2,'Premature end of input file.')
  call WarningMessage(2,'Read error during input preprocessing.')
  call Quit_OnUserError()
  call ABEND()

end subroutine

end subroutine rdinput
