************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine cpinp_(LUnit,iRc)
      implicit integer (a-z)
      character*180 line
      character*1 ch
#include "warnings.fh"
      iRc=_RC_ALL_IS_WELL_
* The following code will open, and return the unit number LUSpool,
* of an ASCII file with a copy of the presently used input.
* The records are strings of 180 characters, conforming to the
* standards set e.g. in src/util/inputil.f. This may change in the
* future, so look out!

      Call SpoolInp(LUSpool)
      Call Disable_Spool()
      Rewind(LUSpool)
* Now open a new file, and copy only the input between the '&RASSCF'
* and the 'End of Input' markers, inclusive, and skipping any commented
* lines. (The latter is put there by sbin/auto.plx, so it is safe to
* assume it is not abbreviated). The copied lines are left adjusted.
* Positioning LUSpool after the '&RASSCF' marker.
      Call RdNLst(LuSpool,'MCPDFT')
* Opening a new file:
      LUnit=99
      LUnit=IsFreeUnit(LUnit)
      Call Molcas_Open (LUnit,'CleanInput')
* Copy only the relevant lines of input:
      line=' '
      line(1:7)='&MCPDFT'
      write(LUnit,'(A180)') line
  10  continue
      read(luspool,'(A180)',err=9910,end=9910) line
      call leftad(line)
      ch=line(1:1)
      if(ch.ne.' ' .and. ch.ne.'*' .and. ch.ne.'!') then
       write(LUnit,'(A180)') line
      end if
      call upcase(line(1:12))
      if (line(1:12).ne.'END OF INPUT') goto 10
      call close_luspool(LUSpool)
      return
 9910 continue
* Something went wrong...Let the caller handle it:
      iRc=_RC_INPUT_ERROR_
      return
      end
