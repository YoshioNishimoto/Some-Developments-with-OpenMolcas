************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Giovanni Li Manni                                      *
************************************************************************
      Subroutine Readinp_expbas()
c
c     Author: G. Li Manni (University of Geneva)
c
      Implicit Real*8(a-h,o-z)
#include "info_expbas.fh"
      Character*180  Line, Blank, key, Get_Ln
      External Get_Ln
c
c Initial values
c
      DoExpbas = .true.
      DoDesy   = .false.
      DoSort   = .false.
      EB_FileOrb  = ' '
c
      LuSpool=18
      LuSpool=isFreeUnit(LuSpool)
      Call SpoolInp(LuSpool)

      Rewind(LuSpool)
      Call RdNLst(LuSpool,'EXPBAS')
      Blank=' '

  999 Continue
*      Read(LuSpool,'(A)',End=9940) Line
      key =Get_Ln(LuSpool)
      Call LeftAd(key)
      Line = key
      If (Line(1:1).eq.'*' ) Goto 999
      If (Line.eq.Blank ) Goto 999
      Call UpCase(Line)
      If (Line(1:4).eq.'NOEX') then
        DoExpbas = .false.
        Go To 999
      end if
      If (Line(1:4).eq.'DESY') then
        DoDesy   = .true.
        Go To 999
      end if
      If (Line(1:4).eq.'SORT') then
        DoSort   = .true.
        Go To 999
      end if
      If (Line(1:4).eq.'FILE') then
        Line=Get_Ln(LuSpool)
        Call FileOrb(Line,EB_FileOrb)
        Go To 999
      end if
      If (Line(1:4).eq.'END ') Go To 99999
      Write (6,*) 'Unidentified key word  : '
      Call FindErrorLine
      Call Quit_OnUserError()
c
c END of Input
c
      WRITE(6,*)' READIN: Premature end of file when reading selected'
      CALL ABEND()

99999 Continue
      End
