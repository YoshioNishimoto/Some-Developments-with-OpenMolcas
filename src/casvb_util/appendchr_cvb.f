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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine appendchr_cvb(c,character,iskip)
      implicit REAL*8 (a-h,o-z)
      character*(*)c,character

      ibegin=len_trim_cvb(c)+1+iskip
      iend=min(len(c),ibegin+len_trim_cvb(character)-1)
      c(ibegin:iend)=character(1:len_trim_cvb(character))
      return
      end
