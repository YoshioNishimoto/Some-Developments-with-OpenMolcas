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
      SUBROUTINE MEMCHK
*
* Check memory allocated  with the memory manager
*
      Call Getmem('FLUSH','CHECK','REAL',idum,idum)
C$$$      CALL MEMMAN(IDUM,IDUM,'CHECK ',IDUM,'IDUM  ')
*
      RETURN
      END