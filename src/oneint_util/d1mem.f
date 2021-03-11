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
      Subroutine D1Mem(nHer,MemD1,la,lb,lr)
      Implicit Real*8 (A-H,O-Z)
*
#include "itmax.fh"
#include "info.fh"
*
      nHer=mCentr
      MemD1 = 3*(la+1)*nHer +
     &        3*(lb+1)*nHer
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
