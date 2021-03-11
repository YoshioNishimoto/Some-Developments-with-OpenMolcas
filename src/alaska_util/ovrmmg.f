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
      Subroutine OvrMmG(nHer,MmOvrG,la,lb,lr)
*
      nHer=(la+lb+1+2)/2
      MmOvrG = 3*nHer*(la+2) +
     &        3*nHer*(lb+2) +
     &        3*nHer +
     &        3*(la+2)*(lb+2) + 2
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lr)
      End
