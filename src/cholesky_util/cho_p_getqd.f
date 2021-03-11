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
      SubRoutine Cho_P_GetQD(QD)
C
C     Purpose: copy qualified diagonal elements from global diagonal to
C              array QD.
C
      Implicit None
      Real*8 QD(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "choglob.fh"
#include "WrkSpc.fh"

      Integer kQD, kD, iSym, iQ, iAB
      Integer IndRed_G, iQuAB, i, j

      IndRed_G(i,j)=iWork(ip_IndRed_G-1+mmBstRT_G*(j-1)+i)
      iQuAB(i,j)=iWork(ip_iQuAB-1+MaxQual*(j-1)+i)

      kQD = 0
      kD = ip_Diag_G - 1
      Do iSym = 1,nSym
         Do iQ = 1,nQual(iSym)
            iAB = IndRed_G(iQuAB(iQ,iSym),2)
            QD(kQD+iQ) = Work(kD+iAB)
         End Do
         kQD = kQD + nQual(iSym)
      End Do

      End
