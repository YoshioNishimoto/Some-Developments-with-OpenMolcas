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
* Copyright (C) 2010, Jonas Bostrom                                    *
************************************************************************
      SubRoutine ChoMP2g_AmpDiag(irc,Diag,EOcc,EVir)
C
C     Jonas Bostrom, Jan. 2010.
C
C     Purpose: Construct diagonal for decomposition of amplitude
C              vectors.
C
#include "implicit.fh"
      Integer irc
      Real*8  Diag(*), EOcc(*), EVir(*)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2g.fh"

      Character(LEN=7), Parameter:: ThisNm = 'AmpDiag'
      Character(LEN=15), Parameter:: SecNam = 'ChoMP2g_AmpDiag'

      MulD2h(k,l)=iEor(k-1,l-1)+1

      irc = 0

C     Initialization.
C     ------------------------

      iVecType = 6
      kD0 = 0

C     Construct Diagonal
C     ------------------------

      Do iSym = 1, nSym
            Do iSymI = 1,nSym
               iSymA = MulD2h(iSymI,iSym)
               kD1 = kD0 + iMoMo(iSymA,iSymI,iVecType)
               Do iI = 1, nOcc(iSymI)
                  kD2 = kD1 + nVir(iSymA)*(iI-1)
                  Ei = EOcc(iOcc(iSymI) +iI)
                  Do iA = 1, nVir(iSymA)
                     iAI = kD2 + iA
                     DE = 2.0d0*(EVir(iVir(iSymA)+iA)-Ei)
                     Diag(iAI) = Diag(iAI)/DE
                  End Do
               End Do
            End Do
            kD0 = kD0 + nMoMo(iSym,iVecType)
         End Do

      Return
      End
