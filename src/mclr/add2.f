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
* Copyright (C) 1991, Anders Bernhardsson                              *
************************************************************************
      SubRoutine Add2(rMat,fact)
      use Arrays, only: SFock
*
*     Purpose:
*             Adds the contribution from the gradient to
*              [2]
*             E   Kappa. This is done to insure us about
*             a beautifull convergence of the PCG,
*             which is just the case if E is symmetric.
*
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"
#include "spin_mclr.fh"

#include "Input.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8 rMat(*)
      Real*8, Allocatable:: Temp(:)

      Do iS=1,nSym
        If (nOrb(is)*nOrb(is)==0) Cycle
        Call mma_allocate(Temp,nBas(is)**2,Label='Temp')
*
*    T=Brillouin matrix
*

        Call DGeSub(SFock(ipCM(is)),nOrb(is),'N',
     &              SFock(ipCM(is)),nOrb(is),'T',
     &              Temp,nOrb(is),
     &              nOrb(is),nOrb(is))
*
*               t           t
*   +1/2 { Kappa T - T kappa  }
*
*
        Call DaXpY_(nOrb(is)**2,-Four*Fact,Temp,1,
     &              rMat(ipMat(is,is)),1)
        Call mma_deallocate(Temp)
      End Do
      Return
      End
