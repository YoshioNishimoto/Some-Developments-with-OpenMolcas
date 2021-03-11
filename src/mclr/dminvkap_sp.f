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
      SubRoutine DMInvKap_sp(rMFact,rin,rout,isym)
*
*     _____     -1
*     Kappa  = M  Kappa
*          ip   pq     iq
*
*
*     In: rMFact        Factorized preconditioner (diagonal part
*                         of the electronic hessian that couple
*                         rotations with one common index)
*     In,Out rOut       Orbital rotaotion
*
*     iSym              Symmetry of rotation
*
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Real*8 rMFact(*),rin(*),rout(*)
*
      Call GetMem('Temp','ALLO','REAL',ipT,ndens2)
      Call Uncompress(rin,Work(ipT),isym)
*
      Call Compress(Work(ipT),rout,isym)
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(rMFact)
      end
