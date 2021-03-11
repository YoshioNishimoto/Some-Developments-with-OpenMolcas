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
      FUNCTION GTIJKL_MCLR(I,J,K,L)
*
* Obtain  integral (I J ! K L )
* where I,J,K and l refers to active orbitals in
* Type ordering
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "detdim.fh"
#include "WrkSpc.fh"
#include "glbbas_mclr.fh"

#include "Input.fh"
#include "orbinp_mclr.fh"
#include "crun_mclr.fh"
#include "Pointers.fh"
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      IABS = IREOTS(I)
*
      JABS = IREOTS(J)
*
      KABS = IREOTS(K)
*
      LABS = IREOTS(L)
*
      IJ=itri(iABS,JABS)
      KL=itri(kABS,lABS)

      GTIJKL_MCLR = WORK(K2INT+itri(IJ,KL)-1)
*
      RETURN
      END
