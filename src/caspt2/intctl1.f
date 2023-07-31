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
      SUBROUTINE INTCTL1(CMO)
      use caspt2_output, only:iPrGlb,debug
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "intgrl.fh"

      DIMENSION CMO(NCMO)

* Compute using conventional integral file:
      IF(IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' INTCTL1 calling TRACTL...'
        CALL XFLUSH(6)
      END IF
      Call TRACTL(0)
      CALL TRAONE(CMO)
      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(6,*)' INTCTL1 back from TRAONE.'
        CALL XFLUSH(6)
      END IF
c Compute FIMO, FAMO, ...  to workspace:
      CALL FOCK_RPT2

      RETURN
      END
