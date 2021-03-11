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
* Copyright (C) 1991,1994, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE ADST(IORB,NORB,ICLS,ISM,IGRP,KMIN,KMAX,I1,XI1S,LI1,
     &                NK,IEND)
*
*
* Obtain mappings
* a+JORB !KSTR> = +/-!ISTR> for orbitals IORB - IORB+NORB-1
*. All orbitals are assumed to belong to the same TS class
* In the form
* I1(KSTR,JORB) =  ISTR if A+JORB !KSTR> = +/-!ISTR> , ISTR is in
* ICLS,ISM,IGRP.
* (numbering relative to TS start)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
* If KMAX = -1, all strings are included
* If all K strings have been searched IEND is set to 1
*
* Jeppe Olsen , Winter of 1991
*               January 1994 : modified to allow for several orbitals
*
* ======
*. Input
* ======
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "detdim.fh"
#include "WrkSpc.fh"
#include "orbinp_mclr.fh"
#include "strinp_mclr.fh"
#include "stinf_mclr.fh"
#include "strbas_mclr.fh"
*
* =======
*. Output
* =======
*
      INTEGER I1(*)
      DIMENSION XI1S(*)
*. Type of mapping
C?     write(6,*) ' ADST: IGRP IORB = ',IGRP, IORB
*
      JGRP = IGRP + 1
      IF(IUNIQMP(JGRP).NE.JGRP) THEN
        JGRP = -IUNIQMP(JGRP)
C        write(6,*) ' Unique string group for mappings ',JGRP
      END IF
*
      IF(ISTAC(JGRP,1).NE.0.AND.ISTAC(JGRP,2).NE.0) THEN
*. Full Map
        IMPF = 1
        LMAP = NACOB
      ELSE
*. Only creation map
        IMPF = 0
        LMAP = -1
      END IF
*
      CALL ADS1(NK,I1,XI1S,LI1,IORB,NORB,
     &           ICLS,ISM,iWORK(KSTSTM(JGRP,1)),
     &           iWORK(KSTSTM(JGRP,2)),iWORK(KSTSTMN(JGRP)),
     &           iWORK(KSTSTMI(JGRP)),IMPF,LMAP,iWORK(KEL1(IGRP)),
     &           iWORK(KEL3(IGRP)),iWORK(KEL1(IGRP+1)),
     &           iWORK(KEL3(IGRP+1)),iWORK(KISTSO(IGRP)),
     &           iWORK(KNSTSO(IGRP)),iWORK(KISTSO(IGRP+1)),
     &           iWORK(KNSTSO(IGRP+1)),NOCTYP(IGRP),NOCTYP(IGRP+1),
     &           NORB1,NORB2,NORB3,ISMFTO,NACOB,KMAX,KMIN,IEND)
*
      RETURN
      END
