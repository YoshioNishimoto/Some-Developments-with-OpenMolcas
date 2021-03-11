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
* Copyright (C) 1991, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE ADADST(IOBTP,IOBSM,IOBOFF,NIOB,JOBTP,JOBSM,JOBOFF,NJOB,
     &                  IJORD,ICLS,ISM,IGRP,KMIN,KMAX,I1,XI1S,NK,
     &                  NKDIM,IEND)
*
*
* Obtain mappings
* a+IORB a+ JORB !KSTR> = +/-!ISTR>
* Where IORB belongs to orbitals IOBTP,IOBSM
* and JORB belongs to JOBTP,JOBSM
* In the form
* I1(KSTR) =  ISTR if a+IORB a+ JORB !KSTR> = +/-!ISTR> , ISTR is in
* ICLS,ISM,IGRP.
* (numbering relative to TS start)
*
* Above +/- is stored in XI1S
* Number of K strings checked is returned in NK
* Only Kstrings with relative numbers from KMIN to KMAX are included
*
* If IEND .ne. 0 last string has been checked
*
* Jeppe Olsen , Winter of 1991
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
      INTEGER I1(NKDIM,*)
      DIMENSION XI1S(NKDIM,*)
*
      JGRP = IGRP + 1
      IF(IUNIQMP(JGRP).NE.JGRP) THEN
         JGRP = -IUNIQMP(JGRP)
C         write(6,*) ' Unique string group for mappings ',JGRP
      END IF
*. Are the creation arrays full or in compact form
*. N-1 => N
      IF(ISTAC(JGRP,1).NE.0.AND.ISTAC(JGRP,2).NE.0) THEN
         I1MPF = 1
         L1MP = NACOB
       ELSE
         I1MPF = 0
         L1MP = 0
       END IF
*
      KGRP = IGRP + 2
      IF(IUNIQMP(KGRP).NE.KGRP) THEN
        KGRP = -IUNIQMP(KGRP)
C        write(6,*) ' Unique string group for mappings ',KGRP
      END IF
*. N-2 => N-1
      IF(ISTAC(KGRP,1).NE.0.AND.ISTAC(KGRP,2).NE.0) THEN
         I2MPF = 1
         L2MP = NACOB
      ELSE
         I2MPF = 0
         L2MP = 0
      END IF
*
      CALL ADADS1(NK,I1,XI1S,IOBSM,IOBTP,IOBOFF,NIOB,
     &           JOBSM,JOBTP,JOBOFF,NJOB,IJORD,NKDIM,
     &           ICLS,ISM,iWORK(KSTSTM(KGRP,1)),
     &           iWORK(KSTSTM(KGRP,2)),I2MPF,L2MP,
     &           iWORK(KSTSTMI(KGRP)), iWORK(KSTSTMN(KGRP)),
     &           iWORK(KSTSTM(JGRP,1)),iWORK(KSTSTM(JGRP,2)),
     &           I1MPF,L1MP,
     &           iWORK(KSTSTMI(JGRP)), iWORK(KSTSTMN(JGRP)),
     &           iWORK(KEL1(IGRP)),
     &           iWORK(KEL3(IGRP)),iWORK(KEL1(IGRP+2)),
     &           iWORK(KEL3(IGRP+2)),iWORK(KISTSO(IGRP)),
     &           iWORK(KNSTSO(IGRP)),iWORK(KISTSO(IGRP+2)),
     &           iWORK(KNSTSO(IGRP+2)),NOCTYP(IGRP),NOCTYP(IGRP+2),
     &           NORB1,NORB2,NORB3,NACOB,KMAX,KMIN,IEND)
*
      RETURN
      END
