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
      SUBROUTINE CHO_WRRSTC(IPASS)
C
C     Purpose: write decomposition restart info for integral pass IPASS.
C
C     NB!!!  The restart files are assumed open on entry.
C
      use ChoArr, only: IntMap
      use ChoSwp, only: InfRed, InfVec
#include "implicit.fh"
#include "cholesky.fh"
#include "choorb.fh"

      PARAMETER (LSCR = 10)
      REAL*8  DSCR(LSCR)
      INTEGER JSCR(LSCR)

C     Start address on file.
C     ----------------------

      IADR = 0

C     Write molecular and configuration info.
C     ---------------------------------------

      IOPT = 1
      NWR  = 4
      JSCR(1) = NSYM
      JSCR(2) = NSHELL
      JSCR(3) = NNSHL
      JSCR(4) = 0
      CALL IDAFILE(LURST,IOPT,JSCR,NWR,IADR)

      IOPT = 1
      NWR  = NSYM
      DO ISYM = 1,NSYM
         JSCR(ISYM) = NBAS(ISYM)
      END DO
      CALL IDAFILE(LURST,IOPT,JSCR,NWR,IADR)

      IOPT = 1
      NWR  = 2
      IF (SCDIAG) THEN
         JSCR(1) = 1
      ELSE
         JSCR(1) = 0
      END IF
      JSCR(2) = CHO_ADRVEC
      CALL IDAFILE(LURST,IOPT,JSCR,NWR,IADR)

      IOPT = 1
      NWR  = 8
      DSCR(1) = THRCOM
      DSCR(2) = THRDIAG
      DSCR(3) = DAMP(1)
      DSCR(4) = DAMP(2)
      DSCR(5) = SPAN
      DSCR(6) = THRNEG
      DSCR(7) = WARNEG
      DSCR(8) = TOONEG
      CALL DDAFILE(LURST,IOPT,DSCR,NWR,IADR)

C     Write vector info.
C     ------------------

      IOPT = 1
      NWR  = 1
      JSCR(1) = IPASS
      CALL IDAFILE(LURST,IOPT,JSCR,NWR,IADR)

      IOPT = 1
      CALL IDAFILE(LURST,IOPT,INFRED,IPASS,IADR)

      DO ISYM = 1,NSYM
         IOPT = 1
         NWR  = 1
         JSCR(1) = NUMCHO(ISYM)
         CALL IDAFILE(LURST,IOPT,JSCR,NWR,IADR)
         IF (NUMCHO(ISYM) .GT. 0) THEN
            DO J = 1,SIZE(INFVEC,2)
               IOPT = 1
               NTOT = NUMCHO(ISYM)
               CALL IDAFILE(LURST,IOPT,InfVec(:,J,ISYM),NTOT,IADR)
            END DO
         END IF
      END DO

C     Write integral shell pair map to disk.
C     --------------------------------------

      IOPT = 1
      NDIM=0
      IF (Allocated(IntMap)) THEN
         NDIM = SIZE(IntMap)
         JADR = 0
         CALL IDAFILE(LUMAP,IOPT,INTMAP,NDIM,JADR)
      END IF

      END
