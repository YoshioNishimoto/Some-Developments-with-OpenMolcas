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
      SUBROUTINE CHO_SETMAXSHL(DIAG,DIASH,ISYSH,IRED)
C
C     Purpose: set max. shell pair data for selection procedure.
C
      use ChoArr, only: iSP2F, iAtomShl
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRed
#include "implicit.fh"
      REAL*8 DIAG(*), DIASH(*)
      INTEGER   ISYSH(*)
#include "real.fh"
#include "cholesky.fh"

      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_SETMAXSHL')

C     Initialize the largest diagonal in each shell pair.
C     ---------------------------------------------------

      DIASH(1:NNSHL)=Zero
      ISYSH(1:NNSHL)=0

C     Find largest diagonal in each shell pair. Loop only
C     over those that are included in the reduced set at hand.
C     --------------------------------------------------------

      IF (IRED .EQ. 1) THEN
         DO ISYMAB = 1,NSYM
            DO ISHLAB = 1,NNSHL
               IAB1 = IIBSTR(ISYMAB,IRED)
     &              + IIBSTRSH(ISYMAB,ISHLAB,IRED) + 1
               IAB2 = IAB1 + NNBSTRSH(ISYMAB,ISHLAB,IRED) - 1
               DO IAB = IAB1,IAB2
                  DIASH(ISHLAB) = Max(DIASH(ISHLAB),DIAG(IAB))
                  If (DIASH(ISHLAB)==DIAG(IAB)) ISYSH(ISHLAB)=ISYMAB
               END DO
            END DO
         END DO
      ELSE IF ((IRED.EQ.2) .OR. (IRED.EQ.3)) THEN
         DO ISYMAB = 1,NSYM
            DO ISHLAB = 1,NNSHL
               JAB1 = IIBSTR(ISYMAB,IRED)
     &              + IIBSTRSH(ISYMAB,ISHLAB,IRED) + 1
               JAB2 = JAB1 + NNBSTRSH(ISYMAB,ISHLAB,IRED) - 1
               DO JAB = JAB1,JAB2
                  IAB = INDRED(JAB,IRED)
                  DIASH(ISHLAB) = Max(DIASH(ISHLAB),DIAG(IAB))
                  If (DIASH(ISHLAB)==DIAG(IAB)) ISYSH(ISHLAB)=ISYMAB
               END DO
            END DO
         END DO
      ELSE
         WRITE(LUPRI,*) SECNAM,': unknown reduced set, IRED = ',IRED
         CALL CHO_QUIT('Unknown reduced set in '//SECNAM,104)
      END IF

C     Exclude 2-center diagonals (if requested).
C     The effect of this is that 2-center diagonals can never be
C     qualified; they may still be included in the vectors, though.
C     If CHO_NO2CENTER=T, the 2-center diagonals are removed from the
C     initial diagonal and we need not worry here.
C     ---------------------------------------------------------------

      IF (CHO_1CENTER .AND. .NOT.CHO_NO2CENTER) THEN
         DO ISAB = 1,NNSHL
            ISHLAB = ISP2F(ISAB)
            CALL CHO_INVPCK(ISHLAB,ISHLA,ISHLB,.TRUE.)
            IF (IATOMSHL(ISHLA) .NE. IATOMSHL(ISHLB)) THEN
               DIASH(ISAB) = 0.0D0
            END IF
         END DO
      END IF

      END
