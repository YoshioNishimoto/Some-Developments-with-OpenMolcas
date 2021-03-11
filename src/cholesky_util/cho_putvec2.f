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
      SUBROUTINE CHO_PUTVEC2(CHOVEC,NUMVEC,IVEC1,ISYM)
C
C     Purpose: write Cholesky vectors IVEC=IVEC1,...,IVEC1+NUMVEC-1
C              of symmetry ISYM to file.
C
C     Version 2: handles several reduced set at a time.
C
#include "implicit.fh"
      DIMENSION CHOVEC(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      external ddot_

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_PUTVEC2')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (N2 = INFVEC_N2)

      INFVEC(I,J,K)=IWORK(ip_INFVEC-1+MAXVEC*N2*(K-1)+MAXVEC*(J-1)+I)
      NDIMRS(I,J)=IWORK(ip_NDIMRS-1+NSYM*(J-1)+I)

#if defined (_DEBUG_)
      CALL QENTER('_PUTVEC2')
#endif

C     Return if no vectors.
C     ---------------------

      IF (NUMVEC .LT. 1) THEN
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': WARNING: no vectors in this call!'
            WRITE(LUPRI,*) SECNAM,': NUMVEC = ',NUMVEC
         END IF
         GO TO 1 ! exit
      END IF

C     Check symmetry.
C     ---------------

      IF ((ISYM.LT.1) .OR. (ISYM.GT.NSYM)) THEN
         WRITE(LUPRI,*) SECNAM,': symmetry out of bounds'
         WRITE(LUPRI,*) 'ISYM = ',ISYM
         CALL CHO_QUIT('Symmetry out of bounds in '//SECNAM,104)
      END IF

C     Check vector index.
C     -------------------

      IVEC2 = IVEC1 + NUMVEC - 1
      IF ((IVEC1.LT.1) .OR. (IVEC1.GT.MAXVEC) .OR.
     &    (IVEC2.LT.1) .OR. (IVEC2.GT.MAXVEC)) THEN
         WRITE(LUPRI,*) SECNAM,': vector index out of bounds'
         WRITE(LUPRI,*) 'IVEC1 = ',IVEC1,' IVEC2 = ',IVEC2
         WRITE(LUPRI,*) '...must be between 1 and ',MAXVEC
         CALL CHO_QUIT('Vector index out of bounds in '
     &                 //SECNAM,104)
      END IF

C     Check for overflow for WA file addressing.
C     ------------------------------------------

      IF (CHO_ADRVEC .EQ. 1) THEN
         IADR2 = INFVEC(IVEC2,4,ISYM)
         IF (INFVEC(IVEC1,4,ISYM) .LT. 0) THEN
            WRITE(LUPRI,*) 'Error in ',SECNAM,':'
            WRITE(LUPRI,*) 'Illegal disk address for first vector: ',
     &                     INFVEC(IVEC1,4,ISYM)
            IF (INFVEC(IVEC1,4,ISYM) .LT. -1) THEN
               WRITE(LUPRI,*) '....is it an overflow?'
            END IF
            WRITE(LUPRI,*) 'IVEC1 = ',IVEC1,' ISYM = ',ISYM
            CALL CHO_QUIT('Illegal disk address in '//SECNAM,104)
         ELSE IF (IADR2 .LT. INFVEC(IVEC1,4,ISYM)) THEN
            WRITE(LUPRI,*) 'Error in ',SECNAM,':'
            WRITE(LUPRI,*) 'Illegal disk address for last vector: ',
     &                     IADR2
            IF (IADR2 .LT. -1) THEN
               WRITE(LUPRI,*) '....is it an overflow?'
            END IF
            WRITE(LUPRI,*) 'IVEC2 = ',IVEC2,' ISYM = ',ISYM
            CALL CHO_QUIT('Illegal disk address in '//SECNAM,104)
         END IF
      END IF

C     Call the low-level I/O routines.
C     CHO_ADRVEC=1: WA files.
C     CHO_ADRVEC=2: DA files.
C     Set (next) disk addresses.
C     --------------------------------


      IF (CHO_ADRVEC .EQ. 1) THEN
         IOPT = 1  ! synchronous write option
         LTOT = 0
         DO IVEC = 1,NUMVEC
            JVEC = IVEC1 + IVEC - 1
            JRED = INFVEC(JVEC,2,ISYM)
            LTOT = LTOT + NDIMRS(ISYM,JRED)
         END DO
         IADR = INFVEC(IVEC1,3,ISYM)
         CALL DDAFILE(LUCHO(ISYM),IOPT,CHOVEC,LTOT,IADR)
         DO IVEC = 1,NUMVEC-1
            JVEC = IVEC1 + IVEC - 1
            JRED = INFVEC(JVEC,2,ISYM)
            IWORK(ip_INFVEC+MAXVEC*N2*(ISYM-1)+MAXVEC*2+JVEC) =
     &      INFVEC(JVEC,3,ISYM) + NDIMRS(ISYM,JRED)
         END DO
         IVEC = NUMVEC
         JVEC = IVEC1 + IVEC - 1
         IF (JVEC .LT. MAXVEC) THEN
            JRED = INFVEC(JVEC,2,ISYM)
            IWORK(ip_INFVEC+MAXVEC*N2*(ISYM-1)+MAXVEC*2+JVEC) =
     &      INFVEC(JVEC,3,ISYM) + NDIMRS(ISYM,JRED)
         END IF
      ELSE IF (CHO_ADRVEC .EQ. 2) THEN
         IOPT = 1  ! synchronous write option
         KOFF = 1
         DO IVEC = 1,NUMVEC-1
            JVEC = IVEC1 + IVEC - 1
            JRED = INFVEC(JVEC,2,ISYM)
            LTOT = NDIMRS(ISYM,JRED)
            IADR = INFVEC(JVEC,3,ISYM)
            CALL DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(KOFF),LTOT,IADR)
            IWORK(ip_INFVEC+MAXVEC*N2*(ISYM-1)+MAXVEC*2+JVEC) = IADR
            KOFF = KOFF + LTOT
         END DO
         IVEC = NUMVEC
         JVEC = IVEC1 + IVEC - 1
         JRED = INFVEC(JVEC,2,ISYM)
         LTOT = NDIMRS(ISYM,JRED)
         IADR = INFVEC(JVEC,3,ISYM)
         CALL DDAFILE(LUCHO(ISYM),IOPT,CHOVEC(KOFF),LTOT,IADR)
         IF (JVEC .LT. MAXVEC) THEN
            IWORK(ip_INFVEC+MAXVEC*N2*(ISYM-1)+MAXVEC*2+JVEC) = IADR
         END IF
      ELSE
         CALL CHO_QUIT('CHO_ADRVEC out of bounds in '//SECNAM,102)
      END IF

C     Debug stuff.
C     ------------

      IF (LOCDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) 'Vectors ',IVEC1,' to ',IVEC1+NUMVEC-1,
     &                  ' of symmetry ',ISYM,' written to unit ',
     &                  LUCHO(ISYM)
         KOFF = 1
         DO IVEC = 1,NUMVEC
            JVEC = IVEC1 + IVEC - 1
            JRED = INFVEC(JVEC,2,ISYM)
            LVEC = NDIMRS(ISYM,JRED)
            JADR = INFVEC(JVEC,3,ISYM)
            XNRM = SQRT(DDOT_(LVEC,CHOVEC(KOFF),1,CHOVEC(KOFF),1))
            WRITE(LUPRI,*) 'Vector:',JVEC,' address: ',JADR,' norm: ',
     &                     XNRM,' dimension: ',LVEC
            KOFF = KOFF + LVEC
         END DO
      END IF

    1 CONTINUE
#if defined (_DEBUG_)
      CALL QEXIT('_PUTVEC2')
#endif

      END
