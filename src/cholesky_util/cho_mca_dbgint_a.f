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
      SUBROUTINE CHO_MCA_DBGINT_A()
C
C     Purpose: regenerate and check all integrals (or the number
C              of columns specified in input).
C
C     NOTE: this is *not* meant for production calculations, only for
C           debugging, as:
C           1) full Cholesky vectors are read
C           2) calculations are performed in full (no use of red. sets
C              apart from first)
C           3) full integral symmetry not used
C              (only partial particle permutation symmetry)
C
      use ChoArr, only: nBstSh, iSP2F
#include "implicit.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      CHARACTER(LEN=16), PARAMETER:: SECNAM = 'CHO_MCA_DBGINT_A'

      Real*8 XLBAS(8)

      Real*8, Allocatable:: Int1(:), Wrk(:)

      MULD2H(I,J)=IEOR(I-1,J-1)+1

C     Force computation of full shell quadruple.
C     ------------------------------------------

      IF (IFCSEW .NE. 1) THEN
         WRITE(LUPRI,*) SECNAM,': WARNING: resetting IFCSEW from ',
     &                  IFCSEW,' to 1.'
         WRITE(LUPRI,*) SECNAM,
     &   ': memory demands are significantly increased by this!'
         IFCSEW = 1
      END IF

C     Initializations.
C     ----------------

      GLMAX = 0.0D0
      GLMIN = 1.0D15
      GLRMS = 0.0D0
      XTCMP = 0.0D0
      XPECT = 0.0D0

C     Make first reduced set the current reduced set.
C     -----------------------------------------------

      CALL CHO_RSCOPY(1,2)

C     Allocate memory for largest integral quadruple.
C     -----------------------------------------------

      LINTMX = MX2SH*MX2SH
      Call mma_allocate(INT1,LINTMX,Label='INT1')

C     Allocate max. memory
C     ----------------------------------------------------------

      Call mma_maxDBLE(LWRK)
      Call mma_allocate(WRK,LWRK/2,Label='WRK')
      CALL XSETMEM_INTS(LWRK/2)

C     Print header.
C     -------------

      CALL CHO_HEAD('Integral Error Analysis','=',80,LUPRI)
      WRITE(LUPRI,'(/,A,/,A)')
     & '    C     D     A     B   Abs. Min.    Abs. Max.      RMS',
     & '--------------------------------------------------------------'

C     Loop over shell quadruples.
C     ---------------------------

      ISAB1 = 1
      IF (NCOL_CHK .GT. 0) THEN
         ISAB2 = MIN(NCOL_CHK,NNSHL)
      ELSE
         ISAB2 = NNSHL
      END IF

      DO ISHLAB = ISAB1,ISAB2

         CALL CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.TRUE.)
         IF (ISHLB .EQ. ISHLA) THEN
            NUMAB = NBSTSH(ISHLA)*(NBSTSH(ISHLA) + 1)/2
         ELSE
            NUMAB = NBSTSH(ISHLA)*NBSTSH(ISHLB)
         END IF

         IF (NCOL_CHK .GT. 0) THEN
            ISCD1 = 1
            ISCD2 = NNSHL
         ELSE
            ISCD1 = ISHLAB
            ISCD2 = NNSHL
         END IF

         DO ISHLCD = ISCD1,ISCD2

            CALL CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.TRUE.)
            IF (ISHLD .EQ. ISHLC) THEN
               NUMCD = NBSTSH(ISHLC)*(NBSTSH(ISHLC) + 1)/2
            ELSE
               NUMCD = NBSTSH(ISHLC)*NBSTSH(ISHLD)
            END IF
            LINT1 = NUMCD*NUMAB ! actual space needed for (CD|AB)

C           Compute expected number of comparisons.
C           ---------------------------------------

            XPECTL = DBLE(LINT1)
            XPECT  = XPECT + XPECTL

C           Calculate shell quadruple (CD|AB).
C           ----------------------------------

            INT1(1:LINT1)=0.0D0
            CALL CHO_MCA_INT_1(ISHLCD,ISHLAB,INT1,LINT1,.FALSE.)

C           Calculate integrals from Cholesky vectors.
C           ------------------------------------------

            CALL CHO_DBGINT_CHO(INT1,NUMCD,NUMAB,WRK,
     &                          LWRK/2,ERRMAX,ERRMIN,ERRRMS,NCMP,
     &                          ISHLCD,ISHLAB)

C           Write report.
C           -------------

            IF (NCMP .LT. 1) THEN
               WRITE(LUPRI,'(4(I5,1X),5X,A)')
     &         ISHLC,ISHLD,ISHLA,ISHLB,' !!! nothing compared !!! '
            ELSE
               RMS = SQRT(ERRRMS/DBLE(NCMP))
               WRITE(LUPRI,'(4(I5,1X),1P,3(D12.4,1X))')
     &         ISHLC,ISHLD,ISHLA,ISHLB,ERRMIN,ERRMAX,RMS
            END IF

            IF (ABS(ERRMAX) .GT. ABS(GLMAX)) THEN
               GLMAX = ERRMAX
            END IF
            IF (ABS(ERRMIN) .LT. ABS(GLMIN)) THEN
               GLMIN = ERRMIN
            END IF
            GLRMS = GLRMS + ERRRMS
            IF (NCMP .GT. 0) XTCMP = XTCMP + DBLE(NCMP)

         END DO

      END DO

C     Print end of table.
C     -------------------

      WRITE(LUPRI,'(A)')
     & '--------------------------------------------------------------'
      IF (XTCMP .LT. 1.0D0) THEN
         WRITE(LUPRI,'(A,23X,A)')
     &   'Total:',' !!! nothing compared !!! '
      ELSE
         GLRMS = SQRT(GLRMS/XTCMP)
         WRITE(LUPRI,'(A,18X,1P,3(D12.4,1X))')
     &   'Total:',GLMIN,GLMAX,GLRMS
      END IF
      WRITE(LUPRI,'(A)')
     & '--------------------------------------------------------------'

C     Release all memory allocated here (and release seward memory).
C     --------------------------------------------------------------

      CALL XRLSMEM_INTS()
      Call mma_deallocate(WRK)
      Call mma_deallocate(INT1)

C     Check total number of comparisons.
C     ----------------------------------

      DO ISYM = 1,NSYM
         XLBAS(ISYM) = DBLE(NBAS(ISYM))
      END DO
      XNINT = 0.0D0
      DO ISYM = 1,NSYM
         XXLBST = 0.0D0
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYM)
            IF (ISYMA .EQ. ISYMB) THEN
               XXLBST = XXLBST + XLBAS(ISYMA)*(XLBAS(ISYMA)+1.0D0)/2.0D0
            ELSE IF (ISYMA .GT. ISYMB) THEN
               XXLBST = XXLBST + XLBAS(ISYMA)*XLBAS(ISYMB)
            END IF
         END DO
         XNINT = XNINT + XXLBST*(XXLBST + 1.0D0)/2.0D0
      END DO

      IF (ABS(XTCMP-XPECT) .GT. 1.0D-15) THEN
         WRITE(LUPRI,'(/,A)')
     &   'WARNING: not all integrals checked:'
      ELSE
         WRITE(LUPRI,*)
      END IF
      WRITE(LUPRI,'(A,1P,D20.10)')
     & 'Total number of integral comparisons    :',XTCMP
      WRITE(LUPRI,'(A,1P,D20.10)')
     & 'Total number expected (full shell pairs):',XPECT
      WRITE(LUPRI,'(A,1P,D20.10)')
     & 'Total number of unique integrals        :',XNINT

      END
