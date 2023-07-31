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
      SUBROUTINE CHO_MCA_INIT(SKIP_PRESCREEN)
C
C     Purpose: initialization of Cholesky decomposition in MOLCAS.
C
      use index_arrays, only: iSO2Sh
      use ChoArr, only: iSOShl, iBasSh, nBasSh, nBstSh, iSP2F, iShlSO,
     &                  iShP2RS, iShP2Q
#include "implicit.fh"
      LOGICAL SKIP_PRESCREEN
#include "cholesky.fh"
#include "choorb.fh"
#include "stdalloc.fh"

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_INIT')

c     INTEGER  CHO_ISAO
c     EXTERNAL CHO_ISAO
c     INTEGER  CHO_ISAOSH
c     EXTERNAL CHO_ISAOSH

C     Check that the number of shells is within limits.
C     -------------------------------------------------

      IF (NSHELL .LT. 1) THEN
         WRITE(LUPRI,*) 'NSHELL out of bounds: ',NSHELL
         CALL CHO_QUIT('NSHELL out of bounds in '//SECNAM,102)
      END IF

C     Compute total #shell pair.
C     --------------------------

      NNSHL_TOT = NSHELL*(NSHELL + 1)/2
      IF (NNSHL_TOT .LT. 1) THEN
         WRITE(LUPRI,*)
     &   'NNSHL_TOT=NSHELL*(NSHELL+1)/2 is non-positive: ',
     &   NNSHL_TOT
         WRITE(LUPRI,*) 'Integer overflow ?'
         CALL CHO_QUIT('NNSHL_TOT out of bounds in '//SECNAM,102)
      END IF

C     Compute contributing #shell pair by prescreening (if requested).
C     iSP2F(k): returns global shell pair of contributing shell pair k.
C               (Allocated and defined in CHO_DIAPS.)
C     -----------------------------------------------------------------

      IF (SKIP_PRESCREEN) THEN
         IF (NNSHL.LT.1 .OR. NNSHL.GT.NNSHL_TOT) THEN
            WRITE(LUPRI,*) SECNAM,': flag SKIP_PRESCREEN is ',
     &                     SKIP_PRESCREEN
            WRITE(LUPRI,*) 'NNSHL is out-of-bounds: ',NNSHL
            WRITE(LUPRI,*) 'Condition: 0 < NNSHL < ',NNSHL_TOT
            CALL CHO_QUIT('Initialization error in '//SECNAM,102)
         END IF
         IF (SIZE(iSP2F) .NE. NNSHL) THEN
            WRITE(LUPRI,*) SECNAM,': flag SKIP_PRESCREEN is ',
     &                     SKIP_PRESCREEN
            WRITE(LUPRI,*) 'NNSHL is: ',NNSHL
            WRITE(LUPRI,*) 'SIZE(iSP2F) must be equal to NNSHL, ',
     &                     'SIZE(iSP2F)= ',SIZE(iSP2F)
            CALL CHO_QUIT('Initialization error in '//SECNAM,102)
         END IF
      ELSE
         CALL CHO_DIASP()
      END IF

C     Get the number of symmetries.
C     -----------------------------

      CALL GET_ISCALAR('nSym',NSYM)  ! Get # irreps.
      IF ((NSYM.LT.1) .OR. (NSYM.GT.8)) THEN
         WRITE(LUPRI,*) 'NSYM out of bounds: ',NSYM
         CALL CHO_QUIT('NSYM out of bounds in '//SECNAM,102)
      END IF

C     NBAS(ISYM): # basis functions (SOs) in symmetry ISYM
C     IBAS(ISYM): offset to basis functions in symmetry ISYM
C     NBAST     : total number of basis functions
C     ------------------------------------------------------

      CALL GET_IARRAY('nBas',NBAS,NSYM)
      IBAS(1) = 0
      NBAST   = NBAS(1)
      DO ISYM = 2,NSYM
         IBAS(ISYM) = NBAST
         NBAST = NBAST + NBAS(ISYM)
      END DO
      IF (NBAST .LT. 1) THEN
         WRITE(LUPRI,*) 'NBAST out of bounds: ',NBAST
         CALL CHO_QUIT('NBAST out of bounds in '//SECNAM,102)
      END IF

C     Allocate shell based index arrays.
C     ----------------------------------

      Call mma_allocate(iBasSh,nSym,nShell,Label='iBasSh')
      Call mma_allocate(nBasSh,nSym,nShell,Label='nBasSh')
      Call mma_allocate(nBstSh,nShell,Label='nBstSh')

C     ISOSHL(I): shell to which SO I belongs
C     --------------------------------------

      Call mma_allocate(iSOShl,NBAST,Label='iSOShl')
      DO ISYM = 1,NSYM
         DO IA = 1,NBAS(ISYM)
            I = IBAS(ISYM) + IA
            iSOShl(I) = ISO2SH(I)
         END DO
      END DO

C     NBASSH(ISYM,ISHL): total dimension of shell ISHL, sym. ISYM
C     NBSTSH(ISHL): total dimension of shell ISHL
C     MXORSH      : max. shell dimension
C     -----------------------------------------------------------

      CALL CHO_SETSH(IBASSH,NBASSH,NBSTSH,
     &               IBAS,NBAS,ISOSHL,NSYM,NSHELL,NBAST)

      MXORSH = NBSTSH(1)
      DO ISHL = 2,NSHELL
         MXORSH = MAX(MXORSH,NBSTSH(ISHL))
      END DO

C     MX2SH: max. dimension of contributing shell pair.
C     -------------------------------------------------

      MX2SH = -1
      DO IJSHL = 1,NNSHL
         IJ = iSP2F(IJSHL)
         CALL CHO_INVPCK(IJ,I,J,.TRUE.)
         IF (I .EQ. J) THEN
            NUMIJ = NBSTSH(I)*(NBSTSH(I) + 1)/2
         ELSE
            NUMIJ = NBSTSH(I)*NBSTSH(J)
         END IF
         MX2SH = MAX(MX2SH,NUMIJ)
      END DO
      IF (MX2SH .LT. 1) THEN
         WRITE(LUPRI,*) 'Max. shell pair dimension non-positive: ',
     &                  MX2SH
         CALL CHO_QUIT('Initialization problem in '//SECNAM,102)
      END IF

C     If needed, allocate memory for extracting qualified columns
C     directly in reduced set from Seward.
C     -----------------------------------------------------------

      IF (IFCSEW .EQ. 2) THEN
         Call mma_allocate(iShP2RS,2,Mx2Sh,Label='iShP2RS')
         Call mma_allocate(iShP2Q ,2,Mx2Sh,Label='iShP2Q ')
      END IF

C     ISHLSO(I): index of SO I within its shell
C     -----------------------------------------

      Call mma_allocate(iShlSO,nBasT,Label='iShlSO')
      CALL CHO_SETSH2(iShlSO,iSOShl,NBSTSH,NBAST,NSHELL)

      END
