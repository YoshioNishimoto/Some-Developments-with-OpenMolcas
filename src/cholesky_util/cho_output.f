!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1971, Nelson H. F. Beebe                               *
!***********************************************************************
      SUBROUTINE CHO_OUTPUT(AMATRX,ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,
     *                      COLDIM,NCTL,LUPRI)
!.......................................................................
!
! OUTPUT PRINTS A REAL MATRIX IN FORMATTED FORM WITH NUMBERED ROWS
! AND COLUMNS.  THE INPUT IS AS FOLLOWS;
!
!        AMATRX(',').........MATRIX TO BE OUTPUT
!
!        ROWLOW..............ROW NUMBER AT WHICH OUTPUT IS TO BEGIN
!
!        ROWHI...............ROW NUMBER AT WHICH OUTPUT IS TO END
!
!        COLLOW..............COLUMN NUMBER AT WHICH OUTPUT IS TO BEGIN
!
!        COLHI...............COLUMN NUMBER AT WHICH OUTPUT IS TO END
!
!        ROWDIM..............ROW DIMENSION OF AMATRX(',')
!
!        COLDIM..............COLUMN DIMENSION OF AMATRX(',')
!
!        NCTL................CARRIAGE CONTROL FLAG; 1 FOR SINGLE SPACE
!                                                   2 FOR DOUBLE SPACE
!                                                   3 FOR TRIPLE SPACE
!
! THE PARAMETERS THAT FOLLOW MATRIX ARE ALL OF TYPE INTEGER*4.  THE
! PROGRAM IS SET UP TO HANDLE 5 COLUMNS/PAGE WITH A 1P,5D24.15 FORMAT
! FOR THE COLUMNS.  IF A DIFFERENT NUMBER OF COLUMNS IS REQUIRED,
! CHANGE FORMATS 1000 AND 2000, AND INITIALIZE KCOL WITH THE NEW NUMBER
! OF COLUMNS.
!
! AUTHOR;  NELSON H.F. BEEBE, QUANTUM THEORY PROJECT, UNIVERSITY OF
!          FLORIDA, GAINESVILLE
! REVISED; FEBRUARY 26, 1971
!
!.......................................................................
!
      Implicit Real*8 (a-h,o-z)
      INTEGER   ROWLOW,ROWHI,COLLOW,COLHI,ROWDIM,COLDIM,BEGIN,KCOL
      DIMENSION AMATRX(ROWDIM,COLDIM)
      CHARACTER*1 ASA(3), BLANK, CTL
      CHARACTER   PFMT*20, COLUMN*8
      PARAMETER (ZERO=0.D00, KCOLP=4, KCOLN=6)
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA COLUMN/'Column  '/, BLANK/' '/, ASA/' ', '0', '-'/
!
      IF (ROWHI.LT.ROWLOW) GO TO 3
      IF (COLHI.LT.COLLOW) GO TO 3
!
      AMAX = ZERO
      DO J = COLLOW,COLHI
         DO I = ROWLOW,ROWHI
            AMAX = MAX( AMAX, ABS(AMATRX(I,J)) )
         END DO
      END DO
      IF (AMAX .EQ. ZERO) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      IF (FFMIN .LE. AMAX .AND. AMAX .LE. FFMAX) THEN
!        use F output format
         PFMT = '(A1,I7,2X,8F15.8)'
      ELSE
!        use 1PD output format
         PFMT = '(A1,I7,2X,1P,8D15.6)'
      END IF
!
      IF (NCTL .LT. 0) THEN
         KCOL = KCOLN
      ELSE
         KCOL = KCOLP
      END IF
      MCTL = ABS(NCTL)
      IF ((MCTL.LE.3).AND.(MCTL.GT.0)) THEN
         CTL = ASA(MCTL)
      ELSE
         CTL = BLANK
      END IF
!
      LAST = MIN(COLHI,COLLOW+KCOL-1)
      DO 2 BEGIN = COLLOW,COLHI,KCOL
         WRITE (LUPRI,1000) (COLUMN,I,I = BEGIN,LAST)
         DO 1 K = ROWLOW,ROWHI
            DO 4 I = BEGIN,LAST
               IF (AMATRX(K,I).NE.ZERO) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,PFMT) CTL,K,(AMATRX(K,I), I = BEGIN,LAST)
    1    CONTINUE
         LAST = MIN(LAST+KCOL,COLHI)
    2 CONTINUE
    3 RETURN
 1000 FORMAT (/12X,6(3X,A6,I4,2X),(3X,A6,I4))
!2000 FORMAT (A1,'Row',I4,2X,1P,8D15.6)
!2000 FORMAT (A1,I7,2X,1P,8D15.6)
      END
