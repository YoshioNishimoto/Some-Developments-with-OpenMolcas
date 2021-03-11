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
      SUBROUTINE PRCMAT(NSS,XMATR,XMATI)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMATR(NSS,NSS),XMATI(NSS,NSS)
C Write out matrix elements over states as a complex matrix
C in square format
      DO JSTA=1,NSS,2
       JEND=MIN(NSS,JSTA+1)
       WRITE(6,*)
       WRITE(6,'(1X,A8,12X,I3,35X,I3)')' STATE  ',(JSS,JSS=JSTA,JEND)
       DO ISS=1,NSS
       WRITE(6,'(1X,I4,2x,2(A1,F10.6,A1,F10.6,A1,3x))')
     &           ISS,('(',XMATR(ISS,JSS),',',XMATI(ISS,JSS),
     &           ')',JSS=JSTA,JEND)
       END DO
      END DO
      RETURN
      END

      SUBROUTINE PRCMAT2(INPUT,NSS,XMATR,XMATI)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMATR(NSS,NSS),XMATI(NSS,NSS)
#include "Molcas.fh"
#include "cntrl.fh"
      CHARACTER(LEN=8) PROPERTY
      CHARACTER(LEN=1) DIRECTION
      CHARACTER(LEN=200) FILENAME
C Write out matrix elements over states as a complex matrix
C in parsable format
      if(INPUT.gt.0) THEN
        PROPERTY = SOPRNM(INPUT)
      ELSE
        PROPERTY = 'EIGVEC'
      ENDIF

      WRITE(DIRECTION,'(I1)') ISOCMP(INPUT)
      IF (PROPERTY(1:4).EQ."MLTP") THEN
          IF (PROPERTY(8:8).EQ.'0') THEN
            FILENAME = 'monopole-'//DIRECTION//'.txt'
          ELSE IF (PROPERTY(8:8).EQ.'1') THEN
            FILENAME = 'dipole-'//DIRECTION//'.txt'
          ELSE IF (PROPERTY(8:8).EQ.'2') THEN
            FILENAME = 'quadrupole-'//DIRECTION//'.txt'
          ELSE
            GO TO 100
          END IF
      ELSE IF (PROPERTY(1:4).EQ."ANGM") THEN
          FILENAME = 'angmom-'//DIRECTION//'.txt'
      ELSE IF (PROPERTY(1:6).EQ."EIGVEC") THEN
          FILENAME = "eigvectors.txt"
      ELSE
          GO TO 100
      END IF
      OPEN(UNIT=88,FILE=FILENAME,STATUS='REPLACE')
      WRITE(88,*) "#NROW NCOL REAL IMAG"
      DO JSTA=1,NSS
        DO ISS=1,NSS
        WRITE(88,'(I4,I4,A1,E25.16,A1,E25.16)') ISS,JSTA,' ',
     &   XMATR(ISS,JSTA),' ',XMATI(ISS,JSTA)
        END DO
      END DO
      CLOSE(88)
 100  CONTINUE
      RETURN
      END

      SUBROUTINE PRCMAT3(NSS,SMATR,SMATI,DIR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SMATR(NSS,NSS), SMATI(NSS,NSS)
      INTEGER DIR
#include "Molcas.fh"
#include "cntrl.fh"
      CHARACTER(LEN=1) DIRECTION
      CHARACTER(LEN=200) FILENAME
C Write out spin matrix elements in parsable format
      WRITE(DIRECTION,'(I1)') DIR
      FILENAME = 'spin-'//DIRECTION//'.txt'
      OPEN(UNIT=88,FILE=FILENAME,STATUS='REPLACE')
      WRITE(88,*) "#NROW NCOL REAL IMAG"
      DO JSTA=1,NSS
        DO ISS=1,NSS
        WRITE(88,'(I4,I4,A1,E25.16,A1,E25.16)') ISS,JSTA,' ',
     &   SMATR(ISS,JSTA),' ',SMATI(ISS,JSTA)
        END DO
      END DO
      CLOSE(88)
      RETURN
      END

      SUBROUTINE MULMAT(NSS,XMATR,XMATI,ee,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 XMATR(NSS,NSS),XMATI(NSS,NSS)
      COMPLEX*16 Z(NSS,NSS)
      ee=0.d0
      DO ISS=1,NSS
      DO JSS=1,NSS
      Z(ISS,JSS)=(0.0d0,0.0d0)
      enddo
      enddo
      DO ISS=1,NSS
      DO JSS=1,NSS
      ee=ee+XMATR(ISS,JSS)*XMATR(ISS,JSS)+
     & XMATI(ISS,JSS)*XMATI(ISS,JSS)
      Z(ISS,JSS)=Z(ISS,JSS)+
     &DCMPLX(XMATR(ISS,JSS),XMATI(ISS,JSS))
      enddo
      enddo
      RETURN
      END

