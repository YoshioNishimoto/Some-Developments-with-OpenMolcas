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
      FUNCTION IOFF_SYM_DIST(ISYM,NGASL,IOFF,MAXVAL,MINVAL)
*
* A ts block of string is given and the individual
* symmetrydisrtributions has been obtained ( for example
* by TS_SYM_PNT)
*
* Obtain offset for symmetrycombination defined by ISYM
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER ISYM(*),IOFF(*),MAXVAL(*),MINVAL(*)
* Address in IOFF is
*     1
*     +  (ISM1-MINVAL(1))
*     +  (ISM2-MINVAL(2))*(MAXVAL(1)-MINVAL(1)+1)
*     +  (ISM3-MINVAL(3))*(MAXVAL(1)-MINVAL(1)+1)*(MAXVAL(2)-MINVAL(2)+1)
*     +
*     +
*     +
*     +  (ISM L-1-MINVAL(L-1))*Prod(i=1,L-2)(MAXVAL(i)-MINVAL(i)+1)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       write(6,*) ' Isym, minval, ioff:'
       call iwrtma(isym,1,ngasl,1,ngasl)
       call iwrtma(minval,1,ngasl,1,ngasl)
       call iwrtma(ioff,1,ngasl,1,ngasl)
      END IF
*. Offset for this symmetry distribution in IOFFI
      I = 1
      IMULT = 1
      DO IGAS = 1, NGASL-1
        I = I + (ISYM(IGAS)-MINVAL(IGAS)) * IMULT
        IMULT = IMULT*(MAXVAL(IGAS)-MINVAL(IGAS)+1)
c        write(6,*) ' igas,i,imult ',igas,i,imult
      END DO
c The following IF block is needed for avoinging going outside the IOFF bounds.
c This is possible for certain GAS setups. Test 897 helped in finding this issue.
      IF (I.le.0) THEN
        IOFF_SYM_DIST=0
      ELSE
        IOFF_SYM_DIST=IOFF(I)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Info from IOFF_SYM_DIST'
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Address and offset ',I,IOFF_SYM_DIST
        WRITE(6,*) ' Symmetry distribution : ', (ISYM(J),J=1,NGASL)
      END IF
*
      RETURN
      END
