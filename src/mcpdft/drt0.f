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
      SUBROUTINE DRT0_m(A0,B0,C0,NVERT,DRT,DOWN,NTMP,TMP)
C
C     PURPOSE: CONSTRUCT THE UNRESTRICTED GUGA TABLE
C
      IMPLICIT INTEGER (A-Z)
C
      DIMENSION DRT(NVERT,5),DOWN(NVERT,0:3),TMP(NTMP)
      PARAMETER(LTAB=1,NTAB=2,ATAB=3,BTAB=4,CTAB=5)
      DIMENSION DA(0:3),DB(0:3),DC(0:3)
      DATA DA /0,0,1,1/, DB /0,1,-1,0/, DC /1,0,1,0/
C
C     SET UP TOP ROW
C
      NACTEL=2*A0+B0
      NLEV=A0+B0+C0
      DRT(1,LTAB)=NLEV
      DRT(1,NTAB)=NACTEL
      DRT(1,ATAB)=A0
      DRT(1,BTAB)=B0
      DRT(1,CTAB)=C0
      VSTA=1
      VEND=1
C
C     LOOP OVER ALL LEVELS
C
      DO LEV=NLEV,1,-1
        MXADDR=((LEV+1)*(LEV+2))/2
        DO I=1,MXADDR
          TMP(I)=0
        END DO
C
C     LOOP OVER VERTICES
C
        DO VUP=VSTA,VEND
          AUP=DRT(VUP,ATAB)
          BUP=DRT(VUP,BTAB)
          CUP=DRT(VUP,CTAB)
C
C     LOOP OVER CASES
C     AND STORE ONLY VALID CASE NUMBERS WITH ADRESSES
C
          DO STEP=0,3
            DOWN(VUP,STEP)=0
            ADWN=AUP-DA(STEP)
            IF(ADWN.LT.0) GOTO 20
            BDWN=BUP-DB(STEP)
            IF(BDWN.LT.0) GOTO 20
            CDWN=CUP-DC(STEP)
            IF(CDWN.LT.0) GOTO 20
            BC=BDWN+CDWN
            ADDR=1+(BC*(BC+1))/2 + CDWN
            TMP(ADDR)=4*VUP+STEP
            DOWN(VUP,STEP)=ADDR
  20        CONTINUE
          END DO
        END DO
        VDWN=VEND
C
C     NOW INSERT VALID CASES INTO DRT TABLE
C
        DO ADDR=1,MXADDR
          VUPS=TMP(ADDR)
          IF(VUPS.EQ.0) GOTO 40
          VUP=VUPS/4
          STEP=MOD(VUPS,4)
          VDWN=VDWN+1
          DRT(VDWN,ATAB)=DRT(VUP,ATAB)-DA(STEP)
          DRT(VDWN,BTAB)=DRT(VUP,BTAB)-DB(STEP)
          DRT(VDWN,CTAB)=DRT(VUP,CTAB)-DC(STEP)
          TMP(ADDR)=VDWN
  40      CONTINUE
        END DO
C
C     CREATE DOWN CHAIN TABLE
C
        DO VUP=VSTA,VEND
          DO STEP=0,3
            DWN=DOWN(VUP,STEP)
            IF(DWN.NE.0) DOWN(VUP,STEP)=TMP(DWN)
          END DO
        END DO
        VSTA=VEND+1
        VEND=VDWN
      END DO
* End of loop over levels.
C
C     ADDING THE ZERO LEVEL TO DRT AND DOWNCHAIN TABLE
C
      DO  I=1,5
        DRT(VEND,I)=0
      END DO
      DO STEP=0,3
        DOWN(VEND,STEP)=0
      END DO
C
C     COMPLETE DRT TABLE BY ADDING NO. OF ORBITALS AND ELECTRONS
C     INTO THE FIRST AND SECOND COLUMN
C
      DO VERT=1,VEND
        DRT(VERT,LTAB)=DRT(VERT,ATAB)+DRT(VERT,BTAB)+DRT(VERT,CTAB)
        DRT(VERT,NTAB)=2*DRT(VERT,ATAB)+DRT(VERT,BTAB)
      END DO
C
      RETURN
      END
