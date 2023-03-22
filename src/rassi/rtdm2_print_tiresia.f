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
* Copyright (C) 2023, Bruno Tenorio                                    *
************************************************************************

C Print the reduced 2-e TDM in ASCII format.
      SUBROUTINE RTDM2_PRINT_TIRESIA(ISTATE, JSTATE, NRT2MAB, RT2M)

      IMPLICIT REAL*8 (A-H,O-Z)

#include "prgm.fh"
      CHARACTER*12 ROUTINE
      PARAMETER (ROUTINE='RTDM2_PRINT_TIRESIA')
#include "rasdim.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "Files.fh"
#include "Struct.fh"
#include "stdalloc.fh"
      INTEGER NRT2MAB
      INTEGER ISTATE, JSTATE, SYM12,ISYI,ISYJ,ISYL
      INTEGER NOI,NAI,NII,NOJ,NAJ,NIJ,NOL,NAL,NIL,KPOS
      INTEGER IA,IO,JA,JO,LA,LO
      Real*8  RT2M(NRT2MAB)
C ------------------------------------------------------------
C Other variables
      DIMENSION IOFFA(8), IOFFO(8)
C IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFA(1)=0
      DO I=1,NSYM-1
        IOFFA(I+1)=IOFFA(I)+NASH(I)
      END DO
C IOFFO=NR OF OCC ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFO(1)=0
      DO J=1,NSYM-1
        IOFFO(J+1)=IOFFO(J)+NOSH(J)
      END DO
C Subroutine starts
      LU=54
      LU=IsFreeUnit(LU)
      !WRITE(NUM1,'(I3.3)') ISTATE
      !WRITE(NUM2,'(I3.3)') JSTATE
      !CALL Molcas_Open(LU,'TIRESIA_r2TM')
      OPEN(LU,file='TIRESIA_r2TM',action='write',position='append')

      SYM12=MUL(LSYM1,LSYM2)
      WRITE(LU,'(A31,8I7)')'# Total Symm of the WF product:',SYM12
      WRITE(LU,'(A9,8I7,8I7)')'# States:',ISTATE, JSTATE
C Write reduced 2-e TDM in CI basis.
      IOFFTD=0
      WRITE(LU,'(A43)')'# 2-e reduced TDM for CI coeff. in MO basis'
      WRITE(LU,'(A26,8I7)')'# Symmetry Block elements:',NRT2MAB
      !WRITE(LU,*) NRT2MAB
      WRITE(LU,'(A45)')'# sub-Block info:Sym(I,J,L), NumOrb SymmBlock'
      DO ISYI=1,NSYM
       NOI=NOSH(ISYI)
       NAI=NASH(ISYI)
       NII=NISH(ISYI)
       IF(NOI.EQ.0) GOTO 270
       DO ISYJ=1,NSYM
        NOJ=NOSH(ISYJ)
        NAJ=NASH(ISYJ)
        NIJ=NISH(ISYJ)
        IF(NOJ.EQ.0) GOTO 370
        DO ISYL=1,NSYM
         NOL=NOSH(ISYL)
         NAL=NASH(ISYL)
         NIL=NISH(ISYL)
         IF(NOL.EQ.0) GOTO 470
          IF(MUL(ISYI,MUL(ISYJ,ISYL)).EQ.SYM12) THEN
            IF(NAI.EQ.0) GOTO 670
            IF(NAJ.EQ.0) GOTO 670
            IF(NAL.EQ.0) GOTO 670
            WRITE(LU,'(A12,8I7,8I7,8I7,8I7)')'# sub-Block:',ISYI,ISYJ,
     &       ISYL,NOI*NOJ*NOL
            DO I=1,NOI
             IA=IOFFA(ISYI)+I-NII
             IO=IOFFO(ISYI)+I
             DO J=1,NOJ
              JA=IOFFA(ISYJ)+J-NIJ
              JO=IOFFO(ISYJ)+J
              DO L=1,NOL
               LA=IOFFA(ISYL)+L-NIL
               LO=IOFFO(ISYL)+L
               IF((IA.LE.0).or.(JA.LE.0).or.(LA.LE.0)) THEN
               write(LU,'(I7,I7,I7,E26.12)') IO,JO,LO,0.0D0
               ELSE
                KPOS=IA+NASHT*((LA+NASHT*(JA-1))-1)
                IF(ABS(RT2M(KPOS)).LT.1.0D-19) THEN
                 RT2M(KPOS) = 0.0D0
                END IF
                write(LU,'(I7,I7,I7,E26.12)') IO,JO,LO,
     &           RT2M(KPOS)
               END IF
              END DO
             END DO
            END DO
670         CONTINUE
          END IF
470     CONTINUE
        END DO
370    CONTINUE
       END DO
270   CONTINUE
      END DO
      CLOSE (LU)
      END SUBROUTINE
