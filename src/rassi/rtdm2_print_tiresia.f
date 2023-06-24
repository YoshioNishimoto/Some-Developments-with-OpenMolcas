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
      INTEGER IZYI,IZYJ,IZYL,I,IZ,IOZ,JOZ,LOZ
      INTEGER NOI,NAI,NII,NOJ,NAJ,NIJ,NOL,NAL,NIL,KPOS
      INTEGER x2sym(8), xorb(NOSHT)
      INTEGER IA,IO,JA,JO,LA,LO
      Real*8  RT2M(NRT2MAB)
C ------------------------------------------------------------
C Other variables
      DIMENSION IOFFA(8), IOFFO(8)

      ! x2sym() translates from Molcas notation to Tiresia
      IF (NSYM.LT.2) THEN
        data x2sym  /  1, 3, 2, 4, 5, 7, 6, 8 /
      ELSE
        data x2sym  /  1, 2, 3, 4, 5, 6, 7, 8 /
      END IF

C IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFA=0
      DO I=1,NSYM-1
        IOFFA(I+1)=IOFFA(I)+NASH(I)
      END DO
C IOFFO=NR OF OCC ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFO=0
      DO I=1,NSYM-1
        IOFFO(I+1)=IOFFO(I)+NOSH(I)
      END DO
C Subroutine starts
      LU=54
      LU=IsFreeUnit(LU)
      !WRITE(NUM1,'(I3.3)') ISTATE
      !WRITE(NUM2,'(I3.3)') JSTATE
      !CALL Molcas_Open(LU,'TIRESIA_r2TM')
      OPEN(LU,file='TIRESIA_r2TM',action='write',position='append')
C xorb() translates from Molcas notation to Tiresia
      K=1
      DO I=1,NSYM
         IZ=x2sym(I)
       DO J=1,NOSH(IZ)
        xorb(k) = J+IOFFO(IZ)
        K=K+1
       END DO
      END DO
      SYM12=x2sym( MUL( LSYM1,LSYM2 ))
      WRITE(LU,'(A27,8I7)')'# Total Symm of WF product:', SYM12
      WRITE(LU,'(A9,8I7,8I7)')'# States:',ISTATE, JSTATE
C Write reduced 2-e TDM in CI basis.
      IOFFTD=0
      WRITE(LU,'(A43)')'# 2-e reduced TDM for CI coeff. in MO basis'
      WRITE(LU,'(A26,8I7)')'# Symmetry Block elements:',NRT2MAB
      !WRITE(LU,*) NRT2MAB
      WRITE(LU,'(A45)')'# sub-Block info:Sym(I,J,L), NumOrb SymmBlock'
      DO ISYI=1,NSYM
       IZYI=x2sym(ISYI)
       NOI=NOSH( IZYI )
       NAI=NASH( IZYI )
       NII=NISH( IZYI )
       IF(NOI.EQ.0) GOTO 270
       DO ISYJ=1,NSYM
        IZYJ=x2sym(ISYJ)
        NOJ=NOSH( IZYJ )
        NAJ=NASH( IZYJ )
        NIJ=NISH( IZYJ )
        IF(NOJ.EQ.0) GOTO 370
        DO ISYL=1,NSYM
         IZYL=x2sym(ISYL)
         NOL=NOSH( IZYL )
         NAL=NASH( IZYL )
         NIL=NISH( IZYL )
         IF(NOL.EQ.0) GOTO 470
          IF(MUL(ISYI, MUL(ISYJ,ISYL)).EQ.SYM12 ) THEN
            IF(NAI.EQ.0) GOTO 670
            IF(NAJ.EQ.0) GOTO 670
            IF(NAL.EQ.0) GOTO 670
            WRITE(LU,'(A12,8I7,8I7,8I7,8I7)')'# sub-Block:',ISYI,ISYJ,
     &   ISYL,NOI*NOJ*NOL
            DO I=1,NOI
             IA=IOFFA(IZYI)+I-NII
             IO=IOFFO(IZYI)+I
             DO J=1,NOJ
              JA=IOFFA(IZYJ)+J-NIJ
              JO=IOFFO(IZYJ)+J
              DO L=1,NOL
               LA=IOFFA(IZYL)+L-NIL
               LO=IOFFO(IZYL)+L
               IOZ=xorb(IO) ! orbital in flipped order
               JOZ=xorb(JO)
               LOZ=xorb(LO)
               IF((IA.LE.0).or.(JA.LE.0).or.(LA.LE.0)) THEN
               write(LU,'(I7,I7,I7,E26.12)') IOZ,JOZ,LOZ,0.0D0
               ELSE
                KPOS=IA+NASHT*((LA+NASHT*(JA-1))-1)
                IF(ABS(RT2M(KPOS)).LT.1.0D-15) THEN
                 RT2M(KPOS) = 0.0D0
                END IF
                write(LU,'(I7,I7,I7,E26.12)') IOZ,JOZ,LOZ,
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
