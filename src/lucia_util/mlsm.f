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
      SUBROUTINE MLSM(IML,IPARI,ISM,TYPE,IWAY)
*
* Transfer between ML,IPARI notation and compound notation ISM
*
* IWAY = 1 : IML,IPARI => Compound
* IWAY = 2 : IML,IPARI <= Compound
*
* TYPE : 'SX','OB','ST','DX','CI'
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.fh"
      CHARACTER*2 TYPE
*./NONAB/
* FIXME: This common block is never initialized or referenced elsewhere
* In the absence of further information, I'm setting the values to 0
*      COMMON /NONAB/ MNMLOB,NMLOB,
*     &               MNMLST,NMLST,
*     &               MNMLSX,NMLSX,
*     &               MNMLCI,NMLCI,
*     &               MNMLDX,NMLDX
       PARAMETER (MNMLOB=0,NMLOB=0,
     &            MNMLST=0,NMLST=0,
     &            MNMLSX=0,NMLSX=0,
     &            MNMLCI=0,NMLCI=0,
     &            MNMLDX=0,NMLDX=0)
*./CSM/
C     COMMON/CSM/NSMSX,NSMDX,NSMST,NSMCI,ITSSX,ITSDX
#include "csm.fh"
*
*.(Tired of warnings from 3090 Compiler so )
* (
      NML = 0
      MNML= 0
*             )
      IF(TYPE.EQ.'OB') THEN
        NML = NMLOB
        MNML = MNMLOB
      ELSE IF(TYPE.EQ.'SX') THEN
        NML = NMLSX
        MNML = MNMLSX
      ELSE IF(TYPE.EQ.'DX') THEN
        NML = NMLDX
        MNML = MNMLDX
      ELSE IF(TYPE.EQ.'ST') THEN
        NML = NMLST
        MNML = MNMLST
      ELSE IF(TYPE.EQ.'CI') THEN
        NML = NMLCI
        MNML = MNMLCI
      END IF
*
      IF(IWAY.EQ.1) THEN
C        ISM = (IPARI-1)*NML + MNML - 1
         ISM = (IPARI-1)*NML + IML - MNML + 1
      ELSE IF(IWAY.EQ.2) THEN
        IF(ISM.GT.NML) THEN
          IPARI = 2
          IML = ISM - NML + MNML - 1
        ELSE
          IPARI = 1
          IML = ISM       + MNML - 1
        END IF
      ELSE
        WRITE(6,*) ' Error in MLSM , IWAY = ' ,IWAY
        WRITE(6,*) ' MLSM stop !!! '
*        STOP 20
        CALL SYSABENDMSG('lucia_util/mlsm','Internal error',' ')
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,A)') ' MLSM speaking ,type= ',TYPE
        WRITE(6,'(A,3I4)') ' IML IPARI ISM ',IML,IPARI,ISM
      END IF
*
      RETURN
      END
