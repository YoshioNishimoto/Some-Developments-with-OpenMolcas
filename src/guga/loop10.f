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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
************************************************************************
      SUBROUTINE LOOP10(KM,ISTOP,IT1,IT2)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      CALL QENTER('LOOP10')
      ISTOP=0
      KM1=KM+1
      J1F=IPO(J1(KM1))
      IDIF=IA(J1(KM1))-IA(J2(KM1))
      IF(IDIF.LT.0.OR.IDIF.GT.1)GO TO 55
      IF(IDIF.EQ.0)GO TO 45
      IWAYKM=IWAY(KM)
      GO TO (39,41,42,43,55),IWAYKM
39    IWAY(KM)=2
C     (CJ,AQ)
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 41
      IF(K1F(J1F).EQ.0)GO TO 141
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP1(KM)=BL1(IB(J2(KM1))+1)*COUP(KM1)
      JM1(KM)=K1F(J1F)
      IF(K2F(J1F).EQ.0)GO TO 40
      GO TO 241
141   IF(K2F(J1F).EQ.0)GO TO 41
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
241   COUP(KM)=COUP(KM1)/IB(J2(KM1))
      JM(KM)=K2F(J1F)
      GO TO 40
C     AK
41    IWAY(KM)=3
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 42
      IF(K2F(J1F).EQ.0)GO TO 42
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP1(KM)=-COUP(KM1)
      JM1(KM)=K2F(J1F)
      GO TO 40
C     BL
42    IWAY(KM)=4
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 43
      IF(K3F(J1F).EQ.0)GO TO 43
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP1(KM)=-BS1(IB(J2(KM1)))*COUP(KM1)
      JM1(KM)=K3F(J1F)
      GO TO 40
C     DL
43    IWAY(KM)=5
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(J1F).EQ.0)GO TO 55
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP1(KM)=-BS2(IB(J2(KM1)))*COUP(KM1)
      JM1(KM)=K3F(J1F)
      GO TO 40
45    IWAYKM=IWAY(KM)
      GO TO (59,61,62,63,55),IWAYKM
59    IWAY(KM)=2
C     (CR,AO)
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K2(IT2+J2(KM1)).EQ.0)GO TO 61
      IF(K1F(J1F).EQ.0)GO TO 161
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
      COUP1(KM)=-COUP(KM1)/(IB(J2(KM1))+2)
      JM1(KM)=K1F(J1F)
      IF(K2F(J1F).EQ.0)GO TO 40
      GO TO 261
161   IF(K2F(J1F).EQ.0)GO TO 61
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K2(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),2)
261   COUP(KM)=BL2(IB(J2(KM1))+1)*COUP(KM1)
      JM(KM)=K2F(J1F)
      GO TO 40
C     BP
61    IWAY(KM)=3
      IF(K1(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 62
      IF(K3F(J1F).EQ.0)GO TO 62
      J1(KM)=K1(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM)=-BS1(IB(J2(KM1))+2)*COUP(KM1)
      JM(KM)=K3F(J1F)
      GO TO 40
C     CN
62    IWAY(KM)=4
      IF(K0(IT1+J1(KM1)).EQ.0.OR.K1(IT2+J2(KM1)).EQ.0)GO TO 63
      IF(K1F(J1F).EQ.0)GO TO 63
      J1(KM)=K0(IT1+J1(KM1))
      J2(KM)=K1(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),1)
      COUP(KM)=-COUP(KM1)
      JM(KM)=K1F(J1F)
      GO TO 40
C     DP
63    IWAY(KM)=5
      IF(K2(IT1+J1(KM1)).EQ.0.OR.K3(IT2+J2(KM1)).EQ.0)GO TO 55
      IF(K3F(J1F).EQ.0)GO TO 55
      J1(KM)=K2(IT1+J1(KM1))
      J2(KM)=K3(IT2+J2(KM1))
      ICOUP1(KM)=ICOUP1(KM1)+IY(IT1+J1(KM1),2)
      ICOUP(KM)=ICOUP(KM1)+IY(IT2+J2(KM1),3)
      COUP(KM)=-BS2(IB(J2(KM1))+2)*COUP(KM1)
      JM(KM)=K3F(J1F)
      GO TO 40
55    ISTOP=1
40      Continue
       CALL QEXIT('LOOP10')
      RETURN
      END
