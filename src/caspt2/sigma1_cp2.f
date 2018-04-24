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
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE SIGMA1_CP2(IP,IQ,CPQ,ISYCI,CI,SGM,NOCSF,IOCSF,NOW,IOW,
     &                 NOCP,IOCP,ICOUP,VTAB,MVL,MVR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NOCSF(NSYM,NMIDV,NSYM),IOCSF(NSYM,NMIDV,NSYM)
      DIMENSION NOW(2,NSYM,NMIDV),IOW(2,NSYM,NMIDV)
      DIMENSION NOCP(MXEO,NSYM,NMIDV),IOCP(MXEO,NSYM,NMIDV)
      DIMENSION VTAB(NVTAB)
      INTEGER ICOUP
      DIMENSION ICOUP(3,NICOUP)
      DIMENSION MVL(NMIDV,2),MVR(NMIDV,2)
      DIMENSION CI(*),SGM(*)


#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"
      INTRINSIC MOD
*****************************************************************
*  GIVEN ACTIVE LEVEL INDICES IP AND IQ, AND INPUT CI ARRAYS
*  CI AND SGM, THIS ROUTINE ADDS TO SGM THE RESULT OF ACTING ON
*  CI WITH THE NUMBER CPQ TIMES THE EXCITATION OPERATOR E(IP,IQ).
*  THE ADDITIONAL ENTRIES IN THE PARAMETER LIST ARE TABLES THAT
*  WERE PREPARED BY GINIT AND ITS SUBROUTINES.
*****************************************************************

      CALL QENTER('SIGMA1')

      IF(ABS(CPQ).LT.1.0D-12) GOTO 999
      iLUW=0
C SYMMETRY OF ORBITALS:
      ISYP=ISM(IP)
      ISYQ=ISM(IQ)
      ISYPQ=MUL(ISYP,ISYQ)
C SYMMETRY OF SIGMA ARRAY:
      ISYSGM=MUL(ISYPQ,ISYCI)

      IF(IP.GT.IQ) THEN
C THEN THIS IS AN EXCITING OPERATOR.
        IF(IP.LE.MIDLEV) GOTO 1600
        IF(IQ.GT.MIDLEV) GOTO 1700
        GOTO 1800
      ELSE IF(IQ.GT.IP) THEN
C THEN THIS IS A DEEXCITING OPERATOR.
        IF(IQ.LE.MIDLEV) GOTO 1300
        IF(IP.GT.MIDLEV) GOTO 1400
        GOTO 1500
      END IF
      IF(IP.GT.MIDLEV) GOTO 1200

C SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
C IP=IQ < MIDLEV.
      DO 100 MVSGM=1,NMIDV
        DO 100 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 100
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOLW=IOW(2,ISYDSG,MVSGM)
          IPSHFT=2*(IP-1)
          LLW=LICASE+IOLW-NIPWLK+IPSHFT/30
cVV next statement is a cheat to forbid optimization
          if(iLUW.eq.666) iLUW=iLUW+LLW
cVV
          IPSHFT=MOD(IPSHFT,30)
          IPPOW=2**IPSHFT
          DO 101 J=1,NDWNSG
            JC=IWORK(LLW+J*NIPWLK)
CPAM96           ICS=IAND(ISHFT(JC,-IPSHFT),3)
            ICS=MOD(JC/IPPOW,4)
            IF(ICS.EQ.0) GOTO 101
            X=CPQ*DBLE((1+ICS)/2)
            JSTA=ISGSTA+NUPSG*(J-1)
            CALL DAXPY_(NUPSG,X,CI(JSTA),1,SGM(JSTA),1)
 101      CONTINUE
 100  CONTINUE

      GOTO 999

 1200 CONTINUE
C SPECIAL CASE: WEIGHT OPERATOR, IP=IQ.
C IP=IQ>MIDLEV
      DO 200 MVSGM=1,NMIDV
        DO 200 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 200
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOUW=IOW(1,ISYUSG,MVSGM)
          IPSHFT=2*(IP-1-MIDLEV)
          LUW=LICASE+IOUW-NIPWLK+IPSHFT/30
cVV next statement is a cheat to forbid optimization
          if(iLUW.eq.666) iLUW=iLUW+LUW
cVV
          IPSHFT=MOD(IPSHFT,30)
          IPPOW=2**IPSHFT
          DO 102 I=1,NUPSG
            IC=IWORK(LUW+I*NIPWLK)
CPAM96            ICS=IAND(ISHFT(IC,-IPSHFT),3)
            ICS=MOD(IC/IPPOW,4)
            IF(ICS.EQ.0) GOTO 102
            X=CPQ*DBLE((1+ICS)/2)
            ISTA=ISGSTA-1+I
            CALL DAXPY_(NDWNSG,X,CI(ISTA),NUPSG,SGM(ISTA),NUPSG)
 102      CONTINUE
 200  CONTINUE

      GOTO 999

 1300 CONTINUE
C DEEXCITING OPERATOR, IP<IQ<=MIDLEV.
      DO 300 MVSGM=1,NMIDV
        DO 300 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 300
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYDC=MUL(ISYPQ,ISYDSG)
          NDWNC=NOW(2,ISYDC,MVSGM)
          IF(NDWNC.EQ.0) GOTO 300
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUSG,MVSGM,ISYCI)
          INDEO=2*NLEV+(IQ*(IQ-1))/2+IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            CALL DEX1_CP2 (CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
          END IF
 300  CONTINUE

      GOTO 999

 1400 CONTINUE
C DEEXCITING OPERATOR, MIDLEV<IP<IQ
      DO 400 MVSGM=1,NMIDV
        DO 400 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 400
          ISYUC=MUL(ISYPQ,ISYUSG)
          NUPC=NOW(1,ISYUC,MVSGM)
          IF (NUPC.EQ.0) GOTO 400
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOC=IOCSF(ISYUC,MVSGM,ISYCI)
          INDEO=2*NLEV+(IQ*(IQ-1))/2+IP
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
C CASE IS: UPPER HALF, DEEXCITE:
            CALL DEX2_CP2 (CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
          END IF
 400  CONTINUE

      GOTO 999

 1500 CONTINUE
C DEEXCITING CASE, IP<=MIDLEV<IQ.
C ALLOCATE TEMPORARY WORK AREA:
      DO 500 MVSGM=1,NMIDV
        MV4=MVR(MVSGM,1)
        MV5=MVR(MVSGM,2)
        DO 500 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 500
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYUC=MUL(ISYQ,ISYUSG)
          ISYDC=MUL(ISYP,ISYDSG)
          IF(MV4.EQ.0) GOTO 499
          NUPC=NOW(1,ISYUC,MV4)
          IF(NUPC.EQ.0) GOTO 499
          NDWNC=NOW(2,ISYDC,MV4)
          IF(NDWNC.EQ.0) GOTO 499
          INDEO=IQ
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.EQ.0) GOTO 499
          NTMP=NUPSG*NDWNC
CMO       vclr seems to be a burden on performace
CMO       replace it by newvclr which should be faster
C          CALL VCLR(WORK(LSGTMP),1,NTMP)
          CALL NEWVCLR(WORK(LSGTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUC,MV4,ISYCI)
C CASE IS: UPPER HALF, DEEXCITE:
          CALL DEX2_CP2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LSGTMP),
     &             NCP,ICOUP(1,LICP),VTAB)
          INDEO=IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            X=1.0D00
            CALL DEX1_CP2 (X,NUPSG,WORK(LSGTMP),SGM(ISGSTA),
     &               NCP,ICOUP(1,LICP),VTAB)
          END IF
 499      CONTINUE
          IF(MV5.EQ.0) GOTO 500
          NUPC=NOW(1,ISYUC,MV5)
          IF(NUPC.EQ.0) GOTO 500
          NDWNC=NOW(2,ISYDC,MV5)
          IF(NDWNC.EQ.0) GOTO 500
          INDEO=NLEV+IQ
          NCP=NOCP(INDEO,ISYUSG,MVSGM)
          IF(NCP.EQ.0) GOTO 500
          NTMP=NUPSG*NDWNC
C          CALL VCLR(WORK(LSGTMP),1,NTMP)
          CALL NEWVCLR(WORK(LSGTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUC,MV5,ISYCI)
C CASE IS: UPPER HALF, DEEXCITE:
          CALL DEX2_CP2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LSGTMP),
     &             NCP,ICOUP(1,LICP),VTAB)
          INDEO=NLEV+IP
          NCP=NOCP(INDEO,ISYDSG,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDSG,MVSGM)
C CASE IS: LOWER HALF, DEEXCITE:
            X=1.0D00
            CALL DEX1_CP2 (X,NUPSG,WORK(LSGTMP),SGM(ISGSTA),
     &               NCP,ICOUP(1,LICP),VTAB)
          END IF
 500  CONTINUE

      GOTO 999

 1600 CONTINUE
C EXCITING CASE, IQ<IP<=MIDLEV.
      DO 600 MVSGM=1,NMIDV
        DO 600 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 600
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYDC=MUL(ISYPQ,ISYDSG)
          NDWNC=NOW(2,ISYDC,MVSGM)
          IF(NDWNC.EQ.0) GOTO 600
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          IOC=IOCSF(ISYUSG,MVSGM,ISYCI)
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          NCP=NOCP(INDEO,ISYDC,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYDC,MVSGM)
C CASE IS: LOWER HALF, EXCITE:
            CALL EXC1_CP2 (CPQ,NUPSG,CI(IOC+1),SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
          END IF
 600  CONTINUE

      GOTO 999

 1700 CONTINUE
C EXCITING CASE, MIDLEV<IQ<IP
      DO 700 MVSGM=1,NMIDV
        DO 700 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 700
          ISYUC=MUL(ISYPQ,ISYUSG)
          NUPC=NOW(1,ISYUC,MVSGM)
          IF (NUPC.EQ.0) GOTO 700
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          NDWNSG=NOW(2,ISYDSG,MVSGM)
          IOC=IOCSF(ISYUC,MVSGM,ISYCI)
          INDEO=2*NLEV+(IP*(IP-1))/2+IQ
          NCP=NOCP(INDEO,ISYUC,MVSGM)
          IF(NCP.GT.0) THEN
            LICP=1+IOCP(INDEO,ISYUC,MVSGM)
C CASE IS: UPPER HALF, EXCITE:
            CALL EXC2_CP2 (CPQ,NDWNSG,NUPC,CI(IOC+1),NUPSG,SGM(ISGSTA),
     &             NCP,ICOUP(1,LICP),VTAB)
          END IF
 700  CONTINUE

      GOTO 999

 1800 CONTINUE
C EXCITING CASE, IQ<=MIDLEV<IP
C ALLOCATE TEMPORARY WORK AREA:
      DO 800 MVSGM=1,NMIDV
        MV1=MVL(MVSGM,2)
        MV2=MVL(MVSGM,1)
        IF((MV1.EQ.0).AND.(MV2.EQ.0)) GOTO 800
        DO 801 ISYUSG=1,NSYM
          NS1=NOCSF(ISYUSG,MVSGM,ISYSGM)
          IF(NS1.EQ.0) GOTO 801
          ISGSTA=1+IOCSF(ISYUSG,MVSGM,ISYSGM)
          NUPSG=NOW(1,ISYUSG,MVSGM)
          ISYDSG=MUL(ISYUSG,ISYSGM)
          ISYUC=MUL(ISYP,ISYUSG)
          ISYDC=MUL(ISYQ,ISYDSG)
          IF(MV2.EQ.0) GOTO 799
          NUPC=NOW(1,ISYUC,MV2)
          IF(NUPC.EQ.0) GOTO 799
          NDWNC=NOW(2,ISYDC,MV2)
          IF(NDWNC.EQ.0) GOTO 799
          INDEO=IP
          NCP=NOCP(INDEO,ISYUC,MV2)
          IF(NCP.EQ.0) GOTO 799
          NTMP=NUPSG*NDWNC
C          CALL VCLR(WORK(LSGTMP),1,NTMP)
          CALL NEWVCLR(WORK(LSGTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUC,MV2)
          IOC=IOCSF(ISYUC,MV2,ISYCI)
C CASE IS: UPPER HALF, EXCITE:
          CALL EXC2_CP2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LSGTMP),
     &             NCP,ICOUP(1,LICP),VTAB)
          INDEO=IQ
          NCP=NOCP(INDEO,ISYDC,MV2)
          IF(NCP.EQ.0) GOTO 799
          LICP=1+IOCP(INDEO,ISYDC,MV2)
C CASE IS: LOWER HALF, EXCITE:
          X=1.0D00
          CALL EXC1_CP2 (X,NUPSG,WORK(LSGTMP),SGM(ISGSTA),
     &               NCP,ICOUP(1,LICP),VTAB)
  799     CONTINUE
          IF(MV1.EQ.0) GOTO 801
          NUPC=NOW(1,ISYUC,MV1)
          IF(NUPC.EQ.0) GOTO 801
          NDWNC=NOW(2,ISYDC,MV1)
          IF(NDWNC.EQ.0) GOTO 801
          INDEO=NLEV+IP
          NCP=NOCP(INDEO,ISYUC,MV1)
          IF(NCP.EQ.0) GOTO 801
          NTMP=NUPSG*NDWNC
C          CALL VCLR(WORK(LSGTMP),1,NTMP)
          CALL NEWVCLR(WORK(LSGTMP),NTMP)
          LICP=1+IOCP(INDEO,ISYUC,MV1)
          IOC=IOCSF(ISYUC,MV1,ISYCI)
C CASE IS: UPPER HALF, EXCITE:
          CALL EXC2_CP2 (CPQ,NDWNC,NUPC,CI(IOC+1),NUPSG,WORK(LSGTMP),
     &             NCP,ICOUP(1,LICP),VTAB)
          INDEO=NLEV+IQ
          NCP=NOCP(INDEO,ISYDC,MV1)
          IF(NCP.EQ.0) GOTO 801
          LICP=1+IOCP(INDEO,ISYDC,MV1)
C CASE IS: LOWER HALF, EXCITE:
          X=1.0D00
          CALL EXC1_CP2 (X,NUPSG,WORK(LSGTMP),SGM(ISGSTA),
     &               NCP,ICOUP(1,LICP),VTAB)
 801  CONTINUE
 800  CONTINUE

 999  CONTINUE
      CALL QEXIT('SIGMA1')
      RETURN
      END
