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
C     SUBROUTINE TRDNS2O(IVEC,JVEC,DPT2)
      SUBROUTINE SIGDER(IVEC,JVEC,SCAL)
      use Fockof
      use caspt2_gradient, only: LUSTD,idSDMat
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "sigma.fh"
#include "SysDef.fh"
#include "caspt2_grad.fh"
#if defined(_MOLCAS_MPP_) && defined(_GA_)
#include "global.fh"
#endif
      COMMON /CPLCAS/ IFCOUP(MXCASE,MXCASE)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      logical :: bStat
#endif
C
C     Work in the MO basis
C     We need both explicit and implicit overlap derivatives. The latter
C     comes from the derivative of the transformation matrix.
C
C     p,q: inactive or secondary
C     y,z: active (t,u)
C     a,b: internally contracted
C     T1_{px} * S1_{xy} f_{yz} * T2_{pz}
C     = T1pa*C1xa * S1xy * fyz * T2pb*C2zb
C     Derivative of S1:
C     = (T1Ct1)px * (T2Ct2*f)py * dS1xy/da
C     Derivative of C1 (or, Lagrangian multiplier in MO basis):
C     = T1pa*dC1xa/da * S1xy * fyz * T2pb*C2zb
C       ...
C     = -1/2 (T1Ct1)pu * dS1tu/da * (T2Ct2*f*S1C1*Ct1)pt
C     Derivative of C2 (or, Lagrangian multiplier in MO basis):
C     = T1pa*C1xa * S1xy * fyz * T2pb*dC2zb/da
C       ...
C     = -1/2 (T1Ct1St1*f*C2*Ct2)pt * (T2Ct2)pu * dS2tu/da
C
C     About IMLTOP for the SGM subroutine
C     With IMLTOP=0: the vector for the second argument has to be
C     contravariant form (T*C),
C     With IMLTOP=1: the vector for the first  argument has to be
C     covariant form (T*SC),
C
C
C     Allocate some matrices for storing overlap and transformation
C     derivatives. Here constructs these derivatives in the MO basis,
C     but not in the internally contracted basis.
C
      MaxLen = 0
      Do iCase = 1, 11
        Do iSym = 1, nSym
          nAS = nASUP(iSym,iCase)
          MaxLen = Max(MaxLen,nAS*nAS)
        End Do
      End Do

      Call GETMEM('WRK','ALLO','REAL',ipWRK,MaxLen)
      Call DCopy_(MaxLen,[0.0D+00],0,Work(ipWRK),1)

      Do iCase = 1, 11
        Do iSym = 1, nSym
          nAS = nASUP(iSym,iCase)
          idSDer = idSDMat(iSym,iCase)
          CALL DDAFILE(LuSTD,1,Work(ipWRK),nAS*nAS,idSDer)
        End Do
      End Do
      Call GETMEM('WRK','FREE','REAL',ipWRK,MaxLen)
C
C Enter coupling cases for non-diagonal blocks:
      DO J=1,NCASES
      DO I=1,NCASES
      IFCOUP(I,J)=0
      END DO
      END DO
      IFCOUP( 2, 1)= 1
      IFCOUP( 3, 1)= 2
      IFCOUP( 5, 1)= 3
      IFCOUP( 6, 1)= 4
      IFCOUP( 7, 1)= 5
      IFCOUP( 6, 2)= 6
      IFCOUP( 7, 3)= 7
      IFCOUP( 5, 4)= 8
      IFCOUP( 8, 4)= 9
      IFCOUP( 9, 4)=10
      IFCOUP(10, 4)=11
      IFCOUP(11, 4)=12
      IFCOUP( 6, 5)=13
      IFCOUP( 7, 5)=14
      IFCOUP(10, 5)=15
      IFCOUP(11, 5)=16
      IFCOUP(12, 5)=23
      IFCOUP(13, 5)=24
      IFCOUP(12, 6)=17
      IFCOUP(13, 7)=18
      IFCOUP(10, 8)=19
      IFCOUP(11, 9)=20
      IFCOUP(12,10)=21
      IFCOUP(13,11)=22

C If the G1 correction to the Fock matrix is used, then the
C inactive/virtual coupling elements (which are non-zero for the
C case of average CASSCF) cannot be used in the CASPT2 equations.
      IF(FOCKTYPE.EQ.'G1      ' .AND. (.NOT. G1SECIN)) THEN
        IFCOUP(12,5)=0
        IFCOUP(13,5)=0
      END IF


C Transform to standard representation:
      CALL PTRTOC(0,IVEC,IVEC) !! T*C (internally contracted -> MO)
      IF(IVEC.NE.JVEC) CALL PTRTOC(0,JVEC,JVEC)

C Set up non-diagonal blocks of Fock matrix:
C SVC: add transposed fock matrix blocks
      NFIT=0
      NFIA=0
      NFTA=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NS=NSSH(ISYM)
        IOFFIT(ISYM)=NFIT
        IOFFIA(ISYM)=NFIA
        IOFFTA(ISYM)=NFTA
        NFIT=NFIT+NA*NI
        NFIA=NFIA+NS*NI
        NFTA=NFTA+NS*NA
      END DO
      NFIT=NFIT+1
      NFIA=NFIA+1
      NFTA=NFTA+1

      Call mma_allocate(FIT_Full,NFIT,Label='FIT_Full')
      Call mma_allocate(FTI_Full,NFIT,Label='FTI_Full')

      Call mma_allocate(FIA_Full,NFIA,Label='FIA_Full')
      Call mma_allocate(FAI_Full,NFIA,Label='FAI_Full')

      Call mma_allocate(FTA_Full,NFTA,Label='FTA_Full')
      Call mma_allocate(FAT_Full,NFTA,Label='FAT_Full')

      IFIFA=0
      DO ISYM=1,NSYM
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NS=NSSH(ISYM)
        NO=NORB(ISYM)

        FIT(ISYM)%A(1:NA*NI) =>
     &     FIT_Full(IOFFIT(ISYM)+1:IOFFIT(ISYM)+NA*NI)
        FTI(ISYM)%A(1:NA*NI) =>
     &     FTI_Full(IOFFIT(ISYM)+1:IOFFIT(ISYM)+NA*NI)

        FIA(ISYM)%A(1:NS*NI) =>
     &     FIA_Full(IOFFIA(ISYM)+1:IOFFIA(ISYM)+NS*NI)
        FAI(ISYM)%A(1:NS*NI) =>
     &     FAI_Full(IOFFIA(ISYM)+1:IOFFIA(ISYM)+NS*NI)

        FTA(ISYM)%A(1:NS*NA) =>
     &     FTA_Full(IOFFTA(ISYM)+1:IOFFTA(ISYM)+NS*NA)
        FAT(ISYM)%A(1:NS*NA) =>
     &     FAT_Full(IOFFTA(ISYM)+1:IOFFTA(ISYM)+NS*NA)

        CALL FBLOCK(WORK(LFIFA+IFIFA),NO,NI,NA,NS,
     &              FIT(ISYM)%A(:),FTI(ISYM)%A(:),
     &              FIA(ISYM)%A(:),FAI(ISYM)%A(:),
     &              FTA(ISYM)%A(:),FAT(ISYM)%A(:))

        IFIFA=IFIFA+(NO*(NO+1))/2

      END DO

      CALL TIMING(CPU0,CPU,TIO0,TIO)
C
C     Is it possible to reduce to one loop? We have to compute bra and
C     ket overlap and bra and ket wavefunctions are not identical, so
C     it seems impossible to reduce?
C
      NLOOP=2
      DO 1000 ILOOP=1,NLOOP
        !! ILOOP1 : <T+lambda|H|T       >
        !! ILOOP2 : <T       |H|T+lambda>

C Loop over types and symmetry block of sigma vector:
      DO 300 ICASE1=1,11
*     DO 300 ICASE1=1,NCASES
        DO 301 ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) GOTO 301
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          NIN1=NINDEP(ISYM1,ICASE1)
          NSGM2=NIS1*NAS1
          IF(NSGM2.EQ.0) GOTO 301

          CALL GETMEM('SGM2','ALLO','REAL',LSGM2,NSGM2)
          CALL DCOPY_(NSGM2,[0.0D0],0,WORK(LSGM2),1)

          NSGM1=0
          LSGM1=1
          IF(ICASE1.EQ.1) THEN
            NSGM1=NASH(ISYM1)*NISH(ISYM1)
          ELSE IF(ICASE1.EQ.4) THEN
            NSGM1=NASH(ISYM1)*NSSH(ISYM1)
          ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
            NSGM1=NIS1
          END IF
          IF(NSGM1.GT.0) THEN
            CALL GETMEM('SGM1','ALLO','REAL',LSGM1,NSGM1)
            CALL DCOPY_(NSGM1,[0.0D0],0,WORK(LSGM1),1)
          END IF

          IMLTOP=0
          DO 200 ICASE2=ICASE1+1,NCASES
            IFC=IFCOUP(ICASE2,ICASE1)
            IF(IFC.EQ.0) GOTO 200
            DO 100 ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) GOTO 100
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NCX=NIS2*NAS2
              IF(NCX.EQ.0) GOTO 100

              CALL RHS_ALLO(NAS2,NIS2,lg_CX)
              CALL RHS_READ(NAS2,NIS2,lg_CX,ICASE2,ISYM2,IVEC)
              IF (IVEC.NE.JVEC .AND. ILOOP.EQ.2) THEN
                !! T = T + \lambda
                If (SCAL.ne.1.0D+00) CALL RHS_SCAL(NAS2,NIS2,lg_CX,SCAL)
                CALL RHS_ALLO(NAS2,NIS2,lg_V1)
                CALL RHS_READ(NAS2,NIS2,lg_V1,ICASE2,ISYM2,JVEC)
                CALL RHS_DAXPY(NAS2,NIS2,1.0D+00,lg_V1,lg_CX)
                CALL RHS_FREE(NAS2,NIS2,lg_V1)
              END IF
C SVC: for case H (12,13) we can now pass the distributed array ID to
C the SGM subroutines
              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                LCX=lg_CX
                XTST=RHS_DDOT(NAS2,NIS2,lg_CX,lg_CX)
              ELSE
                CALL GETMEM('CX','ALLO','REAL',LCX,NCX)
                CALL RHS_GET(NAS2,NIS2,lg_CX,WORK(LCX))
                CALL RHS_FREE(NAS2,NIS2,lg_CX)
                XTST=DDOT_(NCX,WORK(LCX),1,WORK(LCX),1)
              END IF

#ifdef _DEBUGPRINT_
              WRITE(6,*)' ISYM1,ICASE1:',ISYM1,ICASE1
              WRITE(6,*)' ISYM2,ICASE2:',ISYM2,ICASE2
              WRITE(6,*)' SIGMA calling SGM with IMLTOP=',IMLTOP
#endif
C Compute contribution SGM2 <- CX, and SGM1 <- CX  if any
              CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &                 WORK(LSGM1),LSGM2,LCX,iWORK(LLISTS))

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                CALL RHS_FREE(NAS2,NIS2,lg_CX)
              ELSE
                CALL GETMEM('CX','FREE','REAL',LCX,NCX)
              END IF
 100        CONTINUE
 200      CONTINUE

C-SVC: sum the replicate arrays:
          MAX_MESG_SIZE = 2**27
          DO LSGM2_STA=1,NSGM2,MAX_MESG_SIZE
            NSGM2_BLK=MIN(MAX_MESG_SIZE,NSGM2-LSGM2_STA+1)
            CALL GADSUM(WORK(LSGM2+LSGM2_STA-1),NSGM2_BLK)
          END DO

          IF (NSGM1.GT.0) THEN
            CALL GADSUM(WORK(LSGM1),NSGM1)
          END IF

C If there are 1-electron contributions, add them into the 2-el
C part (This requires a non-empty active space.)
          IF(NSGM1.GT.0) THEN
            FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
            IF (ICASE1.EQ.1) THEN
              CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LSGM2),
     &                  WORK(LSGM1))
            ELSE IF(ICASE1.EQ.4) THEN
              CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LSGM2),
     &                  WORK(LSGM1))
            ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
              CALL SPEC1D(IMLTOP,FACT,WORK(LSGM2),WORK(LSGM1))
            END IF

            CALL GETMEM('SGM1','FREE','REAL',LSGM1,NSGM1)
          END IF

C-SVC: no need for the replicate arrays any more, fall back to one array
          CALL RHS_ALLO (NAS1,NIS1,lg_SGM2)
          CALL RHS_PUT (NAS1,NIS1,lg_SGM2,WORK(LSGM2))
          CALL GETMEM('SGM2','FREE','REAL',LSGM2,NSGM2)

C Add to sigma array. Multiply by S to  lower index.
C         CALL RHS_ALLO(NAS1,NIS1,lg_SGMX)
C         CALL RHS_READ(NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
          IF (ICASE1.LE.11) THEN
            CALL RHS_ALLO(NAS1,NIS1,lg_CX)
            CALL RHS_READ(NAS1,NIS1,lg_CX,ICASE1,ISYM1,IVEC)
            If (IVEC.NE.JVEC .AND. ILOOP.EQ.1) Then
              !! T = T + \lambda
              If (SCAL.ne.1.0D+00) CALL RHS_SCAL(NAS1,NIS1,lg_CX,SCAL)
              CALL RHS_ALLO(NAS1,NIS1,lg_V1)
              CALL RHS_READ(NAS1,NIS1,lg_V1,ICASE1,ISYM1,JVEC)
              CALL RHS_DAXPY(NAS1,NIS1,1.0D+00,lg_V1,lg_CX)
              CALL RHS_FREE(NAS1,NIS1,lg_V1)
            End If

            Call GETMEM('SDER1','ALLO','REAL',ipSDER1,NAS1*NAS1)
            idSDer = idSDMat(iSym1,iCase1)
            CALL DDAFILE(LuSTD,2,Work(ipSDER1),nAS1*nAS1,idSDer)
C           Call DCopy_(NAS1*NAS1,[0.0D+00],0,Work(ipSDER1),1)

C           CALL RHS_SCAL(NAS1,NIS1,lg_SGM2,2.0d+00)
            Call C1S1DER(Work(ipSDER1))

            idSDer = idSDMat(iSym1,iCase1)
            CALL DDAFILE(LuSTD,1,Work(ipSDER1),nAS1*nAS1,idSDer)
            Call GETMEM('SDER1','FREE','REAL',ipSDER1,NAS1*NAS1)

            CALL RHS_FREE(NAS1,NIS1,lg_CX)
          END IF

*         IF(ICASE1.NE.12 .AND. ICASE1.NE.13) THEN
C           CALL RHS_STRANS(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX,
C    &                      ICASE1,ISYM1)
*         ELSE
*           CALL RHS_DAXPY(NAS1,NIS1,ALPHA,lg_SGM2,lg_SGMX)
*         END IF
          CALL RHS_FREE (NAS1,NIS1,lg_SGM2)

C Write SGMX to disk.
C         CALL RHS_SAVE (NAS1,NIS1,lg_SGMX,ICASE1,ISYM1,JVEC)
C         CALL RHS_FREE (NAS1,NIS1,lg_SGMX)
 301    CONTINUE
 300  CONTINUE

      IMLTOP=1
C Loop over types and symmetry block of CX vector:
      DO 600 ICASE1=1,11
*     DO 600 ICASE1=1,NCASES
        DO 601 ISYM1=1,NSYM
          IF(NINDEP(ISYM1,ICASE1).EQ.0) GOTO 601
          NIS1=NISUP(ISYM1,ICASE1)
          NAS1=NASUP(ISYM1,ICASE1)
          ND2=NIS1*NAS1
          IF(ND2.EQ.0) GOTO 601

          CALL RHS_ALLO (NAS1,NIS1,lg_D2)
          CALL RHS_SCAL (NAS1,NIS1,lg_D2,0.0D0)
C Contract S*CX to form D2. Also form D1 from D2, if needed.

          NCX=ND2
          CALL RHS_ALLO (NAS1,NIS1,lg_CX)
          CALL RHS_READ (NAS1,NIS1,lg_CX,ICASE1,ISYM1,IVEC)

          IF (IVEC.NE.JVEC .AND. ILOOP.EQ.1) THEN
            !! T = T + \lambda
            If (SCAL.ne.1.0D+00) CALL RHS_SCAL(NAS1,NIS1,lg_CX,SCAL)
            CALL RHS_ALLO(NAS1,NIS1,lg_V1)
            CALL RHS_READ(NAS1,NIS1,lg_V1,ICASE1,ISYM1,JVEC)
            CALL RHS_DAXPY(NAS1,NIS1,1.0D+00,lg_V1,lg_CX)
            CALL RHS_FREE(NAS1,NIS1,lg_V1)
          END IF

          IF(ICASE1.NE.12 .AND. ICASE1.NE.13) THEN
           CALL RHS_STRANS (NAS1,NIS1,1.0D+00,lg_CX,lg_D2,ICASE1,ISYM1)
          ELSE
           CALL RHS_DAXPY(NAS1,NIS1,1.0D+00,lg_CX,lg_D2)
          END IF
          CALL RHS_FREE (NAS1,NIS1,lg_CX)

          CALL GETMEM('D2','ALLO','REAL',LD2,ND2)
          CALL RHS_GET (NAS1,NIS1,lg_D2,WORK(LD2))
          CALL RHS_FREE (NAS1,NIS1,lg_D2)

          ND1=0
          LD1=1
          IMLTOP=1
          FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
          IF(ICASE1.EQ.1) THEN
            ND1=NASH(ISYM1)*NISH(ISYM1)
            IF(ND1.GT.0) THEN
              CALL GETMEM('D1','ALLO','REAL',LD1,ND1)
              CALL DCOPY_(ND1,[0.0D0],0,WORK(LD1),1)
              CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LD2),
     &                    WORK(LD1))
            END IF
          ELSE IF(ICASE1.EQ.4) THEN
            ND1=NASH(ISYM1)*NSSH(ISYM1)
            IF(ND1.GT.0) THEN
              CALL GETMEM('D1','ALLO','REAL',LD1,ND1)
              CALL DCOPY_(ND1,[0.0D0],0,WORK(LD1),1)
              CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LD2),
     &                    WORK(LD1))
            END IF
          ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
            ND1=NIS1
            IF(ND1.GT.0) THEN
              CALL GETMEM('D1','ALLO','REAL',LD1,ND1)
              CALL DCOPY_(ND1,[0.0D0],0,WORK(LD1),1)
              CALL SPEC1D(IMLTOP,FACT,WORK(LD2),WORK(LD1))
            END IF
          END IF

          !! No need to compute for ICASE2 = 12 and 13
          DO 500 ICASE2=ICASE1+1,11 !! NCASES
            IF(IFCOUP(ICASE2,ICASE1).EQ.0) GOTO 500
            DO 400 ISYM2=1,NSYM
              IF(NINDEP(ISYM2,ICASE2).EQ.0) GOTO 400
              NIS2=NISUP(ISYM2,ICASE2)
              NAS2=NASUP(ISYM2,ICASE2)
              NIN2=NINDEP(ISYM2,ICASE2)
              NSGMX=NIS2*NAS2
              IF(NSGMX.EQ.0) GOTO 400

              IF (ICASE2.EQ.12 .OR. ICASE2.EQ.13) THEN
                CALL RHS_ALLO(NAS2,NIS2,lg_SGMX)
                CALL RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
                LSGMX=lg_SGMX
              ELSE
                CALL GETMEM('SGMX','ALLO','REAL',LSGMX,NSGMX)
                CALL DCOPY_(NSGMX,[0.0D0],0,WORK(LSGMX),1)
              END IF

#ifdef _DEBUGPRINT_
              WRITE(6,*)' ISYM1,ICASE1:',ISYM1,ICASE1
              WRITE(6,*)' ISYM2,ICASE2:',ISYM2,ICASE2
              WRITE(6,*)' SIGMA calling SGM with IMLTOP=',IMLTOP
#endif
C Compute contribution SGMX <- D2, and SGMX <- D1  if any
              CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
     &                 WORK(LD1),LD2,LSGMX,iWORK(LLISTS))
C             If (iCase2.LE.11) Then
C             End If

              IF (ICASE2.NE.12 .AND. ICASE2.NE.13) THEN
                MAX_MESG_SIZE = 2**27
                DO LSGMX_STA=1,NSGMX,MAX_MESG_SIZE
                  NSGMX_BLK=MIN(MAX_MESG_SIZE,NSGMX-LSGMX_STA+1)
                  CALL GADSUM(WORK(LSGMX+LSGMX_STA-1),NSGMX_BLK)
                END DO
C               CALL RHS_ALLO(NAS2,NIS2,lg_SGMX)
C               CALL RHS_READ(NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
C               CALL RHS_ADD(NAS2,NIS2,lg_SGMX,WORK(LSGMX))
                !! do C2DER
                Call GETMEM('SDER2','ALLO','REAL',ipSDER2,NAS2*NAS2)
                idSDer = idSDMat(iSym2,iCase2)
                CALL DDAFILE(LuSTD,2,Work(ipSDER2),nAS2*nAS2,idSDer)

                Call C2DER(Work(ipSDER2))

                idSDer = idSDMat(iSym2,iCase2)
                CALL DDAFILE(LuSTD,1,Work(ipSDER2),nAS2*nAS2,idSDer)
                Call GETMEM('SDER2','FREE','REAL',ipSDER2,NAS2*NAS2)
                !!
C               CALL GETMEM('SGMX','FREE','REAL',LSGMX,NSGMX)
              END IF

C-SVC: no need for the replicate arrays any more, fall back to one array
C             CALL RHS_SAVE (NAS2,NIS2,lg_SGMX,ICASE2,ISYM2,JVEC)
              IF (ICASE2.EQ.12 .OR.ICASE2.EQ.13) THEN
                CALL RHS_FREE (NAS2,NIS2,lg_SGMX)
              ELSE
                CALL GETMEM('SGMX','FREE','REAL',LSGMX,NSGMX)
              END IF
 400        CONTINUE
 500      CONTINUE
          CALL GETMEM('D2','FREE','REAL',LD2,ND2)
          IF(ND1.GT.0) CALL GETMEM('D1','FREE','REAL',LD1,ND1)
 601    CONTINUE
 600  CONTINUE

 1000 CONTINUE
C
C
C
      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUSGM=CPUSGM+(CPU1-CPU0)
      TIOSGM=TIOSGM+(TIO1-TIO0)
C
      Call mma_deallocate(FIT_Full)
      Call mma_deallocate(FTI_Full)
      Call mma_deallocate(FIA_Full)
      Call mma_deallocate(FAI_Full)
      Call mma_deallocate(FTA_Full)
      Call mma_deallocate(FAT_Full)
      Do iSym = 1, nSym
         FIT(iSym)%A => Null()
         FTI(iSym)%A => Null()
         FIA(iSym)%A => Null()
         FAI(iSym)%A => Null()
         FTA(iSym)%A => Null()
         FAT(iSym)%A => Null()
      End Do

C Transform contrav C  to eigenbasis of H0(diag):
      CALL PTRTOSR(1,IVEC,IVEC)
      IF(IVEC.NE.JVEC) CALL PTRTOSR(1,JVEC,JVEC)

C 99  CONTINUE
      RETURN

      Contains
C
C-----------------------------------------------------------------------
C
      Subroutine C1S1DER(SDER)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension SDER(*)
C
C     (T2Ct2*f)py * (T1Ct1)pz * dS1yz/da
C     -1/2 (T2Ct2*f*S1*C1*Ct1)pt * (T1Ct1)pu * dS1tu/da
C
      !! initialize
C     CALL GETMEM('TMP2','ALLO','REAL',LTMP2,NVEC1)
C     CALL DCOPY_(NVEC1,[0.0D0],0,WORK(LTMP2),1)
C     CALL GETMEM('TMP1','ALLO','REAL',LTMP1,MAX(1,NWEC1))
C     IF(NWEC1.GT.0) THEN
C       CALL DCOPY_(NWEC1,[0.0D0],0,WORK(LTMP1),1)
C     END IF
C
      !! 1. T2*Ct2*f
C     IMLTOP=0
C     CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
C    &         WORK(LTMP1),LTMP2,LVEC2,iWORK(LLISTS))
C
C     IF(NWEC1.GT.0) THEN
C       FACT=1.0D00/(DBLE(MAX(1,NACTEL)))
C       IF (ICASE1.EQ.1) THEN
C         CALL SPEC1A(IMLTOP,FACT,ISYM1,WORK(LTMP2),
C    &              WORK(LTMP1))
C       ELSE IF(ICASE1.EQ.4) THEN
C         CALL SPEC1C(IMLTOP,FACT,ISYM1,WORK(LTMP2),
C    &              WORK(LTMP1))
C       ELSE IF(ICASE1.EQ.5.AND.ISYM1.EQ.1) THEN
C         CALL SPEC1D(IMLTOP,FACT,WORK(LTMP2),WORK(LTMP1))
C       END IF
C     END IF
C     CALL GETMEM('TMP1','FREE','REAL',LTMP1,MAX(1,NWEC1))
C
      !! Finalize the derivative of S1
      !! 2S. (T2Ct2*f) * T1Ct1
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      if (is_real_par()) then
        CALL GA_CREATE_STRIPED ('H',NAS1,NAS1,'SDER',lg_SDER)
        CALL GA_PUT(lg_SDER,1,NAS1,1,NAS1,SDER,NAS1)
        call GA_DGEMM ('N','T',NAS1,NAS1,NIS1,
     *                 2.0D+00,lg_CX,lg_SGM2,1.0D+00,lg_SDER)
      else
#endif
        Call DGEMM_('N','T',NAS1,NAS1,NIS1,
     *              2.0D+00,WORK(lg_CX),NAS1,WORK(lg_SGM2),NAS1,
     *              1.0D+00,SDER,NAS1)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      end if
#endif
C     do i = 1, nas1*nis1
C       write (*,'(i4,2f20.10)') ,i,work(lg_cx+i-1),work(lg_sgm2+i-1)
C     end do
C
      !! Next, the derivative of C1
      !! 2C. (T2Ct2*f) * S1*C1 (MO -> IC)
      !!     lg_T * lg_V2 -> lg_V1
      CALL RHS_ALLO(NIN1,NIS1,lg_V1)
      ITYPE=1
      !! LSGM2 is local quantity, so put this in GA?
      CALL RHS_SR2C (ITYPE,1,NAS1,NIS1,NIN1,lg_V1,lg_SGM2,ICASE1,ISYM1)
      !! 3C. (T2Ct2*f) * S1*C1 * Ct1 (IC -> MO)
      !!     lg_T * lg_V1 -> lg_V2
      ITYPE=0
      CALL RHS_SR2C (ITYPE,0,NAS1,NIS1,NIN1,lg_V1,lg_SGM2,ICASE1,ISYM1)
      CALL RHS_FREE(NIN1,NIS1,lg_V1)
C
      !! 4C. (T1Ct1*f) * (T2Ct2St2*f*C1*Ct1)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      if (is_real_par()) then
        call GA_DGEMM ('N','T',NAS1,NAS1,NIS1,
     *                -1.0D+00,lg_CX,lg_SGM2,1.0D+00,lg_SDER)
        CALL GA_GET(lg_SDER,1,NAS1,1,NAS1,SDER,NAS1)
        bStat = GA_destroy(lg_SDER)
      else
#endif
        Call DGEMM_('N','T',NAS1,NAS1,NIS1,
     *             -1.0D+00,WORK(lg_CX),NAS1,WORK(lg_SGM2),NAS1,
     *              1.0D+00,SDER,NAS1)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      end if
#endif
C
C     CALL GETMEM('TMP2','FREE','REAL',LTMP2,NVEC1)
C
      End Subroutine C1S1DER
C
C-----------------------------------------------------------------------
C
      Subroutine C2DER(SDER)
C
      Implicit Real*8 (A-H,O-Z)
C
      Dimension SDER(*)
C
C     -1/2 (T2Ct2)pu * dS2tu/da * (T1Ct1St1*f*C2*Ct2)pt
C
      !! initialize
C     CALL RHS_ALLO(NAS2,NIS2,LTMP)
C     CALL RHS_SCAL(NAS2,NIS2,LTMP,0.0D+00)
C
      !! 1. T1*Ct1*St1*f
C     IMLTOP=1
C     CALL SGM(IMLTOP,ISYM1,ICASE1,ISYM2,ICASE2,
C    &         WORK(LWEC1S),LVEC1S,LTMP,iWORK(LLISTS))
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      if (is_real_par()) then
        CALL GA_CREATE_STRIPED ('V',NAS2,NIS2,'SDER',lg_SGMX)
        CALL GA_PUT(lg_SGMX,1,NAS2,1,NIS2,WORK(LSGMX),NAS2)
      else
#endif
       lg_SGMX = LSGMX
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      end if
#endif
C
      !! For icase = 12 or 13, there is no need to transform,
      !! so LTMP is always replicated, but T for A and C are
      !! distributed?, so...
      CALL RHS_ALLO(NIN2,NIS2,lg_V2)
      !! 2. (T1Ct1St1*f) * C2 (MO -> IC; LTMP -> LTMP2)
      ITYPE=0
      CALL RHS_SR2C (ITYPE,1,NAS2,NIS2,NIN2,lg_V2,lg_SGMX,ICASE2,ISYM2)
      !! 3. (T1Ct1St1*f) * C2 * Ct2 (IC -> MO; LTMP2 -> LTMP)
      CALL RHS_SR2C (ITYPE,0,NAS2,NIS2,NIN2,lg_V2,lg_SGMX,ICASE2,ISYM2)
      CALL RHS_FREE(NIN2,NIS2,lg_V2)
C
      !! 4. (T2Ct2*f) * (T1Ct1St1*f*C2*Ct2)
      CALL RHS_ALLO(NAS2,NIS2,lg_SGM)
      CALL RHS_READ(NAS2,NIS2,lg_SGM,ICASE2,ISYM2,IVEC)
          IF (IVEC.NE.JVEC .AND. ILOOP.EQ.2) THEN
            !! T = T + \lambda
            If (SCAL.ne.1.0D+00) CALL RHS_SCAL(NAS2,NIS2,lg_SGM,SCAL)
            CALL RHS_ALLO(NAS2,NIS2,lg_V1)
            CALL RHS_READ(NAS2,NIS2,lg_V1,ICASE2,ISYM2,JVEC)
            CALL RHS_DAXPY(NAS2,NIS2,1.0D+00,lg_V1,lg_SGM)
            CALL RHS_FREE(NAS2,NIS2,lg_V1)
          END IF
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      if (is_real_par()) then
        CALL GA_CREATE_STRIPED ('H',NAS2,NAS2,'SDER',lg_SDER)
        CALL GA_PUT(lg_SDER,1,NAS2,1,NAS2,SDER,NAS2)
        call GA_DGEMM ('N','T',NAS2,NAS2,NIS2,
     *                -1.0D+00,lg_SGM,lg_SGMX,1.0D+00,lg_SDER)
        CALL GA_GET(lg_SDER,1,NAS2,1,NAS2,SDER,NAS2)
        bStat = GA_destroy(lg_SGMX)
        bStat = GA_destroy(lg_SDER)
      else
#endif
        Call DGEMM_('N','T',NAS2,NAS2,NIS2,
     *             -1.0D+00,WORK(lg_SGM),NAS2,WORK(LSGMX),NAS2,
     *              1.0D+00,SDER,NAS2)
#if defined(_MOLCAS_MPP_) && defined(_GA_)
      end if
#endif
      CALL RHS_FREE(NAS2,NIS2,lg_SGM)
C
C     CALL RHS_FREE(NAS2,NIS2,LTMP)
C
      End Subroutine C2DER
C
      End subroutine sigder
