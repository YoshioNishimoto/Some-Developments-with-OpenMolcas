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
* Copyright (C) 2023, BRUNO TENORIO                                    *
************************************************************************
*  SUBROUTINE TRAO_RT2M
*  PURPOSE: TRANSFORM RT2M TO AO BASIS WITH CMO1 AND CMO2.
************************************************************************
      SUBROUTINE TRAO_RT2M(NRT2M, RT2M, CMO1, CMO2, ISTATE, JSTATE,
     & TRAFO, AFLNM)

      IMPLICIT REAL*8 (A-H,O-Z)

      integer SYM12, NRT2M
      Real*8, Allocatable :: SclI(:,:),SclJ(:,:),SclL(:,:)
      Real*8, Allocatable :: RT2MAO(:,:,:)
      Real*8, Allocatable :: Scrv(:,:,:), Scrv1(:,:,:), Scrv2(:,:,:)
      REAL*8 RT2M(nRT2M),CMO1(NCMO),CMO2(NCMO)
      integer :: istcmo(8)
      Integer I,J,IJ,L,LA,LO,isy1,no1,a,b,ist,isci,iscj,iscl
      LOGICAL PRTEST,TRAFO
      INTEGER NOI,NAI,NII,Jx,Ix,NOJ,NAJ,NIJ,NOL,NAL,NIL,IA,IO
      INTEGER ISYI,ISYJ,ISYL,JO,JA,KPOS,SFNT
      INTEGER :: IOFFAZ(8), IOFFOZ(8), IOFFBZ(8)
      CHARACTER*3 NUM1,NUM2
      CHARACTER*10 AFLNM
      CHARACTER*17 FNM
#include "prgm.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "stdalloc.fh"
#include "Files.fh"
#include "Struct.fh"

      PRTEST=.FALSE.!prints the MO basis rT2M matrix

      IF (TRAFO) THEN
       LU=59
       LU=IsFreeUnit(LU)
       WRITE(NUM1,'(I3.3)') ISTATE
       WRITE(NUM2,'(I3.3)') JSTATE
       !AFLNM='rT2DM_TAO_'
       FNM=AFLNM//NUM1//'_'//NUM2
       CALL Molcas_Open(LU,FNM)
      END IF
C IOFFA=NR OF ACTIVE ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFAZ(1)=0
      DO I=1,NSYM-1
        IOFFAZ(I+1)=IOFFAZ(I)+NASH(I)
      END DO
C IOFFO=NR OF OCC ORBITALS IN PREVIOUS SYMMETRY BLOCKS.
      IOFFOZ(1)=0
      DO J=1,NSYM-1
        IOFFOZ(J+1)=IOFFOZ(J)+NOSH(J)
      END DO
C IOFFO=NR OF BASIS FUNC IN PREVIOUS SYMMETRY BLOCKS.
      IOFFBZ(1)=0
      DO K=1,NSYM-1
        IOFFBZ(K+1)=IOFFBZ(K)+NBASF(K)
      END DO

      IST=0
      DO ISY1=1,NSYM
        ISTCMO(ISY1)=IST
        NO1=NOSH(ISY1)
        NB1=NBASF(ISY1)
        IST=IST+NO1*NB1
      END DO
C============================================================

      SYM12=MUL(LSYM1,LSYM2)
      DO ISYI=1,NSYM
       NOI=NOSH(ISYI)
       NBI=NBASF(ISYI)
       NAI=NASH(ISYI)
       NII=NISH(ISYI)
       ISCI=ISTCMO(ISYI)
       IF(NOI.EQ.0) CYCLE 
       Call mma_allocate(SclI,NBI,NOI,Label='SclI')
       SclI=0.0D0
       SFNT=ISCI+1
       DO Ix=1,NOI ! CMO1 on I
           SclI(:,Ix)=CMO1(SFNT:SFNT+NBI-1)
           SFNT=SFNT+NBI
       END DO
       DO ISYJ=1,NSYM
        NOJ=NOSH(ISYJ)
        NBJ=NBASF(ISYJ)
        NAJ=NASH(ISYJ)
        NIJ=NISH(ISYJ)
        ISCJ=ISTCMO(ISYJ)
        IF(NOJ.EQ.0) CYCLE 
        Call mma_allocate(SclJ,NBJ,NOJ,Label='SclJ')
        SclJ=0.0D0
        SFNT=ISCJ+1
        DO Ix=1,NOJ ! CMO2 on J
            SclJ(:,Ix)=CMO2(SFNT:SFNT+NBJ-1)
            SFNT=SFNT+NBJ
        END DO
        DO ISYL=1,NSYM
         NOL=NOSH(ISYL)
         NBL=NBASF(ISYL)
         NAL=NASH(ISYL)
         NIL=NISH(ISYL)
         ISCL=ISTCMO(ISYL)
         IF(NOL.EQ.0) CYCLE 
         Call mma_allocate(SclL,NBL,NOL,Label='SclL')
         SclL=0.0D0
         SFNT=ISCL+1
         DO Ix=1,NOL ! CMO2 on L
           SclL(:,Ix)=CMO2(SFNT:SFNT+NBL-1)
           SFNT=SFNT+NBL
         END DO
          IF(MUL(ISYI,MUL(ISYJ,ISYL)).EQ.SYM12) THEN
            IF(NAI.EQ.0) CYCLE
            IF(NAJ.EQ.0) CYCLE
            IF(NAL.EQ.0) CYCLE

            Call mma_allocate(Scrv,NOI,NOJ,NOL,Label='Scrv')
            Call mma_allocate(Scrv1,NOI,NOJ,NBL,Label='Scrv1')
            Call mma_allocate(Scrv2,NOI,NBJ,NBL,Label='Scrv2')
            Call mma_allocate(rt2mao,NBI,NBJ,NBL,Label='RT2MAO')
            Scrv(:,:,:)=0.0D0
            Scrv1(:,:,:)=0.0D0
            Scrv2(:,:,:)=0.0D0
            rt2mao(:,:,:)=0.0D0

            IF (PRTEST) WRITE(6,'(A10,8I7,8I7,8I7,8I7)')
     & ' # sub-Block:',ISYI,ISYJ, ISYL,NOI*NOJ*NOL
            DO I=1,NOI
             IA=IOFFAZ(ISYI)+I-NII
             IO=IOFFOZ(ISYI)+I
             DO J=1,NOJ
              JA=IOFFAZ(ISYJ)+J-NIJ
              JO=IOFFOZ(ISYJ)+J
              DO L=1,NOL
               LA=IOFFAZ(ISYL)+L-NIL
               LO=IOFFOZ(ISYL)+L
               KPOS=IA+NASHT*((LA+NASHT*(JA-1))-1)
               IF((IA.LE.0).or.(JA.LE.0).or.(LA.LE.0)) THEN
                Scrv(I,J,L)=0.0D0
               ELSE
                Scrv(I,J,L)=RT2M(KPOS)
               END IF
               IF(ABS( Scrv(I,J,L) ).LT.1.0D-19) Scrv(I,J,L)=0.0D0
              IF (PRTEST) write(6,'(I7,I7,I7,E26.12)')
     & IO,JO,LO,RT2M(KPOS)
              END DO
             END DO
            END DO


            ! ijl,lc -> ijc
            do a=1,size(Scrv,1)
             do b=1,size(Scrv,2)
              Scrv1(a,b,:) = matmul(Scrv(a,b,:), Transpose(SclL(:,:)) )
             end do
            end do

            ! ijc,jb -> ibc
            do a=1,size(Scrv1,1)
             do b=1,size(Scrv1,3)
              Scrv2(a,:,b) = matmul(Scrv1(a,:,b), Transpose(SclJ(:,:)) )
             end do
            end do

            ! ai,ibc -> abc
            do a=1,size(Scrv2,2)
             do b=1,size(Scrv2,3)
              rt2mao(:,a,b) = matmul(SclI(:,:), Scrv2(:,a,b) )
             end do
            end do

            ! Now print -> RT2MAO
            IF (TRAFO) THEN
             WRITE(LU,'(A11,8I7,8I7,8I7,8I7)')'#Sub-Block:',
     &       ISYI,ISYJ,ISYL,NBI*NBJ*NBL
              DO I=1,NBI
               IA=IOFFAZ(ISYI)+I-NII
               IB=IOFFBZ(ISYI)+I
               DO J=1,NBJ
                JA=IOFFAZ(ISYJ)+J-NIJ
                JB=IOFFBZ(ISYJ)+J
                DO L=1,NBL
                 LA=IOFFAZ(ISYL)+L-NIL
                 LB=IOFFBZ(ISYL)+L
                 KPOS=IB+NBST*((LB+NBST*(JB-1))-1)
                 WRITE(LU,'(I7,I7,I7,E26.12)')
     &                       IB,JB,LB,rt2mao(I,J,L)
                END DO
               END DO
              END DO
            END IF

CTEST *********************************
            !Print test CMO1
            if (.FALSE.) then
             write(6,*)'CMO1 matrix',ISYI,NBI,NOI
             do j=1,size(SclI,1)
               write(6,*) SclI(j,:)
             end do
             write(6,*)'************************************************
     &********************************'
            end if

            Call mma_deallocate(Scrv2 )
            Call mma_deallocate(Scrv1 )
            Call mma_deallocate(Scrv  )
            Call mma_deallocate(rt2mao)

          END IF
        Call mma_deallocate(SclL )
        END DO
       Call mma_deallocate(SclJ )
       END DO
      Call mma_deallocate(SclI )
      END DO
      IF (TRAFO) CLOSE (LU)
      RETURN
      END
