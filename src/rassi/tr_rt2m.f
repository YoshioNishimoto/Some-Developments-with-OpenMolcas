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
* Copyright (C) 2020, BRUNO TENORIO                                    *
************************************************************************
*  SUBROUTINE TR_RT2M
*  PURPOSE: TRANSFORM RT2M USING L, THE CHOLESKY FACTORIZATION OF CMO2.
*  THE L MATRIX IS OBTAINED FROM THE SUBROUTINE TRORB_LL
************************************************************************
      SUBROUTINE TR_RT2M(NRT2M, RT2M, LMAT, S1MAT)

      IMPLICIT REAL*8 (A-H,O-Z)

      integer SYM12, NRT2M
      Real*8, Allocatable :: SclI(:,:),SclJ(:,:),SclL(:,:)
      Real*8, Allocatable :: S1I(:,:)
      Real*8, Allocatable :: Scrv(:,:,:), Scrv2(:,:,:)
      REAL*8 RT2M(nRT2M)
      REAL*8 LMAT(NTDMAB), S1MAT(NTDMAB)
      integer :: istacc(8)
      Integer I,J,IJ,L,LA,LO,isy1,no1,a,b,iscc,isci,iscj,iscl
      LOGICAL PRTEST
      INTEGER NOI,NAI,NII,Jx,Ix,NOJ,NAJ,NIJ,NOL,NAL,NIL,IA,IO
      INTEGER ISYI,ISYJ,ISYL,JO,JA,KPOS
      INTEGER :: IOFFAZ(8), IOFFOZ(8)
#include "prgm.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "stdalloc.fh"
#include "Files.fh"
#include "Struct.fh"

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

      ISCC=0
      DO ISY1=1,NSYM
        ISTACC(ISY1)=ISCC
        NO1=NOSH(ISY1)
        ISCC=ISCC+NO1*NO1
      END DO
C============================================================
      PRTEST=.FALSE.

      SYM12=MUL(LSYM1,LSYM2)
      !WRITE(6,*)'# sub-Block info:Sym(I,J,L), NumOrb in SymmBlock'

      DO ISYI=1,NSYM
       NOI=NOSH(ISYI)
       NAI=NASH(ISYI)
       NII=NISH(ISYI)
       ISCI=ISTACC(ISYI)
       IF(NOI.EQ.0) CYCLE !GOTO 273
       Call mma_allocate(SclI,NOI,NOI,Label='SclI')
       Call mma_allocate(S1I,NOI,NOI,Label='S1I')
       SclI=0.0D0
       S1I =0.0D0
       DO Ix=1,NOI ! Lmat on I
         DO Jx=1,NOI
           IJ=Ix+(NOI*(Jx-1))
           IF(Jx.LE.Ix) SclI(Ix,Jx)=LMAT(IJ+ISCI)
           S1I(Ix,Jx)=S1MAT(IJ+ISCI) ! S(CMO1) ** -1
         END DO
       END DO
       DO ISYJ=1,NSYM
        NOJ=NOSH(ISYJ)
        NAJ=NASH(ISYJ)
        NIJ=NISH(ISYJ)
        ISCJ=ISTACC(ISYJ)
        IF(NOJ.EQ.0) CYCLE !GOTO 373
        Call mma_allocate(SclJ,NOJ,NOJ,Label='SclJ')
        SclJ=0.0D0
        DO Ix=1,NOJ ! Lmat on J
          DO Jx=1,NOJ
            IJ=Ix+(NOJ*(Jx-1))
            IF(Jx.LE.Ix) SclJ(Ix,Jx)=LMAT(IJ+ISCJ)
          END DO
        END DO
        DO ISYL=1,NSYM
         NOL=NOSH(ISYL)
         NAL=NASH(ISYL)
         NIL=NISH(ISYL)
         ISCL=ISTACC(ISYL)
         IF(NOL.EQ.0) CYCLE !GOTO 473
         Call mma_allocate(SclL,NOL,NOL,Label='SclL')
         SclL=0.0D0
         DO Ix=1,NOL ! Lmat on L
           DO Jx=1,NOL
             IJ=Ix+(NOL*(Jx-1))
             IF(Jx.LE.Ix) SclL(Ix,Jx)=LMAT(IJ+ISCL)
           END DO
         END DO
          IF(MUL(ISYI,MUL(ISYJ,ISYL)).EQ.SYM12) THEN
            IF(NAI.EQ.0) CYCLE! GOTO 673
            IF(NAJ.EQ.0) CYCLE! GOTO 673
            IF(NAL.EQ.0) CYCLE! GOTO 673

            Call mma_allocate(Scrv,NOI,NOJ,NOL,Label='Scrv')
            Call mma_allocate(Scrv2,NOI,NOJ,NOL,Label='Scrv2')
            Scrv(:,:,:)=0.0D0
            Scrv2(:,:,:)=0.0D0

!            WRITE(6,'(A10,8I7,8I7,8I7,8I7)')' # sub-Block:',ISYI,ISYJ,
!     &      ISYL,NOI*NOJ*NOL
            DO I=1,NOI
             IA=IOFFAZ(ISYI)+I-NII
             IO=IOFFOZ(ISYI)+I
             DO J=1,NOJ
              JA=IOFFAZ(ISYJ)+J-NIJ
              JO=IOFFOZ(ISYJ)+J
              DO L=1,NOL
               LA=IOFFAZ(ISYL)+L-NIL
               LO=IOFFOZ(ISYL)+L
               IF((IA.LE.0).or.(JA.LE.0).or.(LA.LE.0)) THEN
!                write(6,'(I7,I7,I7,E26.12)') IO,JO,LO,0.0D0
               ELSE
                KPOS=IA+NASHT*((LA+NASHT*(JA-1))-1)
!                IF(ABS(RT2M(KPOS)).LT.1.0D-19) THEN
!                 RT2M(KPOS) = 0.0D0
!                END IF
!                write(6,'(I7,I7,I7,E26.12)') IO,JO,LO,RT2M(KPOS)
                Scrv(I,J,L)=RT2M(KPOS)
               END IF
              END DO
             END DO
            END DO

            ! S_CMO1^-1 * r2TDM
            ! ai,ibc -> abc
            Scrv2(:,:,:) = 0.0D0
            do a=1,size(Scrv,2)
             do b=1,size(Scrv,3)
              Scrv2(:,a,b) = matmul(Transpose(S1I(:,:)), Scrv(:,a,b) )
             end do
            end do

            ! ijl,lc -> ijc
            Scrv(:,:,:) = 0.0D0
            do a=1,size(Scrv2,1)
             do b=1,size(Scrv2,2)
              Scrv(a,b,:) = matmul(Scrv2(a,b,:), SclL(:,:) )
             end do
            end do

            ! ijc,jb -> ibc
            Scrv2(:,:,:) = 0.0D0
            do a=1,size(Scrv,1)
             do b=1,size(Scrv,3)
              Scrv2(a,:,b) = matmul(Scrv(a,:,b), SclJ(:,:) )
             end do
            end do

            ! ai,ibc -> abc
            Scrv(:,:,:) = 0.0D0
            do a=1,size(Scrv2,2)
             do b=1,size(Scrv2,3)
              Scrv(:,a,b) = matmul(Transpose(SclI(:,:)), Scrv2(:,a,b) )
             end do
            end do

            ! Now Scr2 -> RT2M
            IF (PRTEST) WRITE(6,'(A12,8I7,8I7,8I7,8I7)')'# Sub-Block:',
     &       ISYI,ISYJ,ISYL,NOI*NOJ*NOL
            DO I=1,NOI
             IA=IOFFAZ(ISYI)+I-NII
             IO=IOFFOZ(ISYI)+I
             DO J=1,NOJ
              JA=IOFFAZ(ISYJ)+J-NIJ
              JO=IOFFOZ(ISYJ)+J
              DO L=1,NOL
               LA=IOFFAZ(ISYL)+L-NIL
               LO=IOFFOZ(ISYL)+L
               IF((IA.LE.0).or.(JA.LE.0).or.(LA.LE.0)) THEN
                IF (PRTEST) write(6,'(I7,I7,I7,E26.12)')IO,JO,LO,0.0D0
               ELSE
                KPOS=IA+NASHT*((LA+NASHT*(JA-1))-1)
                RT2M(KPOS)=Scrv(I,J,L)
                !IF(ABS(RT2M(KPOS)).LT.1.0D-19) THEN
                ! RT2M(KPOS) = 0.0D0
                !END IF
                IF (PRTEST) write(6,'(I7,I7,I7,E26.12)')
     &                       IO,JO,LO,RT2M(KPOS)
               END IF
              END DO
             END DO
            END DO

CTEST *********************************
            !Print Lmat
            if (.FALSE.) then
             write(6,*)'Lower Triangular L matrix',ISYl,NOL
             do j=1,size(SclL,1)
               write(6,*) SclL(j,:)
             end do
             write(6,*)'************************************************
     &********************************'
            end if

            Call mma_deallocate(Scrv2 )
            Call mma_deallocate(Scrv  )

          END IF
        Call mma_deallocate(SclL )
        END DO
       Call mma_deallocate(SclJ )
       END DO
      Call mma_deallocate(SclI )
      Call mma_deallocate(S1I )
      END DO

      RETURN
      END
