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
*  SUBROUTINE TRORB_LL
*  PURPOSE: ORTHNORMALIZE CMO2 USING A= L * L**T CHOLESKY FACTORIZATION
*  A=L*L**T WHERE A=C**T * S * C
*  AND U=C*(L**-1)**T IS ORTHONORMAL
************************************************************************
      SUBROUTINE TRORB_LL(CMO1,UMO2,LMAT,S1MAT)

      Implicit none

      integer :: isym
      character(len=8) :: Label
      integer :: nb, nbast, nbast1, nbast2
      real*8, allocatable :: SAO(:), IAO(:), TMPIAO(:)
      real*8, allocatable :: Scr(:), Scr22(:), Scr2(:),Scr1(:)
      real*8, allocatable :: Lmatinv(:), CMOA(:)
      integer :: iOff1, iOff2
      integer :: iOpt,iSyLbl,iRc
      integer :: IC,istca,istcb,ist,ista,istcc,istc
      REAL*8 CMO1(nCMO),UMO2(nCMO),LMAT(NTDMAB),S1MAT(NTDMAB)
      INTEGER :: INFOL
      integer :: istcmo(8), istao(8), istacc(8)
      Integer ::  no1,nb1,isy1,I,J,indij,nldm
      LOGICAL PRTEST
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
#include "rassi.fh"
#include "stdalloc.fh"

C============================================================
      PRTEST=.False.
      nbast=0
      nbast1=0
      nbast2=0
      do isym=1,nsym
        nb=NBASF(isym)
        nbast=nbast+nb
        nbast1=nbast1+(nb*(nb+1))/2
        nbast2=nbast2+nb**2
      end do

      call mma_allocate(CMOA,nCMO)
      CALL DCOPY_(nCMO,UMO2,1,CMOA,1) !CMO2

      call mma_allocate(SAO,NBAST1)
      call mma_allocate(IAO,NBAST2)
      IAO=0.0D0

      iRc=-1
      iOpt=6
      IC=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,IC,SAO,iSyLbl)
      iOff1 = 0
      iOff2 = 0
      Do iSym = 1,nSym
        nb = nBasf(iSym)
        If ( nb.gt.0 ) then
          call mma_allocate(TMPIAO,nb*nb)
          TMPIAO=0.0D0
          Call Square(SAO(1+iOff1),TMPIAO,1,nb,nb)
          CALL DCOPY_(nb*nb,TMPIAO,1,IAO(iOff2+1),1)
          call mma_deallocate(TMPIAO)
        end if
        iOff1 = iOff1 + (nb**2 + nb)/2
        iOff2 = iOff2 + (nb**2)
      end do
      call mma_deallocate(SAO)

C============================================================

      IST=1
      ISTA=1
      ISTCC=1
      DO ISY1=1,NSYM
        ISTCMO(ISY1)=IST
        ISTAO(ISY1)=ISTA
        ISTACC(ISY1)=ISTCC
        NO1=NOSH(ISY1)
        IST=IST+NO1*NBASF(ISY1)
        ISTA=ISTA+(NBASF(ISY1)*NBASF(ISY1) )
        ISTCC=ISTCC+NO1*NO1
      END DO

      DO ISY1=1,NSYM  ! START OF THE SYMM LOOP
        ISTCB=ISTCMO(ISY1)
        ISTCA=ISTAO(ISY1)
        ISTC=ISTACC(ISY1)
        NO1=NOSH(ISY1)
        NB1=NBASF(ISY1)
        NLDM=NO1*NO1
        IF(NB1*NO1.EQ.0) CYCLE !GOTO 15

        call mma_allocate(scr,NB1*NO1)
        call mma_allocate(scr2,NLDM)
        call mma_allocate(scr22,NLDM)
        call mma_allocate(scr1,NO1*(NO1+1)/2)
        Scr(:)=0.0D0
        Scr2(:)=0.0D0
        Scr22(:)=0.0D0
        Scr1(:)=0.0D0

        CALL DGEMM_('N','N', NB1, NO1, NB1, 1.0D0,
     &                 IAO(ISTCA),NB1, CMO1(ISTCB), NB1,
     &         0.0D0, Scr, NB1)

        CALL DGEMM_('T','N', NO1, NO1, NB1, 1.0D0,
     &                 CMO1(ISTCB),NB1, Scr, NB1,
     &         0.0D0, Scr2, NO1)
        CALL DCOPY_(NO1*NO1,Scr2,1,Scr22,1) ! save S in Scr22

        IF (PRTEST) write(6,*)'TEST S(CMO1)'
        IF (PRTEST) call print_matrix(Scr22, no1)

        INFOL=0
        CALL DPOTRF_('L',NO1, Scr2, NO1, INFOL)

        DO I=1,NO1 ! COMPACT LOW TRIANGULAR OF S
         DO J=1,I
          INDIJ= I + ((J-1)*(2*NO1 - J)/2)
          Scr1(INDIJ)=Scr2(I+(NO1*(J-1)))
         END DO
        END DO

        CALL DPPTRI('L', NO1, Scr1, INFOL ) ! -> S ** -1

        Scr2(:)=0.0D0
        DO I=1,NO1 ! NOW SQUARE S ** -1
         DO J=1,I
          INDIJ= I + ((J-1)*(2*NO1 - J)/2)
          Scr2(I+(NO1*(J-1)))=Scr1(INDIJ)
          Scr2(J+(NO1*(I-1)))=Scr1(INDIJ)
         END DO
        END DO

        CALL DCOPY_(NO1*NO1,Scr2,1,S1MAT(ISTC),1)

        IF (PRTEST) THEN
          write(6,*)'TEST S(CMO1) ** -1 matrix'
          call print_matrix(S1MAT(ISTC), no1)
          Scr2(:)=0.0D0
          write(6,*)'TEST S(CMO1) ** -1 * S(CMO1)'
          CALL DGEMM_('N','T', NO1, NO1, NO1, 1.0D0,
     &                 S1MAT(ISTC), NO1, Scr22, NO1,
     &         0.0D0, Scr2, NO1)
          call print_matrix(Scr2, no1)
        END IF
        Call mma_deallocate(Scr22)

! ********************************************

        Scr(:)=0.0D0
        Scr2(:)=0.0D0
        CALL DGEMM_('N','N', NB1, NO1, NB1, 1.0D0,
     &                 IAO(ISTCA),NB1, CMOA(ISTCB), NB1,
     &         0.0D0, Scr, NB1)

        CALL DGEMM_('T','N', NO1, NO1, NB1, 1.0D0,
     &                 CMOA(ISTCB),NB1, Scr, NB1,
     &         0.0D0, Scr2, NO1)

! Proceed to compute A = L*L**T
! WHERE L IS LOWER TRIANGULAR
! DPOTRF computes the Cholesky factorization of a real symmetric
! positive definite matrix A = C**T * S * C

        INFOL=0
        IF (PRTEST) THEN
          write(6,*)'TEST A(= S(CMO2) ) matrix'
          call print_matrix(Scr2, no1)
        END IF

        CALL DPOTRF_('L',NO1, Scr2, NO1, INFOL)
        ! L and L**-1 are square matrices of dimention (NO1,NO1)

        call mma_allocate(Lmatinv,NLDM)
        Lmatinv(:)=0.0D0

        DO I=1,NO1
         DO J=1,NO1
          INDIJ=J+(NO1*(I-1))
          IF(I.LE.J) Lmat(ISTC-1+INDIJ)=Scr2(INDIJ)
         END DO
        END DO
        CALL DCOPY_(NLDM,Lmat(ISTC),1,LMATINV,1)

        IF (PRTEST) THEN
          write(6,*)'TEST L matrix'
          call print_matrix(LMAT(ISTC), no1)

          !BRNCAT TEST A = L * L**T
          Scr2(:)=0.0D0
          CALL DGEMM_('N','T', NO1, NO1, NO1, 1.0D0,
     &                 LMAT(ISTC), NO1, LMAT(ISTC), NO1,
     &         0.0D0, Scr2, NO1)

          WRITE(6,*)'TEST A = L * L**T'
          call print_matrix(Scr2, no1)
        END IF

        ! DEFINE UMO2 = CMO2 * (L**-1)**T
        CALL DTRTRI('L', 'N', NO1, LMATINV, NO1, INFOL ) ! -> LMAT**-1

        IF (PRTEST) THEN
        ! TEST LMAT**-1 * LMAT
          Scr2(:)=0.0D0
          CALL DGEMM_('N','N', NO1, NO1, NO1, 1.0D0,
     &                 LMAT(ISTC), NO1, LMATINV, NO1,
     &         0.0D0, SCR2, NO1)
          WRITE(6,*)'TEST  LMAT**-1 * LMAT'
          call print_matrix(Scr2, no1)
        END IF

        CALL DGEMM_('N','T', NB1, NO1, NO1, 1.0D0,
     &                 CMOA(ISTCB), NB1, LMATINV, NO1,
     &         0.0D0, UMO2(ISTCB), NB1) ! -> CMO2 * LMATINV**T

        IF (PRTEST) THEN
        ! TEST OVERLAP U * U**T
          Scr(:)=0.0D0
          Scr2(:)=0.0D0
          CALL DGEMM_('N','N', NB1, NO1, NB1, 1.0D0,
     &                 IAO(ISTCA),NB1, UMO2(ISTCB), NB1,
     &         0.0D0, Scr, NB1)

          CALL DGEMM_('T','N', NO1, NO1, NB1, 1.0D0,
     &                 UMO2(ISTCB),NB1, Scr, NB1,
     &         0.0D0, Scr2, NO1)

          WRITE(6,*)'TEST  U * U**T'
          call print_matrix(Scr2, no1)
        END IF

        Call mma_deallocate(Scr)
        Call mma_deallocate(Scr2)
        Call mma_deallocate(Scr1)
        Call mma_deallocate(Lmatinv)
!15      CONTINUE
      END DO
      Call mma_deallocate(IAO)
      call mma_deallocate(CMOA)
      RETURN

      END SUBROUTINE TRORB_LL

      subroutine print_matrix(arr, n)
      implicit none
      integer, intent(in) :: n
      real*8 :: arr(n*n), mat(n,n)
      integer :: i
      CHARACTER*10 FWN
      CHARACTER*1 NU
      WRITE(NU,'(I1.1)') n
      FWN='('//NU//'(f16.7))'

      mat(:,:) = reshape(arr,(/ n,n /))
      do i = 1, n
        write(6,FWN) mat(i,:)
      enddo
      end subroutine print_matrix

