!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!subroutine SQPRT(A,N)

!use Constants, only: Zero
!use Definitions, only: wp, iwp, u6

!implicit none
!integer(kind=iwp), intent(in) :: N
!real(kind=wp), intent(in) :: A(N,N)
!character(len=60) :: FRMT
!integer(kind=iwp) :: I, J
!real(kind=wp) :: BIG

!BIG = Zero
!do I=1,N
!  do J=1,N
!    BIG = max(BIG,abs(A(I,J)))
!  end do
!end do
!if ((0.1_wp < BIG) .and. (BIG < 1.0e4_wp)) then
!  FRMT = '(8(1X,F12.6))'
!else
!  FRMT = '(8(1X,E12.6))'
!end if
!do I=1,N
!  write(u6,FRMT) A(I,:)
!end do

!return

!end subroutine SQPRT
      SUBROUTINE SQPRT(V,N)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION V(N,N)
      MAX = 5
      IMAX = 0
  100 IMIN = IMAX+1
      IMAX = IMAX+MAX
      IF (IMAX .GT. N) IMAX = N
      WRITE (6,9008)
      WRITE (6,9028) (I,I = IMIN,IMAX)
      WRITE (6,9008)
      DO 120 J = 1,N
  120 WRITE (6,9048) J,(V(J,I),I = IMIN,IMAX)
      IF (IMAX .LT. N) GO TO 100
      RETURN
 9008 FORMAT(1X)
 9028 FORMAT(6X,10(4X,I4,4X))
 9048 FORMAT(I5,1X,10F12.7)
      END
