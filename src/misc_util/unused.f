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
      SUBROUTINE Unused_real(r)
      IMPLICIT NONE
      REAL*8 r,r2
      IF (.FALSE.) r2=r
      END SUBROUTINE Unused_real

      SUBROUTINE Unused_real_array(r)
      IMPLICIT NONE
      REAL*8 r(*),r2
      IF (.FALSE.) r2=r(1)
      END SUBROUTINE Unused_real_array

      SUBROUTINE Unused_integer(i)
      IMPLICIT NONE
      INTEGER i,i2
      IF (.FALSE.) i2=i
      END SUBROUTINE Unused_integer

      SUBROUTINE Unused_integer_array(i)
      IMPLICIT NONE
      INTEGER i(*),i2
      IF (.FALSE.) i2=i(1)
      END SUBROUTINE Unused_integer_array

      SUBROUTINE Unused_logical(l)
      IMPLICIT NONE
      LOGICAL l,l2
      IF (.FALSE.) l2=l
      END SUBROUTINE Unused_logical

      SUBROUTINE Unused_logical_array(l)
      IMPLICIT NONE
      LOGICAL l(*),l2
      IF (.FALSE.) l2=l(1)
      END SUBROUTINE Unused_logical_array

      SUBROUTINE Unused_character(c)
      IMPLICIT NONE
      CHARACTER c(*),c2
      IF (.FALSE.) c2=c(1)
      END SUBROUTINE Unused_character
