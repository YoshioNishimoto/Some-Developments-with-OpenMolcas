!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019 Marjan Khamesian and Roland Lindh                 *
!***********************************************************************
      Subroutine twlprm(Zeta,P_x,P_y,Alpha,Beta,i_x,i_y,j_x,j_y,l,intg)
      Implicit None
!
!
!     Local Arrays and integers
!
  Integer :: n, m, k, q, i_x, i_y, j_x, j_y, p, t
  Integer :: a, b, v, w, c, d, g, f, l
  Complex*16 intg
  Real*8 :: Lambda, Theta, Zeta
  Real*8 :: a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2
  Real*8 :: g1, g2, h1, h2, k1, k2, l2, l3
  Real*8 :: A_x, A_y, B_x, B_y, Alpha, Beta, P_x, P_y
  Real*8, parameter :: pi = atan(1.0)*4.d0
  Complex*16 cx, Kappa, i, Omega, l1, Gamma

  intg = DCmplx(0.0D0,0.0D0)

  If (i_x.lt.0 .or. i_y.lt.0 .or. j_x.lt.0 .or. j_y.lt.0) Return

  i=DCmplx(0.0D0,1.0D0)

  !- Lambda -----

  Do n = 0, i_x
     a1 = choose(i_x,n)
     a2 = (P_x - A_x)**n
     Do m = 0, j_x
        b1 = choose(j_x, m)
        b2 = (P_x - B_x)**m

        a = i_x + j_x - n - m
        Do k = 0, i_y

           c1 = choose(i_y, k)
           c2 = (P_y - A_y)**k
           Do q = 0, j_y

              d1 = choose(j_y, q)
              d2 = (P_y - B_y)**q

              b = i_y + j_y - k - q

              Lambda = (a1*a2*b1*b2*c1*c2*d1*d2* &
                   exp((-(Alpha*Beta)/Zeta)*((A_x-B_x)**2+(A_y - B_y)**2)))

              !- Theta ------

              Do v = 0, a
                 e1 = choose(a,v)
                 e2 = (P_x)**v
                 Do w = 0, b
                    f1 = choose(b,w)
                    f2 = (P_y)**w

                    Theta = ( e1 * e2 * f1 * f2 )

                    ! - Omega ------

                    f = a-v+b-w+1

                    Do c = 0, a-v
                       g1 = choose(a-v, c)
                       g2 = (1/2)**(a-v)
                       Do d = 0, b-w
                          h1 = choose(b-w, d)
                          h2 = (1/(i))**(b-w)

                          g = a+b+l-v-w-(2*c)-(2*d)

                          cx = (i*g) / ( 2*Zeta*P_y )
                          Kappa = ( exp(-Zeta * P_y**2.0) * (exp(2*pi*i*g) - 1) ) / (2.0*Zeta*P_y)

                          Omega = Kappa *(g1*g2*h1*h2)

                          !- Gamma --------

                          Do t = 0, f
                             Do p = 0, t-1

                                k1 = choose(f,t)
                                k2 = choose(t-1,p)
                                l1 = (-cx)**(f-t) * (cx + P_x)**(t-1-p) * (-1/2)**p *(1/(2*sqrt(2*pi)))
                                l2 = integral_gauss(Zeta,p)
                                Gamma = ( k1 * k2 * l1 * l2 )

                                intg = intg + Lambda * Theta * Omega * Gamma

                             End Do ! p
                          End Do    ! t

                       End Do ! d
                    End Do    ! c

                 End Do ! w
              End Do    ! v

           End Do ! q
        End Do    ! k

     End Do !m
  End Do    !n

!------------------------------------------------------------------
contains


  function factorial (n) result (res)

    implicit none
    integer, intent (in) :: n
    integer :: res
    integer :: i

    res = product ((/(i, i = 1, n)/))

  end function factorial
!!!!!!!!!!!!!!!!!!!!
  function choose (n, k) result (res)

    implicit none
    integer, intent (in) :: n
    integer, intent (in) :: k
    integer :: res

    res = factorial (n) / (factorial (k) * factorial (n - k))

  end function choose

!!!!!!!!!!!!!!!!!!!!
  function integral_gauss (Zeta,p) result (res)

    implicit none
    integer, intent (in) :: p
    real*8, intent (in) :: Zeta
    real*8 :: res
    !    real, intent (out) :: res

    res = (0.5) * (-0.5)**p * Zeta**(0.5-p)


  end function integral_gauss

End Subroutine twlprm
