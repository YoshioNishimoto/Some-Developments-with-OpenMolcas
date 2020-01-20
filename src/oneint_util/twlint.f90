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
      SubRoutine TWLInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,   &
                       Final,nZeta,nIC,nComp,la,lb,A,RB,nRys,          &
                       Array,nArr,kVector,nOrdOp,lOper,iChO,           &
                       iStabM,nStabM,                                  &
                       PtChrg,nGrid,iAddPot)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of integrals for the      *
!         interaction between matter and light with orbital angular    *
!         momentum.                                                    *
!***********************************************************************
      Implicit None
!
!     External Arrays and integers
!
      Integer nZeta, la, lb, nIC, nAlpha, nBeta, nArr, nComp, nOrdOp,  &
              nStabM
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),         &
             Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),     &
             rKappa(nZeta), P(nZeta,3), A(3), RB(3), TC(3),         &
             Array(nZeta*nArr), kvector(3)
      Integer iStabM(0:nStabM-1), iDCRT(0:7),                          &
              iStabO(0:7), lOper(nComp), iChO(nComp)
      Integer nRys, nGrid, iAddPot
      Real*8  PtChrg
!
!     Local Arrays and integers
!
      Integer iZeta, iAlpha, iBeta, i_x, i_y, i_z, j_x, j_y, j_z,      &
              lAng
      Real*8 Value
!
      lAng= 1             ! Temporary set value

!
!     We need here a transformation of the Cartesian coordinates in which
!     the k-vector coinsides with the z-vector direction.
!     In particular, we will transform P to the new coordinate system.
!
!     ... more to come ...
!
!     The integration over the z-subspace is done using the normal code
!     for the exponetial operator.
!
      Do iBeta = 1, nBeta
         Do iAlpha = 1, nAlpha
            iZeta = nAlpha*(iBeta-1) + iAlpha
!
            Do i_x =  0, la
               Do i_y = 0, la-i_x
                  i_z = la - i_x - i_y

            Do j_x =  0, lb
               Do j_y = 0, lb-j_x
                  j_z = lb - j_x - j_y

!                 Compute the x-y part of the integral

                  Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                              Alpha(iAlpha), Beta(iBeta),              &
                              i_x,i_y,                                 &
                              j_x,j_y, lAng, Value)

!                 Compute the z part of the integral

!                 ... more to come ...

!                 Assemble to the complete integral

!                 ... more to come ...


               End Do ! j_y
            End Do    ! j_x

               End Do    ! i_y
            End Do    ! i_x

!           Now when all Cartesian components have been computed we
!           transform back to the coordinate system of the molecule.

!           ... more to come ...

         End Do ! iAlpha
      End Do    ! iBeta
!
      Return

! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nRys)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If

      End Subroutine twlint
!
      Subroutine twlprm(Zeta,P_x,P_y,Alpha,Beta,i_x,i_y,j_x,j_y,l,intg)
      Implicit None
!
!
!     Local Arrays and integers
!
  Integer :: n, m, k, q, i_x, i_y, j_x, j_y, p, t
  Integer :: a, b, v, w, c, d, g, f, l
  Real*8 :: Lambda, Theta, Omega, Gamma, Zeta, intg
  Real*8 :: a1, a2, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2
  Real*8 :: g1, g2, h1, h2, k1, k2, l1, l2, l3
  Real*8 :: A_x, A_y, B_x, B_y, Alpha, Beta, P_x, P_y  ! parameters
  Real*8, parameter :: pi = atan(1.0)*4.d0
  Complex*16 cx, Kappa, i

  i=DCmplx(0.0D0,1.0D0)

  intg = 0.0D0

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
           c2 = (P_y - A_y)**n
           Do q = 0, j_y

              d1 = choose(j_y, q)
              d2 = (P_y - B_y)**n

              b = i_y + j_y - k - q

              Lambda = (a1*a2*b1*b2*c1*c2*d1*d2* &
                   exp((-(Alpha*Beta)/Zeta)*(A_x-B_x)**2+(A_y - B_y)**2))

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
                          h2 = (1/(2*i))**(b-w)

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
