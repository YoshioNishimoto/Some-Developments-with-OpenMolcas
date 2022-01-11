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
* Copyright (C) 2010, Yan Zhao                                         *
************************************************************************
      Subroutine XM06_L(mGrid,CoeffA,iSpin,F_xc)
      Implicit None
      Integer mGrid, iSpin
      Real*8 CoeffA, F_xc(mGrid)
      Integer, parameter :: ijzy=1
      Call XM06(mGrid,CoeffA,iSpin,F_xc,ijzy)
      End Subroutine XM06_L

      Subroutine XM06_HF(mGrid,CoeffA,iSpin,F_xc)
      Implicit None
      Integer mGrid, iSpin
      Real*8 CoeffA, F_xc(mGrid)
      Integer, parameter :: ijzy=2
      Call XM06(mGrid,CoeffA,iSpin,F_xc,ijzy)
      End Subroutine XM06_HF

      Subroutine XM06_2X(mGrid,CoeffA,iSpin,F_xc)
      Implicit None
      Integer mGrid, iSpin
      Real*8 CoeffA, F_xc(mGrid)
      Integer, parameter :: ijzy=4
      Call XM06(mGrid,CoeffA,iSpin,F_xc,ijzy)
      End Subroutine XM06_2X

      Subroutine XM06(mGrid,CoeffA,iSpin,F_xc,ijzy)
************************************************************************
*                                                                      *
*  M06x evaluates the exchange part of the M06 suite of                *
*  functionals on a grid.                                              *
*  !!! Second derivatives are not available yet.                       *
*                                                                      *
*  Ref: (a) Zhao, Y.  and Truhlar, D. G. J. Chem. Phys. 125,           *
*    194101 (2006).                                                    *
*       (b) Y. Zhao and D. G. Truhlar, J. Phys. Chem. A (2006),        *
*    110(49),13126-13130.                                              *
*                                                                      *
*       ijzy - 1 M06-L                                                 *
*       ijzy - 2 M06-HF                                                *
*       ijzy - 3 M06                                                   *
*       ijzy - 4 M06-2X                                                *
*                                                                      *
*  YZ (10/07)                                                          *
************************************************************************
      use nq_Grid, only: Rho, Sigma, Tau
      use nq_Grid, only: vRho, vSigma, vTau
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20
      Integer mGrid

      integer ijzy
      REAL*8 at1, at2, at3, at4, at5, at6, at7, at8, at9
      REAL*8 at, at10, at11, C1, C2, fL, fNL, at0
      REAL*8 rrho, rho43, rho13, rhoo, rho53
      REAL*8 Gamma
      REAL*8 TauUEG, Tsig, Wsig, Fsig
      REAL*8 tauu
      REAL*8 F5o3, F1o3
      REAL*8 F2o3, F3o2, F4o3, F3o5
      REAL*8 F11
      REAL*8  Ax,  x, x2, En, Ed, E, dE, dEn, dEd
      REAL*8 dFdW, dWdT, dTdR, dTdTau, dGGAdR, dFdR
      REAL*8 dFdTau, dGGAdG
      parameter (F1o3=1.d0/3.d0, F2o3=2.d0/3.d0, F3o2=3.d0/2.d0)
      parameter (F4o3=4.d0/3.d0, F3o5=3.d0/5.d0, F5o3=5.d0/3.d0)
      parameter (F11=11.d0)

      if (ijzy.eq.1) then
C     Parameters for M06-L
        at0=    3.987756D-01
        at1=    2.548219D-01
        at2=    3.923994D-01
        at3=    -2.103655D+00
        at4=    -6.302147D+00
        at5=    1.097615D+01
        at6=    3.097273D+01
        at7=    -2.318489D+01
        at8=    -5.673480D+01
        at9=    2.160364D+01
        at10=   3.421814D+01
        at11=   -9.049762D+00
       elseif (ijzy.eq.2) then
C     Parameters for M06-HF
        at0=    1.179732D-01
        at1=    -1.066708D+00
        at2=    -1.462405D-01
        at3=    7.481848D+00
        at4=    3.776679D+00
        at5=    -4.436118D+01
        at6=    -1.830962D+01
        at7=    1.003903D+02
        at8=    3.864360D+01
        at9=    -9.806018D+01
        at10=   -2.557716D+01
        at11=   3.590404D+01
       elseif (ijzy.eq.3) then
C     Parameters for M06
        at0=    5.877943D-01
        at1=    -1.371776D-01
        at2=    2.682367D-01
        at3=    -2.515898D+00
        at4=    -2.978892D+00
        at5=    8.710679D+00
        at6=    1.688195D+01
        at7=    -4.489724D+00
        at8=    -3.299983D+01
        at9=    -1.449050D+01
        at10=   2.043747D+01
        at11=   1.256504D+01
C      elseif (ijzy.eq.4) then
       else
C     Parameters for M06-2X
        at0=    4.600000D-01
        at1=    -2.206052D-01
        at2=    -9.431788D-02
        at3=    2.164494D+00
        at4=    -2.556466D+00
        at5=    -1.422133D+01
        at6=    1.555044D+01
        at7=    3.598078D+01
        at8=    -2.722754D+01
        at9=    -3.924093D+01
        at10=   1.522808D+01
        at11=   1.522227D+01
      endif
*
      at=1.0d0
      C1     = 3.36116D-03
      C2     = 4.49267D-03

      fL  =  1.0d0
      fNL =  1.0d0
*
      Ax = -F3o2*(F4o3*PI)**(-F1o3)
      Ta=0.5D0*T_X
*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
        Do iGrid = 1, mGrid
         if(Rho(1,iGrid).lt.Ta .or.  Tau(1,iGrid).lt.Ta) goto 110
         rhoo=max(1.0D-24,rho(1,igrid),Ta)
         rho43 = rhoo**F4o3
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F5o3
*
         tauu=Tau(1,iGrid)
         TauUEG=F3o5*((Six*pi*pi)**F2o3)*rho53
         Tsig =TauUEG/tauu
         Wsig =(Tsig-One)/(Tsig+One)
         Fsig=at*(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))
         Gamma = Sqrt(Sigma(1,iGrid))
         x = Gamma/rho43
         x2 = x*x
         En = C1*x2
         Ed = One + C2*x2
         E  = -En/Ed
         F_xc(igrid)=F_xc(igrid)+(Ax*fL+fNL*E)*Fsig*rho43*2.0D0
*
*     functional derivatives
*
         dEn   = Two*C1*x
         dEd   = Two*C2*x
         dE    = -(dEn*Ed-En*dEd)/(Ed*Ed)
         dFdW=at*(      at1 + Wsig*(Two  *at2 + Wsig*(Three*at3 + Wsig*(
     &            Four *at4 + Wsig*(Five *at5 + Wsig*(Six  *at6 + Wsig*(
     &            Seven*at7 + Wsig*(Eight*at8 + Wsig*(Nine *at9 + Wsig*(
     &            Ten  *at10+ Wsig*F11*at11))))))))))
         dWdT = Two/((One + Tsig)**2)
         dTdR = ((Six*PI*PI)**F2o3)*(rhoo**F2o3)/tauu
         dTdTau = -TauUEG/tauu**2
         dGGAdR = F4o3*rho13*(fL*Ax+fNL*(E-x*dE))
         dFdR = dFdW*dWdT*dTdR
         dFdTau=dFdW*dWdT*dTdTau
         dGGAdG =(fNL*dE/(Two*Gamma))
*        dF/dRho
         vRho(1,igrid)=vRho(1,igrid)+dGGAdR*Fsig
     &        + (fL*Ax+fNL*E)*rho43*dFdR
*        dF/dGamma
         vSigma(1,igrid)=vSigma(1,iGrid)+ dGGAdG*Fsig
*        dF/dTau
         vTau(1,iGrid)=vTau(1,iGrid)+rho43*(Ax*fL+fNL*E)*dFdTau

110     continue
       Enddo
*                                                                      *
************************************************************************
*                                                                      *
*     ispin .ne. 1, use both alpha and beta components.
*
      Else
*                                                                      *
************************************************************************
*                                                                      *
        Do iGrid = 1, mGrid
*
* alpha component
*
         If(Rho(1,iGrid).lt.Ta .or. Tau(1,iGrid).lt.Ta) goto 210
         rhoo=Rho(1,iGrid)
         rho43 = rhoo**F4o3
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F5o3
*
         tauu=Tau(1,iGrid)
         TauUEG=F3o5*((Six*pi*pi)**F2o3)*rho53
         Tsig =TauUEG/tauu
         Wsig =(TauUEG - tauu)/(TauUEG + tauu)
         Fsig=at*(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))
         Gamma = Sqrt(Sigma(1,iGrid))
         x = Gamma/rho43
         x2 = x*x
         En = C1*x2
         Ed = One + C2*x2
         E  = -En/Ed
         F_xc(iGrid)=F_xc(iGrid)+(Ax*fL+fNL*E)*Fsig*rho43
*
*     functional derivatives
*

         dEn   = Two*C1*x
         dEd   = Two*C2*x
         dE    = -(dEn*Ed-En*dEd)/(Ed*Ed)
         dFdW=at*(      at1 + Wsig*(Two  *at2 + Wsig*(Three*at3 + Wsig*(
     &            Four *at4 + Wsig*(Five *at5 + Wsig*(Six  *at6 + Wsig*(
     &            Seven*at7 + Wsig*(Eight*at8 + Wsig*(Nine *at9 + Wsig*(
     &            Ten  *at10+ Wsig*F11*at11))))))))))
         dWdT = Two/((One + Tsig)**2)
         dTdR = ((Six*PI*PI)**F2o3)*(rhoo**F2o3)/tauu
         dTdTau = -TauUEG/tauu**2
         dGGAdR = F4o3*rho13*(fL*Ax+fNL*(E-x*dE))
         dFdR = dFdW*dWdT*dTdR
         dFdTau=dFdW*dWdT*dTdTau
         dGGAdG =(fNL*dE/(Two*Gamma))
*        dF/dRhoa
         vRho(1,iGrid)=vRho(1,igrid)+dGGAdR*Fsig
     &        + (fL*Ax+fNL*E)*rho43*dFdR
*        dF/dGammaaa
         vSigma(1,iGrid)=vSigma(1,iGrid)+ dGGAdG*Fsig
*        dF/dTaua
         vTau(1,iGrid)=vTau(1,iGrid)+rho43*(Ax*fL+fNL*E)*dFdTau
210      continue
*
* beta component
*
         if(Rho(2,iGrid).lt.Ta .or. Tau(2,iGrid).lt.Ta) goto 310
         rhoo=Rho(2,iGrid)
         rho43 = rhoo**F4o3
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F5o3
*
         tauu=Tau(2,iGrid)
         TauUEG=F3o5*((Six*pi*pi)**F2o3)*rho53
         Tsig =TauUEG/tauu
         Wsig =(Tsig-One)/(Tsig+One)
         Fsig=at*(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))
         Gamma = sqrt(Sigma(3,iGrid))
         x = Gamma/rho43
         x2 = x*x
         En = C1*x2
         Ed = One + C2*x2
         E  = -En/Ed
         F_xc(iGrid)=F_xc(iGrid)+(Ax*fL+fNL*E)*Fsig*rho43
*
*     functional derivatives
*
         dEn   = Two*C1*x
         dEd   = Two*C2*x
         dE    = -(dEn*Ed-En*dEd)/(Ed*Ed)
         dFdW=at*(      at1 + Wsig*(Two  *at2 + Wsig*(Three*at3 + Wsig*(
     &            Four *at4 + Wsig*(Five *at5 + Wsig*(Six  *at6 + Wsig*(
     &            Seven*at7 + Wsig*(Eight*at8 + Wsig*(Nine *at9 + Wsig*(
     &            Ten  *at10+ Wsig*F11*at11))))))))))
         dWdT = Two/((One + Tsig)**2)
         dTdR = ((Six*PI*PI)**F2o3)*(rhoo**F2o3)/tauu
         dTdTau = -TauUEG/tauu**2
         dGGAdR = F4o3*rho13*(fL*Ax+fNL*(E-x*dE))
         dFdR = dFdW*dWdT*dTdR
         dFdTau=dFdW*dWdT*dTdTau
         dGGAdG =(fNL*dE/(Two*Gamma))
*        dF/dRhob
         vRho(2,igrid)=vRho(2,igrid)+dGGAdR*Fsig
     &        + (fL*Ax+fNL*E)*rho43*dFdR
*        dF/dGammabb
         vSigma(3,igrid)=vSigma(3,iGrid)+ dGGAdG*Fsig
*        dF/dTaub
         vTau(2,igrid)=vTau(2,iGrid)+rho43*(Ax*fL+fNL*E)*dFdTau
310     continue
       enddo
*                                                                      *
************************************************************************
*                                                                      *
       End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(CoeffA)
      End
