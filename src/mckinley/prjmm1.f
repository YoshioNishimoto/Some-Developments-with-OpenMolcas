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
* Copyright (C) 1993, Roland Lindh                                     *
************************************************************************
      Subroutine PrjMm1(nHer,MmPrjG,la,lb,lr)
************************************************************************
*                                                                      *
*  Object: to compute the number of real*8 the kernal routine will     *
*          need for the computation of a matrix element between two    *
*          cartesian Gaussin functions with the total angular momentum *
*          of la and lb (la=0 s-function, la=1 p-function, etc.)       *
*          lr is the order of the operator (this is only used when the *
*          integrals are computed with the Hermite-Gauss quadrature).  *
*                                                                      *
*  Called from: OneEl                                                  *
*                                                                      *
************************************************************************
*
#include "itmax.fh"
#include "info.fh"
*
      nElem(i) = (i+1)*(i+2)/2
*
      nOrder = 0
      ld=1
      MmPrjG = 0
      Do 1960 iCnttp = 1, nCnttp
         If (.Not.ECP(iCnttp)) Go To 1960
         Do 1966 iAng = 0, nPrj_Shells(iCnttp)-1
            iShll = ipPrj(iCnttp) + iAng
            If (nExp(iShll).eq.0 .or. nBasis(iShll).eq.0) Go To 1966
*
            ip = 0

            nac = nElem(la)*nElem(iAng)
            ncb = nElem(iAng)*nElem(lb)

            ip=ip+6*nelem(la)*nelem(lb) ! final
            ip=ip+4*nac*nExp(ishll) ! FA1
            ip=ip+4*ncb*nExp(ishll) !FB1
            ip=ip+nexp(ishll)* nExp(ishll) !Tmp

            nHer = (la+1+iAng+1+ld)/2
            nOrder = Max(nHer,nOrder)
            iacore=6+3*nHer*(la+1+ld)+3*nHer*(iAng+1)+
     &           3*nHer*(lr+1)+3*(la+1+ld)*(iAng+1)*(lr+1)+1

            nHer = (lb+1+iAng+1+ld)/2
            nOrder = Max(nHer,nOrder)
            icoreb=6+3*nHer*(lb+1+ld)+3*nHer*(iAng+1)+
     &           3*nHer*(lr+1)+3*(lb+1+ld)*(iAng+1)*(lr+1)+1

            icores = MAX(icoreb,iacore)*nExp(ishll)

            MmPrjG = Max(MmPrjG,ip+icores)
*
 1966    Continue
 1960 Continue
      nHer = nOrder
c      !write(*,*) 'mem',MmPrjG
*
      Return
      End
