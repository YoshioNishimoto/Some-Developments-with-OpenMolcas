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
* Copyright (C) 1990,1992, Roland Lindh                                *
*               1990, IBM                                              *
************************************************************************
      SubRoutine TNAI1(Zeta,Eta,P,Q,nT,T,ZEInv,IsChi,ChiI2)
************************************************************************
*                                                                      *
* Object: to entities for the nucelar attraction integrals which are   *
*         used in the Rys quadrature to evaluate these integrals.      *
*                                                                      *
* Called from: Rys                                                     *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             March '90                                                *
*                                                                      *
*             March '92 modified to gradient calculation.              *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
      Real*8 Zeta(nT), Eta(nT), P(nT,3), Q(nT,3),
     &       ZEInv(nT), T(nT)
*
      iRout = 57
      iPrint = nPrint(iRout)
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt(' Zeta in TNAI1',' ',Zeta,nT,1)
         Call RecPrt(' Eta in TNAI1',' ',Eta,nT,1)
         Call RecPrt(' P in TNAI1',' ',P,nT,3)
         Call RecPrt(' Q in TNAI1',' ',Q,nT,3)
      End If
#endif
      Do iT = 1, nT
         PQ2 = (P(iT,1)-Q(iT,1))**2
     &       + (P(iT,2)-Q(iT,2))**2
     &       + (P(iT,3)-Q(iT,3))**2
         T(iT) = Zeta(iT)*PQ2
         ZEInv(iT) = 1.0D0/Zeta(iT)
      End Do
*
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt('Tvalue',' ',T,nT,1)
      End If
#endif
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Eta)
      If (.False.) Call Unused_integer(IsChi)
      If (.False.) Call Unused_real(ChiI2)
      Return
      End
