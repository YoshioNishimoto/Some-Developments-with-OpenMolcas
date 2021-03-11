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
* Copyright (C) 1991,1993,1998, Roland Lindh                           *
************************************************************************
      SubRoutine CmpctR(abcd,na,nb,nijkl,mijkl,Zeta,Kappab,P,
     &                  Ind_Pair,Con,xZeta,xKapp,xP,IndZ,iOff,Jnd,
     &                  xZInv,CutInt,RadMax,cdMax,EtMax,AeqB,xab,
     &                  xabCon,Alpha,xAlpha,Beta,xBeta)
************************************************************************
*                                                                      *
* Object: to find the largest value of the integrals which will be     *
*         generated by an exponent combination. This information will  *
*         be used to prescreen the data.                               *
*                                                                      *
* Called from: k2Loop                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             June '91                                                 *
*             Modified for direct SCF, January '93                     *
*                                                                      *
*             Modified for k2 prescreening, June '98                   *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 abcd(mijkl,na,nb,na,nb), Zeta(mijkl), xab(nijkl),
     &       KappAB(mijkl), P(nijkl,3), xZeta(nijkl), xKapp(nijkl),
     &       xP(nijkl,3), xZInv(nijkl), Con(mijkl), xabCon(nijkl),
     &       Alpha(mijkl), xAlpha(nijkl), Beta(mijkl), xBeta(nijkl)
      Integer IndZ(nijkl+1), Ind_Pair(mijkl)
      Logical AeqB
*
      iRout = 245
      iPrint = nPrint(iRout)
*     Call QEnter('Cmpct')
      If (iPrint.ge.59) Then
         Write (6,*) ' In CmpctS'
         Write (6,*) AeqB,iOff,Jnd
         Call RecPrt('Zeta',' ',Zeta,mijkl,1)
         Call RecPrt('abcd',' ',abcd,mijkl,(na*nb)**2)
      End If
*
*     Move data
*
      If (AeqB) Then
         Call ICopy(mijkl,Ind_Pair,1,IndZ(iOff+1),  1)
         call dcopy_(mijkl,Zeta,    1,xZeta(iOff+1), 1)
         call dcopy_(mijkl,KappAB,  1,xKapp(iOff+1), 1)
         call dcopy_(mijkl,P(1,1),  1,xP(iOff+1,1),  1)
         call dcopy_(mijkl,P(1,2),  1,xP(iOff+1,2),  1)
         call dcopy_(mijkl,P(1,3),  1,xP(iOff+1,3),  1)
         call dcopy_(mijkl,Alpha,   1,xAlpha(iOff+1),1)
         call dcopy_(mijkl,Beta,    1,xBeta(iOff+1), 1)
         Do ijkl = 1, mijkl
            ijkl_=Ind_Pair(ijkl)
            xZInv(iOff+ijkl)= One/Zeta(ijkl)
            Temp = Zero
            Do ia = 1, na
               Do ib = 1, nb
                  Temp=Max(Temp,Abs(abcd(ijkl,ia,ib,ia,ib)))
               End Do
            End Do
            xab   (ijkl+iOff)= Sqrt(Temp)
            xabCon(ijkl+iOff)= Sqrt(Temp)*Con(ijkl_)
         End Do
         Jnd = Jnd + mijkl
      Else
         Do ijkl = 1, mijkl
            ijkl_=Ind_Pair(ijkl)
*
*-----------Select the largest element (over components)
            Temp = Zero
            Do ia = 1, na
               Do ib = 1, nb
                  Temp=Max(Temp,Abs(abcd(ijkl,ia,ib,ia,ib)))
               End Do
            End Do
            Temp1= Sqrt(Temp)
            Temp2= Sqrt(Temp)* Con(ijkl_)
            If (KappAB(ijkl)*Con(ijkl)*RadMax.ge.CutInt) Then
               Jnd = Jnd + 1
               xZeta(Jnd) =Zeta(ijkl)
               xKapp(Jnd) =KappAB(ijkl)
               xP(Jnd,1)  =P(ijkl,1)
               xP(Jnd,2)  =P(ijkl,2)
               xP(Jnd,3)  =P(ijkl,3)
               xZInv(Jnd) =One/Zeta(ijkl)
               IndZ(Jnd)  =Ind_Pair(ijkl)
               xab   (Jnd)=Temp1
               xabCon(Jnd)=Temp2
               xAlpha(Jnd)=Alpha(ijkl)
               xBeta(Jnd)=Beta(ijkl)
            End If
         End Do
      End If
      IndZ(nijkl+1)=Jnd
*
      If (iPrint.ge.99) Then
         Write (6,*) 'AeqB=',AeqB
         Write (6,*) 'IndZ=',IndZ
         Call RecPrt('xZeta ',' ',xZeta,  1,nijkl)
         Call RecPrt('xKapp ',' ',xKapp,  1,nijkl)
         Call RecPrt('xP(x) ',' ',xP(1,1),1,nijkl)
         Call RecPrt('xP(y) ',' ',xP(1,2),1,nijkl)
         Call RecPrt('xP(z) ',' ',xP(1,3),1,nijkl)
         Call RecPrt('xZInv ',' ',xZInv,  1,nijkl)
         Call RecPrt('xab   ',' ',xab,    1,nijkl)
         Call RecPrt('xabCon',' ',xabCon, 1,nijkl)
         Call RecPrt('xAlpha',' ',Alpha,  1,nijkl)
         Call RecPrt('xBeta ',' ',Beta,   1,nijkl)
      End If
*
*     Call GetMem(' Exit Cmpct','Check','Real',iDum,iDum)
*     Call QExit('Cmpct')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(cdMax)
         Call Unused_real(EtMax)
      End If
      End
