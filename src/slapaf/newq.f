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
* Copyright (C) 1994, Roland Lindh                                     *
************************************************************************
      SubRoutine Newq(q,nInter,nIter,dq,H,g,error,B,RHS,
     &                Scrt1,nScrt1,
     &                Beta,nFix,iP,Energy,
     &                Step_Trunc,Thr_RS)
************************************************************************
*                                                                      *
* Object: Driver for optimization procedures.                          *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             December '94                                             *
************************************************************************
      use Slapaf_Parameters, only: iOptC, UpMeth, Line_Search, E_Delta
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
#include "stdalloc.fh"
C     Logical Print
      Real*8 q(nInter,nIter+1), dq(nInter,nIter), g(nInter,nIter),
     &       error(nInter,nIter+1), B((nIter+1)*(nIter+1)),
     &       RHS(nIter+1),
     &       Scrt1(nScrt1), Energy(nIter), H(nInter,nInter)
      Integer   iP(nIter)
      Character*1 Step_Trunc
      Real*8, Allocatable:: t_q(:), t_g(:), t_dq(:)
*
      Lu=6
      iRout = 113
      iPrint = nPrint(iRout)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (Lu,*) ' Newq: nIter,Beta=',nIter,Beta
      Call RecPrt(' Newq: q',' ',q,nInter,nIter+1)
      Call RecPrt(' Newq: dq',' ',dq,nInter,nIter)
      Call RecPrt(' Newq: g',' ',g,nInter,nIter)
      Call RecPrt(' Newq: H   ',' ',H   ,nInter,nInter)
#endif
*
*---- Print out of the Hessian and determination of the Hessian index.
*
C     print=.true.
C     Call View(H,nInter,print)
*
      E_Delta=Zero
*
*---- Perform first a linear search for the last two points
*     to find minimum along the direction of the last step.
*
*     The new point will temporarily replace the last point!
*
      If (iPrint.ge.6) Write (Lu,*)
      If (Line_Search) Then
         If (nIter.ge.2) Then
            Call mma_allocate(t_q ,nInter,Label='t_q')
            Call mma_allocate(t_g ,nInter,Label='t_g')
            Call mma_allocate(t_dq,nInter,Label='t_dq')
*
            call dcopy_(nInter,dq(1,nIter-1),1,t_dq,1)
            call dcopy_(nInter, q(1,nIter  ),1,t_q ,1)
            call dcopy_(nInter, g(1,nIter  ),1,t_g ,1)
*
            Call LnSrch(Energy,q,dq,g,nInter,nIter,E_Delta)
         Else
            If (iPrint.ge.6) Write (Lu,*)
     &          '-- First iteration no line search'
         End If
      Else
         If (iPrint.ge.6) Write (Lu,*) '-- Line search is disabled'
      End If
      If (iPrint.ge.6) Write (Lu,*)
*
*-----Invoke the quadratic optimization procedure
*
*     iOptC=0001   : quasi Newton-Raphson
*     iOptC=0010   : C1-DIIS
*     iOptC=0100   : C2-DIIS
*     iOptC=1000   : Rational Function
*
      call dcopy_(nInter,[Zero],0,Scrt1,1)
*                                                                      *
************************************************************************
*                                                                      *
      If (iOptC.eq.0) Then
*
*------- No update of the geometry
*
         Call FZero(dq,nInter)
*                                                                      *
************************************************************************
*                                                                      *
      Else If (iAnd(iOptC,1).eq.1) Then
*
*------- quasi Newton-Raphson
*
         UpMeth=' qNR  '
         Call QNR(nInter,nIter,dq,H,g)
         Beta_New=Sqrt(DDot_(nInter,dq(1,nIter),1,dq(1,nIter),1))
         If (Beta_New.gt.Beta) Then
            Call DScal_(nInter,Beta/Beta_New,dq(1,nIter),1)
            Step_Trunc='*'
         End If
*
*                                                                      *
************************************************************************
*                                                                      *
      Else If (iAnd(iOptC,2).eq.2) Then
*
*------- C1-DIIS
*
         UpMeth='c1DIIS'
         MinWdw=5
         Call C1DIIS(q,nInter,nIter,dq,H,g,error,B,RHS,nFix,iP,MinWdw)
         Beta_New=Sqrt(DDot_(nInter,dq(1,nIter),1,dq(1,nIter),1))
         If (Beta_New.gt.Beta) Then
            Call DScal_(nInter,Beta/Beta_New,dq(1,nIter),1)
            Step_Trunc='*'
         End If
*
*                                                                      *
************************************************************************
*                                                                      *
      Else If (iAnd(iOptC,4).eq.4) Then
*
*------- C2-DIIS
*
         UpMeth='c2DIIS'
         Call C2DIIS(q,nInter,nIter,dq,H,g,error,B,RHS,
     &               Scrt1,nScrt1,nFix,iP)
         Beta_New=Sqrt(DDot_(nInter,dq(1,nIter),1,dq(1,nIter),1))
         If (Beta_New.gt.Beta) Then
            Call DScal_(nInter,Beta/Beta_New,dq(1,nIter),1)
            Step_Trunc='*'
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else If (iAnd(iOptC,8).eq.8) Then
*
*------- Rational Function Optimization (RFO)
*
         UpMeth='  RF  '
*
         If (iAnd(iOptC,128).ne.128) Then
*
*---------- TS optimization
*
            If (iAnd(iOptC,512).ne.512) Then
*
*------------- Restricted Step Partitioned RFO
*
               Call RS_P_RFO(H,g(1,nIter),nInter,dq(1,nIter),
     &                       UpMeth,E_Delta,Beta,Step_Trunc)
            Else
*
*------------- Restricted Step Image RFO
*
               Call RS_I_RFO(H,g(1,nIter),nInter,dq(1,nIter),
     &                       UpMeth,E_Delta,Beta,Step_Trunc,Thr_RS)
            End If
*
         Else
*
*---------- Restricted Step RFO
*
*
            Call RS_RFO(H,g(1,nIter),nInter,dq(1,nIter),
     &                  UpMeth,E_Delta,Beta,Step_Trunc,Thr_RS)
*
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Call WarningMessage(2,'Error in NewQ')
         Write (6,*) ' Newq: Illegal setting of iOptC!'
         Write (6,*) '  iOptC=',iOptC
         Call Abend
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- In case of a line search restore some data and add the
*     replacements.
*
*
      If (Line_Search.and.nIter.ge.2) Then
         If (iPrint.ge.99) Then
            Call RecPrt(' Newq: q ',' ',q, nInter,nIter+1)
            Call RecPrt(' Newq: dq',' ',dq,nInter,nIter  )
            Call RecPrt(' Newq: g ',' ',g, nInter,nIter)
         End If
         call dcopy_(nInter,q(1,nIter),1,q(1,nIter+1),1)
         Call DaXpY_(nInter,One,dq(1,nIter),1,q(1,nIter+1),1)
         call dcopy_(nInter,t_q ,1, q(1,nIter  ),1)
         call dcopy_(nInter,q(1,nIter+1),1,dq(1,nIter),1)
         Call DaXpY_(nInter,-One,q(1,nIter),1,dq(1,nIter),1)
*
         call dcopy_(nInter,t_dq, 1,dq(1,nIter-1),1)
         call dcopy_(nInter,t_g , 1, g(1,nIter  ),1)
         If (iPrint.ge.99) Then
            Call RecPrt(' Newq: q ',' ',q, nInter,nIter+1)
            Call RecPrt(' Newq: dq',' ',dq,nInter,nIter  )
            Call RecPrt(' Newq: g ',' ',g, nInter,nIter)
         End If
*
         Call mma_deallocate(t_q)
         Call mma_deallocate(t_g)
         Call mma_deallocate(t_dq)
      End If
*
*-----Estimate energy at the relaxed geometry
*
      If (iAnd(iOptC,8).ne.8) Then
*
*..... g(r)
*
         call dcopy_(nInter,g(1,nIter),1,Scrt1,1)
         Call DScal_(nInter,-One,Scrt1,1)
*
*.....1/2 H(r) (r -r)
*                0
*
         Call dGeMV_('N',nInter,nInter,
     &              Half,H,nInter,dq(1,nIter),1,
     &              One,Scrt1,1)
         E_Delta=DDot_(nInter,Scrt1,1,dq(1,nIter),1)
      End If
*
      Do i = 1, nInter
        q(i,nIter+1) = q(i,nIter) + dq(i,nIter)
      End Do
*
#ifdef _DEBUGPRINT_
      Write (Lu,*) ' E_Delta=',E_Delta
      Call RecPrt('Newq: q',' ',q,nInter,nIter+1)
      Call RecPrt('Newq: dq',' ',dq,nInter,nIter)
      Call RecPrt('Newq: g',' ',g,nInter,nIter)
#endif
      Return
      End
