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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine P_Int (Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,
     &                  Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,
     &                  Array,nArr,Ccoor,nOrdOp,lOper,iChO,
     &                  iStabM,nStabM,
     &                  PtChrg,nGrid,iAddPot)
************************************************************************
*                                                                      *
* Object: to compute the multipole moments integrals with the          *
*         Gauss-Hermite quadrature.                                    *
*                                                                      *
* Called from: OneEl                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              RecPrt                                                  *
*              CrtCmp                                                  *
*              SOS                                                     *
*              DCR                                                     *
*              Assmbl                                                  *
*              GetMem                                                  *
*              DCopy   (ESSL)                                          *
*              CmbnMP                                                  *
*              SymAdO                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
*             Modified to multipole moments November '90               *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "oneswi.fh"
#include "print.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),
     &       Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),
     &       rKappa(nZeta), P(nZeta,3), A(3), RB(3),
     &       Array(nZeta*nArr), Ccoor(3)
      Character*80 Label
      Integer lOper(nComp), iStabM(0:nStabM-1), iChO(nComp)
*
*     Statement function for Cartesian index
*
      nElem(i) = (i+1)*(i+2)/2
*
      iRout = 122
      iPrint = nPrint(iRout)
*     Call qEnter('P_Int')
*
*---- Observe that this code does not make any sense in case of symmetry!
*
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
*
      If (iPrint.ge.99) Then
         Write (6,*) ' Result in P_Int'
         Do ia = 1, nElem(la)
            Do ib = 1, nElem(lb)
               Do iIC = 1, nIC
                  Write (Label,'(A,I2,A,I2,A,I2,A)')
     &               ' Final(a=',ia,',b=',ib,',iIC=',iIC,')'
                  Call RecPrt(Label,' ',Final(1,ia,ib,iIC),nAlpha,nBeta)
               End Do
            End Do
         End Do
      End If
*
*     Call GetMem(' Exit P_Int','LIST','REAL',iDum,iDum)
*     Call qExit('p_Int')
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Alpha)
         Call Unused_real_array(Beta)
         Call Unused_real_array(Zeta)
         Call Unused_real_array(ZInv)
         Call Unused_real_array(rKappa)
         Call Unused_real_array(P)
         Call Unused_real_array(A)
         Call Unused_real_array(RB)
         Call Unused_integer(nHer)
         Call Unused_real_array(Array)
         Call Unused_real_array(Ccoor)
         Call Unused_integer(nOrdOp)
         Call Unused_integer_array(lOper)
         Call Unused_integer_array(iChO)
         Call Unused_integer_array(iStabM)
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
      End
