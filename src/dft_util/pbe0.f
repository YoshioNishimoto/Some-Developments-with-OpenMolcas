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
* Copyright (C) 2005, Per Ake Malmqvist                                *
************************************************************************
      Subroutine PBE0(mGrid,Rho,nRho,P2_ontop,
     &                nP2_ontop,iSpin,F_xc,
     &                dF_dRho,ndF_dRho,dF_dP2ontop,ndF_dP2ontop,T_X)
************************************************************************
*                                                                      *
* Object: To compute the PBE0 density functional and its first         *
* derivatives. This functional is defined as a hybrid functional       *
* having 25% HF exchange, 75% PBE exchange, and 100% PBE correlation.  *
*   See   J.P. Perdew, M. Ernzerhof, and K. Burke,                     *
*      J. Chem. Phys. 105, 9982 (1996).                                *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. December 2005                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       P2_ontop(nP2_ontop,mGrid), F_xc(mGrid),
     &       dF_dP2ontop(ndF_dP2ontop,mGrid)
*                                                                      *
************************************************************************
*                                                                      *
      CoeffA=1.0D0*CoefR
      Call CPBE(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &          CoeffA,iSpin,F_xc,T_X)

      CoeffB=0.75D0*CoefX
      Call XPBE(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &          CoeffB,iSpin,F_xc,T_X)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(P2_ontop)
         Call Unused_real_array(dF_dP2ontop)
      End If
      End
