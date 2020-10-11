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
*               2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine PLF_LDF_JK_2P_2(TInt,nTInt,Map,
     &                           AOint,ijkl,iCmp,jCmp,kCmp,lCmp,
     &                           iAO,iAOst,iBas,jBas,kBas,lBas,kOp)
************************************************************************
*                                                                      *
*  object: to sift and index the petite list format integrals.         *
*                                                                      *
*          the indices has been scrambled before calling this routine. *
*          Hence we must take special care in order to regain the can- *
*          onical order.                                               *
*                                                                      *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, Ca     *
*          May '90                                                     *
*          Modified for Local DF, Thomas Bondo Pedersen, June 2010     *
*                                                                      *
************************************************************************
      use SOAO_Info, only: iAOtSO
      Implicit None
      Integer nTInt
      Real*8  TInt(nTInt)
      Integer Map(4)
      Integer ijkl, iCmp, jCmp, kCmp, lCmp
      Real*8  AOInt(ijkl,iCmp,jCmp,kCmp,lCmp)
      Integer iAO(4), iAOst(4)
      Integer iBas, jBas, kBas, lBas
      Integer kOp(4)
#include "localdf_bas.fh"
#include "localdf_int2.fh"
#include "WrkSpc.fh"

      Character*15 SecNam
      Parameter (SecNam='PLF_LDF_JK_2P_2')

      Integer i1, i2, i3, i4
      Integer iShlJ, iShlK, iShlL
      Integer jSO, lSO, kSO
      Integer jSOj, lSOl, kSOk
      Integer ll, jj, kk
      Integer iShlKL, kl, kl0
      Integer nijkl
      Integer Jkl

      Integer i, j
      Integer nBasSh, iShlSO, iRow, iCol
#if defined (_DEBUGPRINT_)
      Integer iSOShl
      iSOShl(i)=iWork(ip_iSOShl-1+i)
#endif
      nBasSh(i)=iWork(ip_nBasSh-1+i)
      iShlSO(i)=iWork(ip_iShlSO-1+i)
      iRow(i,j)=iWork(ip_AB_IndxG-1+l_AB_IndxG_1*(j-1)+i)
      iCol(i,j)=iWork(ip_CD_IndxG2-1+l_CD_IndxG2_1*(j-1)+i)

      If (Map(1).eq.1 .and. Map(2).eq.2 .and.
     &    Map(3).eq.3 .and. Map(4).eq.4) Then
#if defined (_DEBUGPRINT_)
         If (iSOShl(iAOtSO(iAO(1)+1,kOp(1))+iAOst(1)).ne.SHA) Then
            Call WarningMessage(2,SecNam//': Shell problem [1]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(2)+1,kOp(2))+iAOst(2)).ne.SHB) Then
            Call WarningMessage(2,SecNam//': Shell problem [2]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(3)+1,kOp(3))+iAOst(3)).ne.SHC) Then
            Call WarningMessage(2,SecNam//': Shell problem [3]')
            Call LDF_Quit(1)
         End If
         If (iSOShl(iAOtSO(iAO(4)+1,kOp(4))+iAOst(4)).ne.SHD) Then
            Call WarningMessage(2,SecNam//': Shell problem [4]')
            Call LDF_Quit(1)
         End If
         If (iCmp*iBas.ne.1) Then
            Call WarningMessage(2,
     &                          SecNam//': Shell dimension problem [1]')
            Call LDF_Quit(1)
         End If
#endif
         iShlJ=SHB
         iShlK=SHC
         iShlL=SHD
         iShlKL=SPCD
         i1=1
         Do i4=1,lCmp
            lSO=iAOtSO(iAO(4)+i4,kOp(4))+iAOst(4)
            Do i3=1,kCmp
               kSO=iAOtSO(iAO(3)+i3,kOp(3))+iAOst(3)
               Do i2=1,jCmp
                  jSO=iAOtSO(iAO(2)+i2,kOp(2))+iAOst(2)
                  nijkl=0
                  Do lSOl=lSO,lSO+lBas-1
                     ll=iShlSO(lSOl)
                     kl0=nBasSh(iShlK)*(ll-1)
                     Do kSOk=kSO,kSO+kBas-1
                        kk=iShlSO(kSOk)
                        kl=kl0+kk
                        If (iCol(kl,iShlKL).gt.0) Then
                           Do jSOj=jSO,jSO+jBas-1
                              jj=iShlSO(jSOj)
                              nijkl=nijkl+1
                              If (iRow(jj,iShlJ).gt.0) Then
                                 Jkl=nAB*(iCol(kl,iShlKL)-1)
     &                              +iRow(jj,iShlJ)
                                 TInt(Jkl)=AOInt(nijkl,i1,i2,i3,i4)
                              End If
                           End Do
                        Else
                           nijkl=nijkl+jBas
                        End If
                     End Do
                  End Do
               End Do
            End Do
         End Do
      Else
         Call WarningMessage(2,
     &                    SecNam//': shell reordering not implemented!')
         Call LDF_Quit(1)
      End If

#ifndef _DEBUGPRINT_
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iBas)
#endif
      End
