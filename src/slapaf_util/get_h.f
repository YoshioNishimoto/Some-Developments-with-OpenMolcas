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
* Copyright (C) 2010, Roland Lindh                                     *
************************************************************************
      Subroutine Get_H(ip_H)
************************************************************************
*                                                                      *
*     Object: to get the force constant matrix in cartesians           *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh                                             *
*             Department of Quantum Chemistry                          *
*             Uppsala University, Sweden                               *
*             October 2010                                             *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "info_slapaf.fh"
#include "WrkSpc.fh"
      Logical Found
*
*define _DEBUG_
      nX=3*nsAtom
      mInter=mInt + mTROld
      nInter=mInt
*
      Call Allocate_Work(ip_H,nX**2)
      Call Allocate_Work(ipTmp2,nX**2)
      Call Allocate_Work(ipH,nInter**2)
*---- If there is an updated Hessian in the runfile, use it;
*     otherwise use the last computed one.
*     (Hss_upd must be removed every time Hss_Q is added)
      Call Qpg_dArray('Hss_upd',Found,nHss)
      If (Found.and.(nHss.eq.nInter**2)) Then
         Call Get_dArray('Hss_upd',Work(ipH),nInter**2)
      Else
         Call Get_dArray('Hss_Q',Work(ipH),nInter**2)
      End If
      Call Allocate_Work(ipBOld,nX*nInter)
*---- If there is an old BMx stored, use it;
*     otherwise use the current BMx
      Call Qpg_dArray('BMxOld',Found,nBMx)
      If (Found.and.(nBMx.eq.nX*nInter)) Then
         Call Get_dArray('BMxOld',Work(ipBOld),nX*nInter)
      Else
         Call Get_dArray('BMtrx',Work(ipBOld),nX*nInter)
      End If
*
      Call Get_H_(nX,Work(ipBOld),mInter,nInter,Work(ipH),Degen,
     &            Work(ipTmp2),Work(ip_H),Smmtrc,nsAtom)
*
      Call Free_Work(ipBOld)
      Call Free_Work(ipH)
      Call Free_Work(ipTmp2)
*
      Return
      End
      Subroutine Get_H_(nX,BMtrx,mInter,nInter,H,Degen,
     &                 Tmp2,Tmp3,Smmtrc,nAtom)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
      Logical Smmtrc(3,nAtom)
      Real*8 BMtrx(nX,nInter), H(nInter,nInter),Degen(3,nAtom),
     &       Tmp2(nX**2), Tmp3(nX**2)
*                                                                      *
************************************************************************
*                                                                      *
*     The Hessian matrix (zero elements over translations and
*     rotations are explicitly excluded).
*
#ifdef _DEBUG_
      Call RecPrt('BMtrx',' ',BMtrx,nX,nInter)
#endif
      ii = 0
      Do i = 1, nAtom
         Do ix = 1, 3
*
            If (Smmtrc(ix,i)) Then
               iix = (i-1)*3 + ix
               ii = ii + 1
               Do j = 1, nInter
                  ij = (j-1)*mInter + ii

                  tmp_ij=Zero
                  Do k = 1, nInter
                     tmp_ij=tmp_ij+BMtrx(iix,k)*H(k,j)
                  End Do
                  Tmp2(ij)=tmp_ij

               End Do
            End If
*
         End Do
      End Do
*
      Do i = 1, mInter
*
         jj = 0
         Do j = 1, nAtom
            Do jx = 1, 3
*
               If (Smmtrc(jx,j)) Then
                  jjx=(j-1)*3+jx
                  jj = jj + 1
                  ij = (jj-1)*mInter + i
                  tmp_ij=Zero
                  Do k = 1, nInter
                     ik = (k-1)*mInter + i
                     tmp_ij=tmp_ij+Tmp2(ik)*BMtrx(jjx,k)
                  End Do
                  Tmp3(ij)=tmp_ij
               End If
*
            End Do
         End Do
*
      End Do
*
#ifdef _DEBUG_
      Call RecPrt('Hessian (cartesian)',' ',Tmp3,mInter,mInter)
#endif
      Call Put_dArray('FC-Matrix',Tmp3,mInter**2)
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real_array(Degen)
      End
