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
* Copyright (C) 1998, Roland Lindh                                     *
************************************************************************
      Subroutine Nr_Shells(nSkal)
************************************************************************
*                                                                      *
*     Object: to compute the number of unique shells in the input.     *
*                                                                      *
*     Author: Roland Lindh, Chemical Physics, University of Lund,      *
*             Sweden. January '98.                                     *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "Basis_Mode_Parameters.fh"
#include "Basis_Mode.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Determine the number of shells
*
      nSkal=0
      If (Basis_Mode.ne.Valence_Mode .and.
     &    Basis_Mode.ne.Auxiliary_Mode .and.
     &    Basis_Mode.ne.Fragment_Mode .and.
     &    Basis_Mode.ne.With_Auxiliary_Mode .and.
     &    Basis_Mode.ne.With_Fragment_Mode .and.
     &    Basis_Mode.ne.All_Mode) Then
         Call WarningMessage(2,'Nr_Shells: illegal Basis_Mode')
         Call Abend()
      End If
*
C      write (*,*) "in nr_shells" 
C      write (*,*) "atomic = ", atomic
      If (Atomic) Go To 300
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Molecular set up                                                 *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
C     write (*,*) "nCnttp = ", nCnttp
      Do iCnttp = 1, nCnttp
         nTest = nVal_Shells(iCnttp)-1
C        write (*,*) "ntest,nCntr(iCnttp) = ", ntest, nCntr(iCnttp)
         Do iCnt = 1, nCntr(iCnttp)
*
            Do 200 iAng=0, nTest
               iShll = ipVal(iCnttp) + iAng
               If (nExp(iShll).eq.0) Go To 200
               If (nBasis(iShll).eq.0) Go To 200
*
               If (Basis_Mode.eq.Valence_Mode .and.
     &             (AuxShell(iShll).or.FragShell(iShll))) Go To 200
               If (Basis_Mode.eq.Auxiliary_Mode .and.
     &             .Not.AuxShell(iShll)) Go To 200
               If (Basis_Mode.eq.Fragment_Mode .and.
     &             .Not.FragShell(iShll)) Go To 200
               If (Basis_Mode.eq.With_Auxiliary_Mode .and.
     &             FragShell(iShll)) Go To 200
               If (Basis_Mode.eq.With_Fragment_Mode .and.
     &             AuxShell(iShll)) Go To 200
               nSkal = nSkal + 1
*
 200        Continue                     ! iAng
         End Do                          ! iCnt
      End Do                             ! iCnttp
*
      Return
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Atomic set up                                                    *
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
 300  Continue
*
      Do iCnttp = kCnttp, lCnttp
      nTest = nVal_Shells(iCnttp)-1
      Do 400 iAng=0, nTest
         iShll = ipVal(iCnttp) + iAng
         If (nExp(iShll).eq.0) Go To 400
         If (nBasis(iShll).eq.0) Go To 400
*
         If (FragShell(iShll)) Go To 400
         nSkal = nSkal + 1
*
 400  Continue                     ! iAng
      End Do
      If (AuxCnttp(kCnttp)) nSkal=nSkal+1 ! Add dummy shell
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
