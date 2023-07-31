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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      Subroutine Rd2Int_SCF()
************************************************************************
*                                                                      *
*     purpose: Read basis set informations from two-electron file      *
*              and compare them with those read from 2-el. file        *
*                                                                      *
************************************************************************
      use InfSCF, only: nSym, nBas, nSkip
      Implicit None
#include "Molcas.fh"
*
*---- Define local variables
      Integer iRC, nSymX, iSym
      Integer nBasX(MxSym)
      Logical SqI2

      iRc=-1
      Call GetOrd(iRc,SqI2,nSymX,nBasX,nSkip)
      If (iRc.ne.0) Then
         Write (6,*) 'The program failed to read the header of ORDINT.'
         Call Abend
      End If
*
      If (nSymX.ne.nSym) Then
         Write (6,*) 'nSymX.ne.nSym, nSymX, nSym=',nSymX,nSym
         Call Abend
      End If
      Do iSym = 1, nSym
         If (nBas(iSym).ne.nBasX(iSym)) Then
            Write (6,*) 'nBas(iSym).ne.nBasX(iSym)'
            Write (6,*) 'nBas=',nBas
            Write (6,*) 'nBasX=',nBasX
            Call Abend
         End If
      End Do

      End subroutine Rd2Int_SCF
