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
* Copyright (C) 1993, Markus P. Fuelscher                              *
************************************************************************
      Subroutine Rd2Int_RASSCF()
************************************************************************
*                                                                      *
*     Read header of the two-electron integral file                    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      use Definitions, only: LF => u6
      use UnixInfo, only: ProgName
      use General_data
      Implicit None
#include "rasscf.fh"
      Integer nSymX,nBasX(8),iRc,i,IERR,iSym
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      iRc=-1
      Call GetOrd(iRc,lSquare,nSymX,nBasX,nSkipX)
      If (iRc/=0) Then
        Write(LF,*)'RD2INT Error: Failed to read from ORDINT file.'
        Write(LF,*) Progname,
     &            ' tried to read two-electron integrals from'
        Write(LF,*)'the ORDINT file, but failed. Something is wrong'
        Write(LF,*)'with the file. Perhaps it is missing?'
        Call Quit_OnUserError()
      End If
      If ( nSymX/=nSym ) Then
        Write(LF,*)'RD2INT Error: Wrong size of symmetry group.'
        Write(LF,*)ProgName,' tried to use two-electron integrals from'
        Write(LF,*)'a file that was evidently created for some other'
        Write(LF,*)'program run.'
        Write(LF,'(1x,a,2i8)')'nSymX,nSym:',nSymX,nSym
        Call Quit_OnUserError()
      End If
      IERR=0
      Do iSym=1,nSym
         If ( nBas(iSym).ne.nBasX(iSym) ) IERR=1
      End Do
      If ( IERR==1 ) Then
        Write(LF,*)'RD2INT Error: Wrong nr of basis functions.'
        Write(LF,*)'RASSCF tried to use two-electron integrals from'
        Write(LF,*)'a file that was evidently created for some other'
        Write(LF,*)'program run.'
        Write(LF,'(1x,a,8i8)')'nBas :',(nBas(i),i=1,nSym)
        Write(LF,'(1x,a,8i8)')'nBasX:',(nBasX(i),i=1,nSym)
        Call Quit_OnUserError()
      End If
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      End Subroutine Rd2Int_RASSCF
