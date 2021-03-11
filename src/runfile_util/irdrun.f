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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine read a record from the runfile.                         *
* Data type is Real*8.                                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine iRdRun(Label,Data,nData)
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Character*(*) Label
      Integer       Data(*)
      Integer       nData
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Character*64  ErrMsg
      Integer       iRc
      Integer       iOpt
*----------------------------------------------------------------------*
* Call extended reading routine.                                       *
*----------------------------------------------------------------------*
      iRc=0
      iOpt=0
      Call ixRdRun(iRc,Label,Data,nData,iOpt)
      If(iRc.ne.0) Then
         Write(ErrMsg,'(3a)') 'Error reading field "',
     &                        Label,
     &                        '" from runfile'
         Call SysAbendMsg('iRdRun',ErrMsg,' ')
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
