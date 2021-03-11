************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Init_PPList
      Implicit Real*8 (a-h,o-z)
#include "tlist.fh"
#include "para_info.fh"
#include "WrkSpc.fh"
#include "status.fh"
      Logical Debug
      Data Debug/.False./
*
      If (Debug) Then
         If (PP_Status.eq.Active) Then
            Write (6,*) 'Init_PPList: Active'
         Else If (PP_Status.eq.InActive) Then
            Write (6,*) 'Init_PPList: Active'
         Else
            Write (6,*) 'Init_PPList: undefined!'
         End If
      End If
      If (PP_Status.eq.Active) Return
      PP_Status=Active
*
c     Write (*,*) 'Init_PPList'
      iTskCan=0
      mTasks=0
      iStrt_TList=0
      iEnd_TList=nTasks+1
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
chjw  now allocated by Alloc_tlist (init_tlist.f)
chjw  Call GetMem('TskLst','ALLO','INTE',ipTskL,nTasks*2)
      call izero(iWork(ipTskL),nTasks)
      Do iTsk = 0, nTasks-1
        iWork(ipTskL+iTsk)=MOD(iTsk+MyRank,nTasks)+1
      End Do
c     Write (*,*) (iWork(ipTskL+iTsk),iTsk = 0, nTasks-1)
*
*---- Copy list in inverse order to ensure the proper order if
*     reinitiated at once.
*
      iE = nTasks-1
      call izero(iWork(ipTskL+nTasks),nTasks)
      Do i = 0, nTasks-1
         iWork(ipTskL+nTasks+iE) = iWork(ipTskL+i)
         iE = iE - 1
      End Do
*
      QLast(1)=Not_Used
      QLast(2)=Not_Used
*
      Return
      End
      Subroutine ReInit_PPList(Semi_Direct)
      Implicit Real*8 (a-h,o-z)
#include "tlist.fh"
#include "para_info.fh"
#include "WrkSpc.fh"
#include "status.fh"
      Logical Debug,Semi_Direct
      Data Debug/.False./
*
c     Write (*,*) 'ReInit_PPList'
*
      If (Debug) Then
         If (PP_Status.eq.Active) Then
            Write (6,*) 'ReInit_PPList: Active'
         Else If (PP_Status.eq.InActive) Then
            Write (6,*) 'ReInit_PPList: Active'
         Else
            Write (6,*) 'ReInit_PPList: undefined!'
         End If
      End If
      If (PP_Status.ne.Active) Then
         Write (6,*) 'ReInit_PPList: List is not active!'
         Call Abend
      End If
      iTskCan=0
      mTasks=iStrt_TList
      If (nProcs.eq.1) Then
         iStrt_TList=0
         iEnd_TList=nTasks+1
         Return
      End If
*
      If (Semi_Direct) Then
*
*---- Copy first the task indices of tasks that were exectuted
         Call ICopy(mTasks,iWork(ipTskL+nTasks),1,iWork(ipTskL),1)
c     Write (*,*) 'mTasks=',mTasks
*
*---- Now copy task indices of tasks which were not executed by this node.
*     Change the order so that the first task it the largest in the list.
*
         iE = mTasks
         iCount = 1
         Do i = mTasks, nTasks-1
            If (iCount .gt. myRank) Then
               iWork(ipTskL+i) = iWork(ipTskL+nTasks+i)
            Else
               iWork(ipTskL+i) = iWork(ipTskL+nTasks+iE)
               iE = iE - 1
               iCount = iCount + 1
            Endif
         Enddo

c        Write (*,*) (iWork(ipTskL+iTsk),iTsk = 0, nTasks-1)
c        Write (*,*) 'mTasks=',mTasks
*
      End If
*
      iStrt_TList=0
      iEnd_TList=nTasks+1
*
      QLast(1)=Not_Used
      QLast(2)=Not_Used
*
      Return
      End
*
      Subroutine Free_PPList
#include "tlist.fh"
#include "para_info.fh"
#include "WrkSpc.fh"
#include "status.fh"
*
      If (PP_Status.ne.Active) Return
      PP_Status=Inactive
*
      If (.Not. Is_Real_Par() .OR. nProcs.eq.1) Return
C     Call GetMem('TskLst','Free','INTE',ipTskL,nTasks_alloc*2)
      Call Free_iWork(ipTskL)
*
      Return
      End
