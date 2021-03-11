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
      SUBROUTINE CHO_FINAL(WriteBookmarks)
C
C     Purpose: Cholesky finalizations.
C
      Implicit None
      Logical WriteBookmarks
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "choini.fh"
#include "chobkm.fh"
#include "WrkSpc.fh"

      INTEGER CHOISINI, IREO
      INTEGER NUMV(8)
      Integer ip, l
#if defined (_DEBUG_)
      Integer is1CCD
#endif

#if defined (_DEBUG_)
      CALL QENTER('_FINAL')
#endif

C     Write NUMCHO array, shell indices, and threshold to runfile.
C     ------------------------------------------------------------

      CALL CHO_P_GETGV(NUMV,NSYM)
      CALL PUT_IARRAY('NUMCHO',NUMV,NSYM)
      CALL PUT_IARRAY('iSOShl',IWORK(ip_ISOSHL),NBAST)
      CALL PUT_DSCALAR('Cholesky Threshold',THRCOM)
#if defined (_DEBUG_)
      ! This is needed in order for bookmark tests in cho_x_init to work
      If (WriteBookmarks) Then
         If (Cho_1Center) Then
            is1CCD=1
         Else
            is1CCD=0
         End If
         Call Put_iScalar('1C-CD',is1CCD)
      End If
#endif

C     Write bookmarks to runfile.
C     First, transpose array.
C     ---------------------------

      If (WriteBookmarks) Then
         l=4
         Call GetMem('BkmDim','Allo','Inte',ip,l)
         iWork(ip)=nCol_BkmVec
         iWork(ip+1)=nRow_BkmVec
         iWork(ip+2)=nCol_BkmThr
         iWork(ip+3)=nRow_BkmThr
         Call Put_iArray('Cholesky BkmDim',iWork(ip),l)
         Call GetMem('BkmDim','Free','Inte',ip,l)
         If (nRow_BkmVec.gt.0 .and. nCol_BkmVec.gt.0 .and.
     &       nRow_BkmThr.gt.0 .and. nCol_BkmThr.gt.0) Then
            l=nRow_BkmVec*nCol_BkmVec
            Call GetMem('Scratch','Allo','Inte',ip,l)
            Call iTrnsps(nSym,nCol_BkmVec,iWork(ip_BkmVec),iWork(ip))
            Call Put_iArray('Cholesky BkmVec',iWork(ip),l)
            Call GetMem('Scratch','Free','Inte',ip,l)
            Call GetMem('BkmVec','Free','Inte',ip_BkmVec,l_BkmVec)
            ip_BkmVec=0
            l_BkmVec=0
            nRow_BkmVec=0
            nCol_BkmVec=0
            l=nRow_BkmThr*nCol_BkmThr
            Call GetMem('Scratch','Allo','Real',ip,l)
            Call Trnsps(nSym,nCol_BkmThr,Work(ip_BkmThr),Work(ip))
            Call Put_dArray('Cholesky BkmThr',Work(ip),l)
            Call GetMem('Scratch','Free','Real',ip,l)
            Call GetMem('BkmVec','Free','Real',ip_BkmThr,l_BkmThr)
            ip_BkmThr=0
            l_BkmThr=0
            nRow_BkmThr=0
            nCol_BkmThr=0
         End If
      End If
      If (l_BkmVec.gt.0) Then
         Call GetMem('BkmVec','Free','Inte',ip_BkmVec,l_BkmVec)
         ip_BkmVec=0
         l_BkmVec=0
         nRow_BkmVec=0
         nCol_BkmVec=0
      End If
      If (l_BkmThr.gt.0) Then
         Call GetMem('BkmThr','Free','Real',ip_BkmThr,l_BkmThr)
         ip_BkmThr=0
         l_BkmThr=0
         nRow_BkmThr=0
         nCol_BkmThr=0
      End If

C     Write vector file address mode to runfile.
C     ------------------------------------------

      CALL PUT_ISCALAR('ChoVec Address',CHO_ADRVEC)

C     Write reorder mark to runfile.
C     ------------------------------

      IF (CHO_REORD) THEN
         IREO = 1
      ELSE
         IREO = 0
      END IF
      CALL PUT_ISCALAR('Cholesky Reorder',IREO)

C     Set initialization integer flag to "not set".
C     ---------------------------------------------

      CHOISINI = CHOINICHECK + 1
      CALL PUT_ISCALAR('ChoIni',CHOISINI)

#if defined (_DEBUG_)
      CALL QEXIT('_FINAL')
#endif

      END
