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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  CD_Tester
*
*> @brief
*>   Test the decomposition modules
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Runs and tests the output from the general Cholesky decomposers
*> ::ChoDec (out-of-core) and ::CD_InCore (in-core). The positive
*> definite matrix to which \p ip_PDM points should be stored as a
*> lower triangle.
*>
*> Error codes:
*> - \p irc = ``-1``: \p n non-positive (&rarr; nothing done)
*> - \p irc =  ``0``: all OK
*> - \p irc =  ``1``: error in ::ChoDec
*> - \p irc =  ``2``: error in ::CD_InCore
*> - \p irc =  ``3``: error in both
*>
*> @param[out] irc     Return code
*> @param[in]  ip_PDM  Pointer to matrix in \p Work
*> @param[in]  n       Dimension of matrix (\p n &times; \p n)
*> @param[in]  Verbose Print flag
************************************************************************
      SubRoutine CD_Tester(irc,ip_PDM,n,Verbose)
      Implicit Real*8 (a-h,o-z)

      External CD_Tester_Col
      External CD_Tester_Vec
      Logical  Verbose

      Character*9 SecNam
      Parameter (SecNam = 'CD_Tester')

      Logical Restart

      Common / CDTHLP / ip_Mat, l_Mat, ip_Vec, l_Vec

#include "WrkSpc.fh"

      Call qEnter(SecNam)
      irc = 0
      If (n .lt. 1) Then
         If (Verbose) Then
            Write(6,*) SecNam,': nothing tested! Dimension is: n = ',n
            Call xFlush(6)
         End If
         Go To 1
      End If

C     Test ChoDec.
C     ============

      Restart = .false.
      Thr     = 1.0d-12
      Span    = 1.0d-2
      MxQual  = max(n/10,1)
      NumCho  = 0

      irc_sav = irc
      If (Verbose) Then
         Write(6,*)
         Write(6,*) '          >>>>>>>>>>>><<<<<<<<<<<<'
         Write(6,*) '          >>>> Testing ChoDec <<<<'
         Write(6,*) '          >>>>>>>>>>>><<<<<<<<<<<<'
         Write(6,*)
         Call xFlush(6)
      End If

      l_Mat   = n*n
      l_Vec   = n*n
      Call GetMem('Matrix','Allo','Real',ip_Mat,l_Mat)
      Call GetMem('Vec','Allo','Real',ip_Vec,l_Vec)

      l_Buf   = n*(n+MxQual)
      l_Diag  = n
      l_Qual  = n*(MxQual+1)
      l_ES    = 6
      l_Pivot = n
      l_iQual = MxQual
      Call GetMem('Buf','Allo','Real',ip_Buf,l_Buf)
      Call GetMem('Diag','Allo','Real',ip_Diag,l_Diag)
      Call GetMem('Qual','Allo','Real',ip_Qual,l_Qual)
      Call GetMem('Err','Allo','Real',ip_ES,l_ES)
      Call GetMem('Pivot','Allo','Inte',ip_Pivot,l_Pivot)
      Call GetMem('iQual','Allo','Inte',ip_iQual,l_iQual)

      Call CD_Tester_CPPF(Work(ip_PDM),Work(ip_Mat),n)
      Call CD_Tester_Diag(Work(ip_PDM),Work(ip_Diag),n)
      jrc = 0
      Call ChoDec(CD_Tester_Col,CD_Tester_Vec,
     &            Restart,Thr,Span,MxQual,
     &            Work(ip_Diag),Work(ip_Qual),Work(ip_Buf),
     &            iWork(ip_Pivot),iWork(ip_iQual),
     &            n,l_Buf,
     &            Work(ip_ES),NumCho,
     &            jrc)
      If (jrc .eq. 0) Then
         Call DGEMM_('N','T',n,n,NumCho,
     &              -1.0d0,Work(ip_Vec),n,Work(ip_Vec),n,
     &              1.0d0,Work(ip_Mat),n)
         Call CD_Tester_ES(Work(ip_Mat),n,Work(ip_ES))
         Call CD_Tester_Diff(Work(ip_Mat),n,Work(ip_ES+3))
         Call CD_Tester_Final(jrc,NumCho,n,Thr,Work(ip_ES),Verbose)
         If (jrc .ne. 0) Then
            irc = irc + 1
         End If
      Else
         If (Verbose) Then
            Write(6,*) SecNam,': ChoDec returned error code ',jrc
            Write(6,*) '<=> decomposition failed, no test performed!'
            Call xFlush(6)
         End If
         irc = irc + 1
      End If
      If (Verbose) Then
         If (irc .eq. irc_sav) Then
            Write(6,*) 'ChoDec succeeded....'
         Else
            Write(6,*) 'ChoDec failed....'
         End If
         Call xFlush(6)
      End If

C     Test in-core decomposition.
C     ===========================

      If (Verbose) Then
         irc_sav = irc
         Write(6,*)
         Write(6,*) '          >>>>>>>>>>>><<<<<<<<<<<<<<<'
         Write(6,*) '          >>>> Testing CD_InCore <<<<'
         Write(6,*) '          >>>>>>>>>>>><<<<<<<<<<<<<<<'
         Write(6,*)
         Call xFlush(6)
      End If

      Call CD_Tester_CPPF(Work(ip_PDM),Work(ip_Mat),n)
      jrc    = 0
      NumCho = 0
      Call CD_InCore(Work(ip_Mat),n,Work(ip_Vec),n,NumCho,Thr,jrc)
      If (jrc .eq. 0) Then
         Call CD_Tester_CPPF(Work(ip_PDM),Work(ip_Mat),n)
         Call DGEMM_('N','T',n,n,NumCho,
     &              -1.0d0,Work(ip_Vec),n,Work(ip_Vec),n,
     &              1.0d0,Work(ip_Mat),n)
         Call CD_Tester_ES(Work(ip_Mat),n,Work(ip_ES))
         Call CD_Tester_Diff(Work(ip_Mat),n,Work(ip_ES+3))
         Call CD_Tester_Final(jrc,NumCho,n,Thr,Work(ip_ES),Verbose)
         If (jrc .ne. 0) Then
            irc = irc + 2
         End If
      Else
         If (Verbose) Then
            Write(6,*) SecNam,': CD_InCore returned error code ',jrc
            Write(6,*) '<=> decomposition failed, no test performed!'
            Call xFlush(6)
         End If
         irc = irc + 2
      End If
      If (Verbose) Then
         If (irc .eq. irc_sav) Then
            Write(6,*) 'CD_InCore succeeded....'
         Else
            Write(6,*) 'CD_InCore failed....'
         End If
         Write(6,*)
         Call xFlush(6)
      End If

      Call GetMem('iQual','Free','Inte',ip_iQual,l_iQual)
      Call GetMem('Pivot','Free','Inte',ip_Pivot,l_Pivot)
      Call GetMem('Err','Free','Real',ip_ES,l_ES)
      Call GetMem('Qual','Free','Real',ip_Qual,l_Qual)
      Call GetMem('Diag','Free','Real',ip_Diag,l_Diag)
      Call GetMem('Buf','Free','Real',ip_Buf,l_Buf)
      Call GetMem('Vec','Free','Real',ip_Vec,l_Vec)
      Call GetMem('Matrix','Free','Real',ip_Mat,l_Mat)

    1   Continue
       Call qExit(SecNam)
      End
C
      SubRoutine CD_Tester_Col(Col,nDim,iCol,nCol,Buf,lBuf)
      Implicit Real*8 (a-h,o-z)
      Dimension Col(nDim,nCol), Buf(lBuf)
      Integer   iCol(nCol)

#include "WrkSpc.fh"

      Common / CDTHLP / ip_Mat, l_Mat, ip_Vec, l_Vec

      Do i = 1,nCol
         kOff = ip_Mat + nDim*(iCol(i) - 1)
         Call dCopy_(nDim,Work(kOff),1,Col(1,i),1)
      End Do

      If (lBuf .gt. 0) x = Buf(1) ! to make some compilers happy

      End
C
      SubRoutine CD_Tester_Vec(iVec1,nVec,Buf,lBuf,nDim,iOpt)
      Implicit Real*8 (a-h,o-z)
      Dimension Buf(lBuf)

      Character*13 SecNam
      Parameter (SecNam = 'CD_Tester_Vec')

#include "WrkSpc.fh"

      Common / CDTHLP / ip_Mat, l_Mat, ip_Vec, l_Vec

      If (iOpt .eq. 1) Then
         kOff = ip_Vec + nDim*(iVec1 - 1)
         lTot = nDim*nVec
         Call dCopy_(lTot,Buf,1,Work(kOff),1)
      Else If (iOpt .eq. 2) Then
         kOff = ip_Vec + nDim*(iVec1 - 1)
         lTot = nDim*nVec
         Call dCopy_(lTot,Work(kOff),1,Buf,1)
      Else
         Write(6,*)
         Write(6,*) 'WARNING! WARNING! WARNING! WARNING!'
         Write(6,*) SecNam,': illegal option: iOpt = ',iOpt
         Write(6,*) 'WARNING! WARNING! WARNING! WARNING!'
         Write(6,*)
         Call xFlush(6)
      End If

      End
C
      SubRoutine CD_Tester_CPPF(PDM,X,n)
      Implicit None
      Integer n
      Real*8  PDM(n*(n+1)/2), X(n,n)

      Integer i, j, iTri
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j

      Do j = 1,n
         X(j,j) = PDM(iTri(j,j))
         Do i = j+1,n
            X(i,j) = PDM(iTri(i,j))
            X(j,i) = X(i,j)
         End Do
      End Do

      End
C
      SubRoutine CD_Tester_Diag(PDM,Diag,n)
      Implicit None
      Integer n
      Real*8  PDM(n*(n+1)/2), Diag(n)

      Integer i, ii

      Do i = 1,n
         ii = i*(i-3)/2 + 2*i
         Diag(i) = PDM(ii)
      End Do

      End
C
      SubRoutine CD_Tester_Diff(X,n,Err)
      Implicit None
      Integer n
      Real*8  X(n*n), Err(3)

      Integer i
      Real*8  xn2
      Real*8  dble

      If (n .lt. 1) Then
         Err(1) =  9.876543210d15
         Err(2) = -9.876543210d15
         Err(3) =  9.876543210d15
      Else
         Err(1) = x(1)
         Err(2) = x(1)
         Err(3) = x(1)*x(1)
         Do i = 2,n*n
            Err(1) = min(Err(1),x(i))
            Err(2) = max(Err(2),x(i))
            Err(3) = Err(3) + x(i)*x(i)
         End Do
         xn2 = dble(n)*dble(n)
         Err(3) = Err(3)/xn2
      End If

      End
C
      SubRoutine CD_Tester_Final(irc,NumCho,n,Thr,Err,Verbose)
      Implicit None
      Integer irc, NumCho, n
      Real*8  Thr, Err(6)
      Logical Verbose

      Character*15 SecNam
      Parameter (SecNam = 'CD_Tester_Final')

      irc = 0

      If (Verbose) Then
         Write(6,*)
         Write(6,*) 'Final results from ',SecNam,':'
         Write(6,*) 'Matrix dimension: ',n
         Write(6,*) 'Number of vecs. : ',NumCho
         Write(6,*) 'Threshold       : ',Thr
         Write(6,*) 'Min. Diag. err. : ',Err(1)
         Write(6,*) 'Max. Diag. err. : ',Err(2)
         Write(6,*) 'RMS  Diag. err. : ',Err(3)
         Write(6,*) 'Min. Matr. err. : ',Err(4)
         Write(6,*) 'Max. Matr. err. : ',Err(5)
         Write(6,*) 'RMS  Matr. err. : ',Err(6)
      End If

      If (NumCho.lt.0 .or. NumCho.gt.n) Then
         irc = -1
         If (Verbose) Then
            Write(6,*) '>>> NumCho out of bounds!'
         End If
         Return
      End If

      If (abs(Err(1)) .gt. Thr) Then
         irc = irc + 1
         If (Verbose) Then
            Write(6,*) '>>> LARGE MINIMUM DIAGONAL ERROR: ',Err(1)
         End If
      End If
      If (abs(Err(2)) .gt. Thr) Then
         irc = irc + 1
         If (Verbose) Then
            Write(6,*) '>>> LARGE MAXIMUM DIAGONAL ERROR: ',Err(2)
         End If
      End If
      If (abs(Err(3)) .gt. Thr) Then
         irc = irc + 1
         If (Verbose) Then
            Write(6,*) '>>> LARGE RMS     DIAGONAL ERROR: ',Err(3)
         End If
      End If
      If (abs(Err(4)) .gt. Thr) Then
         irc = irc + 1
         If (Verbose) Then
            Write(6,*) '>>> LARGE MINIMUM MATRIX   ERROR: ',Err(4)
         End If
      End If
      If (abs(Err(5)) .gt. Thr) Then
         irc = irc + 1
         If (Verbose) Then
            Write(6,*) '>>> LARGE MAXIMUM MATRIX   ERROR: ',Err(5)
         End If
      End If
      If (abs(Err(6)) .gt. Thr) Then
         irc = irc + 1
         If (Verbose) Then
            Write(6,*) '>>> LARGE RMS     MATRIX   ERROR: ',Err(6)
         End If
      End If

      If (Verbose) Call xFlush(6)

      End
C
      SubRoutine CD_Tester_ES(X,n,Err)
      Implicit None
      Integer n
      Real*8  X(n,n), Err(3)
      Integer i

      If (n .lt. 1) Then
         Err(1) =  9.876543210d15
         Err(2) = -9.876543210d15
         Err(3) =  9.876543210d15
      Else
         Err(1) = X(1,1)
         Err(2) = X(1,1)
         Err(3) = X(1,1)*X(1,1)
         Do i = 1,n
            Err(1) = min(Err(1),X(i,i))
            Err(2) = max(Err(2),X(i,i))
            Err(3) = Err(3) + X(i,i)*X(i,i)
         End Do
         Err(3) = sqrt(Err(3)/dble(n))
      End If

      End
