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
* Copyright (C) 2020, Roland Lindski                                   *
************************************************************************
      SubRoutine goDelocalize(EVal,EVec,n,nB)
************************************************************************
*                                                                      *
*     purpose: Rotate degenerate eigenvectors to optimize              *
*              delocalization.                                         *
*                                                                      *
*                                                                      *
*     input:                                                           *
*       EVal    : the set of eigenvalues in random order               *
*       EVec    : the set of eigenvectors in random order              *
*       n,nB    : dimensions                                           *
*                                                                      *
*     output:                                                          *
*       EVal    : sorted set of eigenvalues                            *
*       EVec    : sorted set of eigenvectors                           *
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     Roland Lindh                                                     *
*     University of Uppsala, Sweden 2020                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit None
*
      Integer nB, N
      Real*8 EVal(n),EVec(nB,n)
*
      Integer nAtom, i, j, k
      Real*8 E1, E2, Test
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call Get_nAtoms_All(nAtom)
      If (nAtom.eq.1) Return
      Call RecPrt('Eval',' ',Eval,1,n)
#ifdef _OLD_CODE_
      i = 1
 666  Continue
      E1=EVal(i)
!
!     Do not try to delocalize these cases, probably fake degeneracy
      If (Abs(E1).lt.1.0D-10.or.E1.gt.5.0D0) Then
         i=i+1
         Go To 666
      End If
      Do j = i+1, n
         E2=EVal(j)
         Test=2.0D0* Abs(E1-E2)/Abs(E1+E2)
         If (Test.gt.1.0D-12) Then
            If (j.eq.i+1) Then
               i=j
               Go To 666
            Else
               Write (6,*) 'Eigenvalues',(EVal(k),k=i,j-1)
               k=j-i
               Call canonical_degenerate_eigenvectors(EVec(1,i),nB,k)
               i=j
               Go To 666
            End if
         End If
      End Do
!
#else
      i=1
      Do
         If (i.eq.n) Exit
         E1=EVal(i)
!
!        Do not try to delocalize these cases, probably fake degeneracy
         If (Abs(E1).lt.1.0D-10.or.E1.gt.5.0D0) Then
            i=i+1
            Cycle
         End If
         Do j = i+1, n
            E2=EVal(j)
            Test=2.0D0* Abs(E1-E2)/Abs(E1+E2)
            If (Test.gt.1.0D-12) Then
               If (j.eq.i+1) Then
                  i=j
                  Exit
               Else
                  Write (6,*) 'Eigenvalues',(EVal(k),k=i,j-1)
                  k=j-i
                  Call canonical_degenerate_eigenvectors(EVec(1,i),nB,k)
                  i=j
                  Exit
               End if
            End If
         End Do
!
      End Do
#endif
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
