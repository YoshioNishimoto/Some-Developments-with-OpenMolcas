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
      SubRoutine FckMat
************************************************************
*                                                          *
*   Driver for calculation of optimized fock matrix.       *
*                                                          *
************************************************************
      use Arrays, only: FAMO, FIMO, INT2
      implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "machine.fh"
*
*                                                                      *
************************************************************************
*                                                                      *
*...  Read density matrix
*
      nm=0
      nmm=0
      nmmm=0
      Do iS=1,nSym
         nm=nAsh(is)+nm
         nMM=Max(nMM,nAsh(is)+nIsh(iS))
         nMMM=Max(nmmM,nBas(is))
      End Do
      nAtri=nm*(nm+1)/2
      nAtri=nAtri*(nAtri+1)/2
      nmmm=((nmmm-1)/nRec+1)*nRec
      nmm=nmm*nMMM
      nmm=nmm**2
      Call GetMem('f0sqMO','Allo','Real',ipf0sqMO,ndens2)
      Call mma_allocate(FIMO,ndens2,Label='FIMO')
      If (iMethod.eq.2) Then
         Call mma_allocate(Int2,nAtri,Label='Int2')
      Else
         Call mma_allocate(Int2,1,Label='Int2')
      End If
      Int2(:)=0.0d0
      Call mma_allocate(FAMO,nDens2,Label='FAMO')
      Call GetMem('Temp4','ALLO','Real',ipQ,nDens2)
      Call GetMem('Temp2','Allo','Real',ipTmp2,2*ndens2)
      ipScr=ipTmp2+ndens2
      Call GetMem('Temp5','Allo','Real',ipT3,ndens2)
*                                                                      *
************************************************************************
*                                                                      *
*     Calculate two-electron contribution
*
      Call Read22_2(Int2,
     &              Work(ipF0SqMo),Work(ipQ),
     &              FIMO,FAMO,
     &              Work(ipTmp2),Work(ipScr),
     &              Work(ipT3))
*                                                                      *
************************************************************************
*                                                                      *
      Call GetMem('Temp5','FREE','Real',ipT3,ndens2)
      Call GetMem('Temp2','Free','Real',iptmp2,ndens2)
      Call GetMem('Temp4','Free','Real',ipQ,ndens2)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
