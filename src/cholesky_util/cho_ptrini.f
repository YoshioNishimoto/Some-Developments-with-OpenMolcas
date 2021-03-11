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
      SUBROUTINE CHO_PTRINI(irc)
C
C     Purpose: set all entries in /CHMIND/ zero.
C
      IMPLICIT NONE
      Integer irc
#include "choptr.fh"
      Integer nAlloc

      nAlloc = 0  ! allocation counter

      ip_INFRED = 0
      l_INFRED  = 0
      nAlloc    = nAlloc + 1

      ip_INFVEC = 0
      l_INFVEC  = 0
      nAlloc    = nAlloc + 1

      ip_INDRED = 0
      l_INDRED  = 0
      nAlloc    = nAlloc + 1

      ip_INDRSH = 0
      l_INDRSH  = 0
      nAlloc    = nAlloc + 1

      ip_ISCR = 0
      l_ISCR  = 0
      nAlloc  = nAlloc + 1

      ip_IIBSTRSH = 0
      l_IIBSTRSH  = 0
      nAlloc      = nAlloc + 1

      ip_NNBSTRSH = 0
      l_NNBSTRSH  = 0
      nAlloc      = nAlloc + 1

      ip_INTMAP = 0
      l_INTMAP  = 0
      nAlloc    = nAlloc + 1

      ip_NDIMRS = 0
      l_NDIMRS  = 0
      nAlloc    = nAlloc + 1

      ip_IRS2F = 0
      l_IRS2F  = 0
      nAlloc   = nAlloc + 1

      ip_iSOShl = 0
      l_iSOShl  = 0
      nAlloc    = nAlloc + 1

      ip_iShlSO = 0
      l_iShlSO  = 0
      nAlloc    = nAlloc + 1

      ip_iQuab = 0
      l_iQuab  = 0
      nAlloc   = nAlloc + 1

      ip_iBasSh = 0
      l_iBasSh  = 0
      nAlloc    = nAlloc + 1

      ip_nBasSh = 0
      l_nBasSh  = 0
      nAlloc    = nAlloc + 1

      ip_nBstSh = 0
      l_nBstSh  = 0
      nAlloc    = nAlloc + 1

      ip_iAtomShl = 0
      l_iAtomShl  = 0
      nAlloc      = nAlloc + 1

      ip_iSP2F = 0
      l_iSP2F  = 0
      nAlloc   = nAlloc + 1

      irc = CHO_NALLOC - nAlloc

      End
