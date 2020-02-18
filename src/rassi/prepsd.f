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
      SUBROUTINE PREPSD(WFTP,ISGSTR,ICISTR,LSYM,
     &                  ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,
     &                  NCONF,CI,DET,rdetcoeff,rdetocc,ndet,norb)
      IMPLICIT NONE
      INTEGER ISGSTR(*),ICISTR(*)
      INTEGER ICNFTAB(*),ISPNTAB(*),ISSTAB(*),IFSBTAB(*)
      INTEGER LSYM,NCONF
      REAL*8 CI(*),DET(*)
      INTEGER IMODE,LCTMP
      CHARACTER*8 WFTP
CC VK/GG 2020 CC
      real*8 :: rdetcoeff(10000000)
      character*99 :: rdetocc(10000000)
      integer :: ndet, norb
CC CC
#include "WrkSpc.fh"
C Purpose: Given a RASSCF wave function in Split-GUGA format
C and an orbital transformation matrix for the purpose of
C getting biorthonormal orbitals, prepare a wave function
C in the general SD format, using transformed orbitals.

      IF(WFTP.EQ.'GENERAL ') THEN
C Transform SGUGA to SymmG:
        CALL GETMEM('PREPSD','ALLO','REAL',LCTMP,NCONF)
        IMODE=1
        CALL SYG2SGU(IMODE,ISGSTR,ICISTR,LSYM,ICNFTAB,ISPNTAB,
     &                  CI,WORK(LCTMP))
C Transform SymmG to Slater Dets:
        CALL SYGTOSD(ICNFTAB,ISPNTAB,ISSTAB,IFSBTAB,WORK(LCTMP),DET,
     &               rdetcoeff,rdetocc,ndet,norb)
        CALL GETMEM('PREPSD','FREE','REAL',LCTMP,NCONF)
      ELSE
        DET(1)=CI(1)
CC VK/GG 2020 CC in case of HISPIN type of WF we have to set the number of determinants to one here
        NDET=1
CC CC
      END IF
      RETURN
      END
