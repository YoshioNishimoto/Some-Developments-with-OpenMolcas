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
* Copyright (C) 1990,1994,1995, Jeppe Olsen                            *
************************************************************************
      SUBROUTINE FREESTR
*
*
* Free pointers for saving information about strings and
* their mappings
*
*========
* Input :
*========
* Number and types of strings defined by /STRINP/
* Symmetry information stored in         /CSM/
* String information stored in           /STINF/
*=========
* Output
*=========
* Pointers stored in common block /STRBAS/
*
* Jeppe Olsen , Winter of 1990
*
* Last Revision , Dec 24 1990 , Almaden
*
* Updated with iuniqtp, dec 11, 1994
* Modified for deallocation, Sept. 25, 2005.
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "detdim.fh"
#include "WrkSpc.fh"
#include "orbinp_mclr.fh"
#include "strinp_mclr.fh"
#include "strbas_mclr.fh"
#include "csm.fh"
#include "stinf_mclr.fh"
*. Start of string information
*     CALL MEMMAN(KSTINF,IDUMMY,'FREE  ',IDUMMY,'DUMMY ')
      NTEST = 0
      IF(NTEST.NE.0)
     &WRITE(6,*) ' First word with string information',KSTINF
* =====================================================================
*
* 1 : String information
*
* =====================================================================
*
      nDum=1
      IIITEST = 1
      DO 10 ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
*.  Offsets for occupation of strings and reordering array
          Call GetMem('OCSTR ','Free','INTEGER',KOCSTR(ITYP),nDum)
          Call GetMem('STREO','Free','INTEGER',KSTREO(ITYP),nDum)
*. Symmetry and class of each string
          Call GetMem('STSM  ','Free','INTEGER',KSTSM(ITYP),nDum)
          Call GetMem('STCL  ','Free','INTEGER',KSTCL(ITYP),nDum)
        ELSE
          IITYP = - IUNIQTP(ITYP)
          KSTREO(ITYP) = KSTREO(IITYP)
          KSTSM(ITYP)  = KSTSM(IITYP)
          KSTCL(ITYP)  = KSTCL(IITYP)
        END IF
   10 CONTINUE
*. Number of strings per symmetry and occupation
      DO ITYP = 1, NSTTYP
        IF(IUNIQTP(ITYP).EQ.ITYP) THEN
        Call GetMem('NSTSO ','Free','INTEGER',KNSTSO(ITYP),nDum)
*. Offset of strings per symmetry and occupation
        Call GetMem('ISTSO ','Free','INTEGER',KISTSO(ITYP),nDum)
*. Number of electrons in RAS1 and RAS3 per sub type, is sub-type active
        Call GetMem('IEL1  ','Free','INTEGER',KEL1(ITYP),nDum)
        Call GetMem('IEL3  ','Free','INTEGER',KEL3(ITYP),nDum)
        Call GetMem('ACTP ','Free','INTEGER',KACTP(ITYP),nDum)
CMS: New array introduced according to Jeppes new strinfo representation
        Call GetMem('KEL123','Free','INTEGER',KEL123(ITYP),nDum)
**. Lexical adressing of arrays: NB! Not allocated here in Jeppes new version!
        Call GetMem('Zmat  ','Free','INTEGER',KZ(ITYP),nDum)
        ELSE
*. redirect
          IITYP = - IUNIQTP(ITYP)
          KNSTSO(ITYP) = KNSTSO(IITYP)
          KISTSO(ITYP) = KISTSO(IITYP)
          KEL1(ITYP)   = KEL1(IITYP)
          KEL3(ITYP)   = KEL3(IITYP)
          KACTP(ITYP)  = KACTP(IITYP)
          KZ(ITYP)     = KZ(IITYP)
          KEL123(ITYP) = KEL123(IITYP)
        END IF
      END DO
*. Mappings between different string types
      DO ITYP = 1, NSTTYP
          NSTRIN = NSTFTP(ITYP)
          IF(ISTAC(ITYP,2).NE.0.AND.ISTAC(ITYP,1).NE.0) THEN
*.creation on string allowed , use full orbital notation
            LENGTH = NACOB*NSTRIN
*. No explicit offset or length. NEW:
            KSTSTMI(ITYP) = 0
            KSTSTMN(ITYP) = 0
          ELSE IF(ISTAC(ITYP,1).NE.0.AND.ISTAC(ITYP,2).EQ.0) THEN

*. only annihilation allowed, use compact scheme
            LENGTH = NELEC(ITYP)*NSTRIN
*. No explicit offset or length. NEW:
            KSTSTMI(ITYP) = 0
            KSTSTMN(ITYP) =  0
CMS: New else block
          ELSE IF (ISTAC(ITYP,1).EQ.0.AND.ISTAC(ITYP,2).NE.0) THEN
*. Only creation allowed, use compact scheme with offsets
*
*. Explicit offsets and lengths
            Call GetMem('STSTMI','Free','INTEGER',KSTSTMI(ITYP),nDum)
            Call GetMem('STSTMN','Free','INTEGER',KSTSTMN(ITYP),nDum)
          END IF
*. has this map been constructed before ?
          IIIITEST = 0
          IF(IUNIQTP(ITYP).EQ.ITYP.OR.IIIITEST.EQ.1) THEN
            IMNEW = 1
            IUNIQMP(ITYP) = ITYP
          ELSE
*. check type of previous map
            DO JJTYP = 1, ITYP-1
            IITYP = -IUNIQTP(ITYP)
            IF(ABS(IUNIQTP(JJTYP)).EQ.IITYP.AND.
     &          IUNIQMP(JJTYP).EQ.JJTYP) THEN
            IF((ISTAC(ITYP,1).EQ.0.AND.ISTAC(JJTYP,1).EQ.0).OR.
     &         (ISTAC(ITYP,1).NE.0.AND.ISTAC(JJTYP,1).NE.0.AND.
     &          ABS(IUNIQTP(ISTAC(ITYP,1))).EQ.
     &          ABS(IUNIQTP(ISTAC(JJTYP,1)))) ) THEN
                IANEQ = 1
            ELSE
                IANEQ = 0
            END IF
            IF((ISTAC(ITYP,2).EQ.0.AND.ISTAC(JJTYP,2).EQ.0).OR.
     &         (ISTAC(ITYP,2).NE.0.AND.ISTAC(JJTYP,2).NE.0.AND.
     &          ABS(IUNIQTP(ISTAC(ITYP,2))).EQ.
     &          ABS(IUNIQTP(ISTAC(JJTYP,2)))) ) THEN
                ICREQ = 1
            ELSE
                ICREQ = 0
            END IF
            IF(IANEQ.EQ.1.AND.ICREQ.EQ.1) THEN
              IMNEW = 0
              IUNIQMP(ITYP) = -JJTYP
              GOTO 1211
            END IF
*
            END IF
            END DO
*. Normal exit from DO loop only if no identical map was found
            IMNEW = 1
            IUNIQMP(ITYP) = ITYP
 1211       CONTINUE
          END IF
          IF(IMNEW.EQ.1) THEN
            Call GetMem('CREMAP','Free','INTE',KSTSTM(ITYP,1),nDum)
            Call GetMem('ANNMAP','Free','INTE',KSTSTM(ITYP,2),nDum)
          ELSE
            KSTSTM(ITYP,1) = KSTSTM(-IUNIQMP(ITYP),1)
            KSTSTM(ITYP,2) = KSTSTM(-IUNIQMP(ITYP),2)
          END IF
      END DO
*. Symmetry of conjugated orbitals and orbital excitations
*     KCOBSM,KNIFSJ,KIFSJ,KIFSJO
      Call GetMem('Cobsm ','Free','INTEGER',KCOBSM,nDum)
      Call GetMem('Nifsj ','Free','INTEGER',KNIFSJ,nDum)
      Call GetMem('Ifsj  ','Free','INTEGER',KIFSJ,nDum)
      Call GetMem('Ifsjo ','Free','INTEGER',KIFSJO,nDum)
*. Symmetry of excitation connecting  strings of given symmetry
      Call GetMem('Ststx ','Free','INTEGER',KSTSTX,nDum)
*
**. Up and down mappings of strings containing the same number of electrons
*
      DO 70 ITYP = 1, NSTTYP
       IF(INUMAP(ITYP).NE.0)
     &Call GetMem('Numup ','Free','INTEGER',KNUMAP(ITYP),nDum)
       IF(INDMAP(ITYP).NE.0)
     &Call GetMem('Ndmup ','Free','INTEGER',KNDMAP(ITYP),nDum)
   70 CONTINUE
*
      RETURN
      END
