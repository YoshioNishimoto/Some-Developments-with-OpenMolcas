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
* Copyright (C) Francesco Aquilante                                    *
*               2007, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE CHO_VTRA2(irc,scr,lscr,jVref,JVEC1,JNUM,NUMV,JSYM,
     &        IREDC,iSwap,nDen,kDen,CMO,ISTART,NUSE,ipChoT,iSkip)
* Compare original:
*      SUBROUTINE CHO_VTRA(irc,scr,lscr,jVref,JVEC1,JNUM,NUMV,JSYM,IREDC,
*     &                   iSwap,nDen,kDen,ipMOs,nPorb,ipChoT,iSkip)

*********************************************************
*   Author: F. Aquilante
*   Modified PAM 2007: Use ordinary CMO array without restructuring
*   This amounts to (1) having symmetry blocks of MO:s accessed
*   untransposed, (2) not through pointers to workspace array Work()
*   but giving the index of start orbital in each symmetry instead,
*   and (3) having CMO array in call parameter list.
*
*   Purpose:  SCR(lscr) contains JNUM cholesky vectors
*             starting from JVEC1 and stored in reduced
*             sets. The routine performs an MOs half transformation
*             of these elements in a set of target
*             arrays identified by the pointers ipChoT.
*             In the target arrays, the vectors are
*             stored in full dimension and as a
*             subset of a a given NUMV number of vectors.
*             Each pointer should thereby point to a
*             location where the corresponding Cholesky
*             vector of a given unique symmetry pair
*             of indices has to be stored
*
*   Input:
*       jVref =  index of the first vector to be transformed
*                computed wrt to the first vector in the
*                target arrays (see calling routine)
*
*       Jvec1 =  first vector to be transformed
*       JNum  =  # of vectors to be transformed
*
*       NumV  =  total # of vectors in the target arrays
*
*       nDen  =  total # of densities to which MOs refer
*       kDen  =  first density to be treated
*
*       iSwap :   = 0   L(k,b,J) is returned
*                 = 1   L(a,k,J) is returned
*                 = 2   L(k,J,b) is returned
*                 = 3   L(a,J,k) is returned
*
*       iSkip(syma)=0 : skip the symmetry block a.
*                    Any vector L(ak) with syma x symk=JSYM
*                    won't be returned in the target array
*
*       IREDC :  reduced set in core at the moment of
*                the first call to the routine.
*                Can be set to -1 by the calling routine
*
*********************************************************

      Implicit Real*8 (a-h,o-z)
      Real*8  Scr(lscr)
      Logical Debug
      Integer nDen,kDen
*      Integer ipChoT(8,*),ipMOs(8,*),iSkip(*),nPorb(8,*)
      Integer ipChoT(8,*),iSkip(*)
      Real*8 CMO(*)
      Integer IOFFC(8),ISTART(8,*),NUSE(8,*)


      Integer  cho_isao
      External cho_isao

      Character*8  SECNAM
      Parameter (SECNAM = 'CHO_VTRA')

      Real*8  Fac(0:1)
      Data Fac /0.5D0,1.0D0/

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      Parameter (N2 = InfVec_N2)

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      IndRed(i,k) = iWork(ip_IndRed-1+nnBstrT(1)*(k-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)
******
      iRS2F(i,j)  = iWork(ip_iRS2F-1+2*(j-1)+i)
************************************************************************

#ifdef _DEBUG_
      Debug=.true.
#else
      Debug=.false.
#endif

      Call qEnter(SECNAM)

**********************************************************
C
C    From Reduced sets to half-MOs full storage
C    ------------------------------------------
C
C    iSwap = 0
C
C     L{a,b,J} ---> L(p,b,J)   ! stride-1 transformation
C
C    iSwap = 1
C
C     L{a,b,J} ---> L(a,q,J)
C
C    iSwap = 2
C
C     L{a,b,J} ---> L(p,J,b)   ! stride-1 transformation
C
C    iSwap = 3
C
C     L{a,b,J} ---> L(a,J,q)
**********************************************************

      iLoc = 3 ! use scratch location in reduced index arrays

* PAM: Offset counter into CMO array:
      IOC=0
      DO ISYM=1,NSYM
       IOFFC(ISYM)=IOC
       IOC=IOC+NBAS(ISYM)**2
      END DO
* PAM Added to allow call directly from TRACHO
      Do iSymp=1,nSym
         IF (iSkip(iSymp).ne.0) THEN
            iSymb = muld2h(JSYM,iSymp)
            Do jDen=kDen,nDen
               If (iSwap.eq.0 .or. iSwap.eq.2) then ! Lpb,J or LpJ,b
                  Call FZero(Work(ipChoT(iSymp,jDen)),
     &                 nUse(iSymp,jDen)*nBas(iSymb)*NUMV)
               ElseIf (iSwap.eq.1 .or. iSwap.eq.3) then !Laq,J or LaJ,q
                  Call FZero(Work(ipChoT(iSymp,jDen)),
     &                 nUse(iSymb,jDen)*nBas(iSymp)*NUMV)
               EndIf
            End Do
         ENDIF
      End Do


      IF (iSwap.eq.0) THEN     ! L(pb,J) storage

         NREAD = 0
         kchot = 0
         DO JVEC=1,JNUM

            JRED = InfVec(JVEC1-1+JVEC,2,jSym)

            IF (JRED .NE. IREDC) THEN ! JRED is not the rs in core
               Call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
               IREDC = JRED
            END IF

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,JRED)

            IF (JSYM.eq.1) THEN  ! L(a,b,J)=L(b,a,J); only a.ge.b presnt

               Do jRab=1,nnBstR(jSym,iLoc)

                  kRab = iiBstr(jSym,iLoc) + jRab
                  iRab = IndRed(kRab,iLoc)

                  iag   = iRS2F(1,iRab)  !global address
                  ibg   = iRS2F(2,iRab)

                  iSyma = cho_isao(iag) !symmetry block; S(b)=S(a)=S(p)

                  kscr  = kscr + 1

                  IF (iSkip(iSyma).ne.0) THEN

                     ias   = iag - ibas(iSyma) !address within that sym
                     ibs   = ibg - ibas(iSyma)
                     xfd   = Fac(min(abs(ias-ibs),1)) !fac for diag
C
C     L(p,b,J) = sum_a  xfd* L(a,b,J) * C(p,a)
C     ----------------------------------------
                     DO jDen=kDen,nDen
                        ! pointer to C(1,a)
*                        ISMOSA = ipMOs(iSyma,jDen)
*     &                         + nPorb(iSyma,jDen)*(ias-1)
*                        ! pointer to C(1,b)
*                        ISMOSB = ipMOs(iSyma,jDen)
*     &                         + nPorb(iSyma,jDen)*(ibs-1)
* Replaced by:
                        ISYMB=ISYMA
                        NBA=NBAS(ISYMA)
                        ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA,JDEN)-1)
                        NBB=NBAS(ISYMB)
                        ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB,JDEN)-1)

*                        ichot = nPorb(iSyma,jDen)*NBAS(iSyma)*(JVEC-1)
                        ichot = nUse(iSyma,jDen)*NBA*(JVEC-1)
     &                        + ipChoT(iSyma,jDen)

*                        kchot = ichot + nPorb(iSyma,jDen)*(ias-1)
                        kchot = ichot + nUse(iSyma,jDen)*(ias-1)

*                        CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
*     &                             Work(ISMOSB),1,Work(kchot),1)
                        CALL DAXPY_(nUse(iSyma,jDen),xfd*Scr(kscr),
     &                             CMO(ISCB),NBB,Work(kchot),1)

*                        kchot = ichot + nPorb(iSyma,jDen)*(ibs-1)
                        kchot = ichot + nUse(iSyma,jDen)*(ibs-1)

*                        CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
*     &                             Work(ISMOSA),1,Work(kchot),1)
                        CALL DAXPY_(nUse(iSyma,jDen),xfd*Scr(kscr),
     &                             CMO(ISCA),NBA,Work(kchot),1)

                     END DO  ! loop over densities

                  ENDIF  ! skipping blocks check

               End Do  ! jRab loop

            ELSE  ! jSym.ne.1

               Do jRab=1,nnBstR(jSym,iLoc)

                  kRab = iiBstr(jSym,iLoc) + jRab
                  iRab = IndRed(kRab,iLoc)

                  iag   = iRS2F(1,iRab)  !global address
                  ibg   = iRS2F(2,iRab)

                  iSyma = cho_isao(iag)  !symmetry block
                  iSymb = muld2h(jSym,iSyma)

                  kscr  = kscr + 1

                  ias   = iag - ibas(iSyma)  !address within that sym
                  ibs   = ibg - ibas(iSymb)

                  IF (iSkip(iSyma).ne.0) THEN
C
C     L(p,b,J) = sum_a  L(a,b,J) * C(p,a)
C     -----------------------------------
                     DO jDen=kDen,nDen

*                        ISMOSA = ipMOs(iSyma,jDen)
*     &                         + nPorb(iSyma,jDen)*(ias-1)
* Replaced by:
                        NBA=NBAS(ISYMA)
                        ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA,JDEN)-1)
                        NBB=NBAS(ISYMB)

*                        kchot = nPorb(iSyma,jDen)*NBAS(iSymb)*(JVEC-1)
*     &                        + nPorb(iSyma,jDen)*(ibs-1)
*     &                        + ipChoT(iSyma,jDen)
                        kchot = nUse(iSyma,jDen)*NBB*(JVEC-1)
     &                        + nUse(iSyma,jDen)*(ibs-1)
     &                        + ipChoT(iSyma,jDen)

*                        CALL DAXPY_(nPorb(iSyma,jDen),Scr(kscr),
*     &                           Work(ISMOSA),1,Work(kchot),1)
                        CALL DAXPY_(nUse(iSyma,jDen),Scr(kscr),
     &                           CMO(ISCA),NBA,Work(kchot),1)

                     END DO

                  ENDIF  ! skipping block

                  IF (iSkip(iSymb).ne.0) THEN
C
C     L(p,a,J) = sum_b  L(a,b,J) * C(p,b)
C     -----------------------------------
                     DO jDen=kDen,nDen

*                        ISMOSB = ipMOs(iSymb,jDen)
*     &                         + nPorb(iSymb,jDen)*(ibs-1)
* Replaced by:
                        NBA=NBAS(ISYMA)
                        NBB=NBAS(ISYMB)
                        ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB,JDEN)-1)

*                        kchot = nPorb(iSymb,jDen)*NBAS(isyma)*(JVEC-1)
*     &                        + nPorb(iSymb,jDen)*(ias-1)
*     &                        + ipChoT(iSymb,jDen)
                        kchot = nUse(iSymb,jDen)*NBA*(JVEC-1)
     &                        + nUse(iSymb,jDen)*(ias-1)
     &                        + ipChoT(iSymb,jDen)

*                        CALL DAXPY_(nPorb(iSymb,jDen),Scr(kscr),
*     &                           Work(ISMOSB),1,Work(kchot),1)
                        CALL DAXPY_(nUse(iSymb,jDen),Scr(kscr),
     &                           CMO(ISCB),NBB,Work(kchot),1)

                     END DO

                  ENDIF  ! skipping blocks check

               End Do  ! jRab loop

            ENDIF ! total symmetric vectors check

         END DO

      ELSEIF (iSwap.eq.1) THEN     ! L(ap,J) storage

      NREAD = 0
      kchot = 0
      DO JVEC=1,JNUM

         JRED = InfVec(JVEC1-1+JVEC,2,jSym)

         IF (JRED .NE. IREDC) THEN ! JRED is not the reduced set in core
            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            IREDC = JRED
         END IF

         kscr = NREAD
         NREAD= NREAD + nDimRS(jSym,JRED)

        IF (JSYM.eq.1) THEN  ! L(a,b,J)=L(b,a,J); only a.ge.b present

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)=Sym(p)

            kscr  = kscr + 1

            IF (iSkip(iSyma).ne.0) THEN

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSyma)
            xfd   = Fac(min(abs(ias-ibs),1)) !fac for diagonal elements
C
C     L(a,p,J) = sum_b  xfd* L(a,b,J) * C(p,b)
C     ----------------------------------------
            DO jDen=kDen,nDen

*               ISMOSA = ipMOs(iSyma,jDen)
*     &                + nPorb(iSyma,jDen)*(ias-1)  ! pointer to C(1,a)

*               ISMOSB = ipMOs(iSyma,jDen)
*     &                + nPorb(iSyma,jDen)*(ibs-1)  ! pointer to C(1,b)
* Replaced by:
                ISYMB=ISYMA
                NBA=NBAS(ISYMA)
                ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA,JDEN)-1)
                NBB=NBAS(ISYMB)
                ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB,JDEN)-1)

*               kchot = NBAS(iSyma)*nPorb(iSyma,jDen)*(JVEC-1)
               kchot = NBA*nUse(iSyma,jDen)*(JVEC-1)
     &               + ipChoT(iSyma,jDen) - 1

*               CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),Work(ISMOSB),
*     &                    1,Work(kchot+ias),NBAS(iSyma))
               CALL DAXPY_(nUse(iSyma,jDen),xfd*Scr(kscr),CMO(ISCB),
     &                    NBB,Work(kchot+ias),NBA)

*               CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),Work(ISMOSA),
*     &                    1,Work(kchot+ibs),NBAS(iSyma))
               CALL DAXPY_(nUse(iSyma,jDen),xfd*Scr(kscr),CMO(ISCA),
     &                    NBA,Work(kchot+ibs),NBA)

            END DO  ! loop over densities

            ENDIF  ! skipping blocks check

         End Do  ! jRab loop

        ELSE  ! jSym.ne.1

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block
            iSymb = muld2h(jSym,iSyma)

            kscr  = kscr + 1

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSymb)

            IF (iSkip(iSyma).ne.0) THEN
C
C     L(a,q,J) = sum_b  L(a,b,J) * C(q,b)
C     -----------------------------------
             DO jDen=kDen,nDen

*               ISMOSB = ipMOs(iSymb,jDen)
*     &                + nPorb(iSymb,jDen)*(ibs-1)  ! pointer to C(1,b)
* Replaced by:
               NBA=NBAS(ISYMA)
               NBB=NBAS(ISYMB)
               ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB,JDEN)-1)

*               kchot = NBAS(iSyma)*nPorb(iSymb,jDen)*(JVEC-1)
               kchot = NBA*nUse(iSymb,jDen)*(JVEC-1)
     &               + ias
     &               + ipChoT(iSyma,jDen) - 1

*               CALL DAXPY_(nPorb(iSymb,jDen),Scr(kscr),Work(ISMOSB),1,
               CALL DAXPY_(nUse(iSymb,jDen),Scr(kscr),CMO(ISCB),NBB,
     &                    Work(kchot),NBA)

             END DO

            ENDIF  ! skipping block

            IF (iSkip(iSymb).ne.0) THEN
C
C     L(b,q,J) = sum_a  L(a,b,J) * C(q,a)
C     -----------------------------------
             DO jDen=kDen,nDen

*               ISMOSA = ipMOs(iSyma,jDen)
*     &                + nPorb(iSyma,jDen)*(ias-1)  ! pointer to C(1,a)
* Replaced by:
                NBA=NBAS(ISYMA)
                ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA,JDEN)-1)
                NBB=NBAS(ISYMB)

*               kchot = NBAS(iSymb)*nPorb(iSyma,jDen)*(JVEC-1)
               kchot = NBB*nUse(iSyma,jDen)*(JVEC-1)
     &               + ibs
     &               + ipChoT(iSymb,jDen) - 1

*               CALL DAXPY_(nPorb(iSyma,jDen),Scr(kscr),Work(ISMOSA),1,
               CALL DAXPY_(nUse(iSyma,jDen),Scr(kscr),CMO(ISCA),NBA,
     &                    Work(kchot),NBA)

             END DO

            ENDIF  ! skipping blocks check

         End Do  ! jRab loop

        ENDIF ! total symmetric vectors check

      END DO


      ELSEIF (iSwap.eq.2) THEN     ! L(pJ,b) storage

             NREAD = 0
             kchot = 0
             DO JVEC=1,JNUM

                JRED = InfVec(JVEC1-1+JVEC,2,jSym)

                IF (JRED .NE. IREDC) THEN
                   Call Cho_X_SetRed(irc,iLoc,JRED)
                   IREDC = JRED
                END IF

             kscr = NREAD
             NREAD= NREAD + nDimRS(jSym,JRED)

             IF (JSYM.eq.1) THEN  ! L(a,b,J)=L(b,a,J); only a.ge.b

                Do jRab=1,nnBstR(jSym,iLoc)

                   kRab = iiBstr(jSym,iLoc) + jRab
                   iRab = IndRed(kRab,iLoc)

                   iag   = iRS2F(1,iRab)  !global address
                   ibg   = iRS2F(2,iRab)

                   iSyma = cho_isao(iag)

                   kscr  = kscr + 1

                   IF (iSkip(iSyma).ne.0) THEN

                      ias   = iag - ibas(iSyma)
                      ibs   = ibg - ibas(iSyma)
                      xfd   = Fac(min(abs(ias-ibs),1))
C
C     L(p,J,b) = sum_a  xfd* L(a,b,J) * C(p,a)
C     ----------------------------------------
                      DO jDen=kDen,nDen

*                         ISMOSA = ipMOs(iSyma,jDen)
*     &                          + nPorb(iSyma,jDen)*(ias-1)

*                         ISMOSB = ipMOs(iSyma,jDen)
*     &                          + nPorb(iSyma,jDen)*(ibs-1)
* Replaced by:
                        ISYMB=ISYMA
                        NBA=NBAS(ISYMA)
                        ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA,JDEN)-1)
                        NBB=NBAS(ISYMB)
                        ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB,JDEN)-1)

*                         ichot = nPorb(iSyma,jDen)*(jVref+JVEC-2)
                         ichot = nUse(iSyma,jDen)*(jVref+JVEC-2)
     &                         + ipChoT(iSyma,jDen)

*                         kchot = ichot + nPorb(iSyma,jDen)*NUMV*(ias-1)
                         kchot = ichot + nUse(iSyma,jDen)*NUMV*(ias-1)

*                         CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
*     &                           Work(ISMOSB),1,Work(kchot),1)
                         CALL DAXPY_(nUse(iSyma,jDen),xfd*Scr(kscr),
     &                           CMO(ISCB),NBB,Work(kchot),1)

*                         kchot = ichot + nPorb(iSyma,jDen)*NUMV*(ibs-1)
                         kchot = ichot + nUse(iSyma,jDen)*NUMV*(ibs-1)

*                         CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),
*     &                           Work(ISMOSA),1,Work(kchot),1)
                         CALL DAXPY_(nUse(iSyma,jDen),xfd*Scr(kscr),
     &                           CMO(ISCA),NBA,Work(kchot),1)

                      END DO  ! loop over densities

                   ENDIF  ! skipping blocks check

                End Do  ! jRab loop

             ELSE  ! jSym.ne.1

                Do jRab=1,nnBstR(jSym,iLoc)

                   kRab = iiBstr(jSym,iLoc) + jRab
                   iRab = IndRed(kRab,iLoc)

                   iag   = iRS2F(1,iRab)  !global address
                   ibg   = iRS2F(2,iRab)

                   iSyma = cho_isao(iag)  !symmetry block
                   iSymb = muld2h(jSym,iSyma)  !(syma>symb)

                   kscr  = kscr + 1

                   ias   = iag - ibas(iSyma)
                   ibs   = ibg - ibas(iSymb)

                   IF (iSkip(iSyma).ne.0) THEN
C
C     L(p,J,b) = sum_a  L(a,b,J) * C(p,a)
C     -----------------------------------
                      DO jDen=kDen,nDen

*                         ISMOSA = ipMOs(iSyma,jDen)
*     &                          + nPorb(iSyma,jDen)*(ias-1)
* Replaced by:
                        NBA=NBAS(ISYMA)
                        ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA,JDEN)-1)

*                         kchot = nPorb(iSyma,jDen)*NUMV*(ibs-1)
*     &                         + nPorb(iSyma,jDen)*(jVref+JVEC-2)
*     &                         + ipChoT(iSyma,jDen)
                        kchot = nUse(iSyma,jDen)*NUMV*(ibs-1)
     &                         + nUse(iSyma,jDen)*(jVref+JVEC-2)
     &                         + ipChoT(iSyma,jDen)

*                         CALL DAXPY_(nPorb(iSyma,jDen),Scr(kscr),
*     &                              Work(ISMOSA),1,Work(kchot),1)
                        CALL DAXPY_(nUse(iSyma,jDen),Scr(kscr),
     &                              CMO(ISCA),NBA,Work(kchot),1)

                      END DO

                   ENDIF  ! skipping block

                   IF (iSkip(iSymb).ne.0) THEN
C
C     L(p,J,a) = sum_b  L(a,b,J) * C(p,b)
C     -----------------------------------
                      DO jDen=kDen,nDen

*                         ISMOSB = ipMOs(iSymb,jDen)
*     &                          + nPorb(iSymb,jDen)*(ibs-1)
* Replaced by:
                        NBB=NBAS(ISYMB)
                        ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB,JDEN)-1)

*                         kchot = nPorb(iSymb,jDen)*NUMV*(ias-1)
*     &                         + nPorb(iSymb,jDen)*(jVref+JVEC-2)
*     &                         + ipChoT(iSymb,jDen)
                         kchot = nUse(iSymb,jDen)*NUMV*(ias-1)
     &                         + nUse(iSymb,jDen)*(jVref+JVEC-2)
     &                         + ipChoT(iSymb,jDen)

*                         CALL DAXPY_(nPorb(iSymb,jDen),Scr(kscr),
*     &                              Work(ISMOSB),1,Work(kchot),1)
                         CALL DAXPY_(nUse(iSymb,jDen),Scr(kscr),
     &                              CMO(ISCB),NBB,Work(kchot),1)

                      END DO

                   ENDIF  ! skipping blocks check

                End Do  ! jRab loop

              ENDIF ! total symmetric vectors check

           END DO

      ELSEIF (iSwap.eq.3) THEN     ! L(aJ,p) storage

      NREAD = 0
      kchot = 0
      DO JVEC=1,JNUM

         JRED = InfVec(JVEC1-1+JVEC,2,jSym)

         IF (JRED .NE. IREDC) THEN ! JRED is not the reduced set in core
            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            IREDC = JRED
         END IF

         kscr = NREAD
         NREAD= NREAD + nDimRS(jSym,JRED)

        IF (JSYM.eq.1) THEN  ! L(a,b,J)=L(b,a,J); only a.ge.b present

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)=Sym(p)

            kscr  = kscr + 1

            IF (iSkip(iSyma).ne.0) THEN

            ias   = iag - ibas(iSyma)  !address within that sym block
            ibs   = ibg - ibas(iSyma)
            xfd   = Fac(min(abs(ias-ibs),1)) !scale fac for diag
C
C     L(a,J,p) = sum_b  xfd* L(a,b,J) * C(p,b)
C     ----------------------------------------
            DO jDen=kDen,nDen

*               ISMOSA = ipMOs(iSyma,jDen)
*     &                + nPorb(iSyma,jDen)*(ias-1)  ! pointer to C(1,a)

*               ISMOSB = ipMOs(iSyma,jDen)
*     &                + nPorb(iSyma,jDen)*(ibs-1)  ! pointer to C(1,b)
* Replaced by:
                ISYMB=ISYMA
                NBA=NBAS(ISYMA)
                ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA,JDEN)-1)
                NBB=NBAS(ISYMB)
                ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB,JDEN)-1)

               kchot = NBA*NUMV
     &               + NBA*(jVref+JVEC-2)
     &               + ipChoT(iSyma,jDen) - 1

*               CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),Work(ISMOSB),
*     &                    1,Work(kchot+ias),NBAS(iSyma)*NUMV)
               CALL DAXPY_(nUse(iSyma,jDen),xfd*Scr(kscr),CMO(ISCB),
     &                    NBB,Work(kchot+ias),NBA*NUMV)

*               CALL DAXPY_(nPorb(iSyma,jDen),xfd*Scr(kscr),Work(ISMOSA),
*     &                    1,Work(kchot+ibs),NBAS(iSyma))
               CALL DAXPY_(nUse(iSyma,jDen),xfd*Scr(kscr),CMO(ISCA),
     &                    NBA,Work(kchot+ibs),NBA)

            END DO  ! loop over densities

            ENDIF  ! skipping blocks check

         End Do  ! jRab loop

        ELSE  ! jSym.ne.1

         Do jRab=1,nnBstR(jSym,iLoc)

            kRab = iiBstr(jSym,iLoc) + jRab
            iRab = IndRed(kRab,iLoc)

            iag   = iRS2F(1,iRab)  !global address
            ibg   = iRS2F(2,iRab)

            iSyma = cho_isao(iag)  !symmetry block
            iSymb = muld2h(jSym,iSyma)

            kscr  = kscr + 1

            ias   = iag - ibas(iSyma)  !address within that symm block
            ibs   = ibg - ibas(iSymb)

            IF (iSkip(iSyma).ne.0) THEN
C
C     L(a,J,q) = sum_b  L(a,b,J) * C(q,b)
C     -----------------------------------
             DO jDen=kDen,nDen

*               ISMOSB = ipMOs(iSymb,jDen)
*     &                + nPorb(iSymb,jDen)*(ibs-1)  ! pointer to C(1,b)
* Replaced by:
               NBB=NBAS(ISYMA)
               NBB=NBAS(ISYMB)
               ISCB=IOFFC(ISYMB)+IBS+NBB*(ISTART(ISYMB,JDEN)-1)

               kchot = NBA*NUMV
     &               + NBA*(jVref+JVEC-2) + ias
     &               + ipChoT(iSyma,jDen) - 1

*               CALL DAXPY_(nPorb(iSymb,jDen),Scr(kscr),Work(ISMOSB),1,
               CALL DAXPY_(nUse(iSymb,jDen),Scr(kscr),CMO(ISCB),NBB,
     &                    Work(kchot),NBA*NUMV)

             END DO

            ENDIF  ! skipping block

            IF (iSkip(iSymb).ne.0) THEN
C
C     L(b,J,q) = sum_a  L(a,b,J) * C(q,a)
C     -----------------------------------
             DO jDen=kDen,nDen

*               ISMOSA = ipMOs(iSyma,jDen)
*     &                + nPorb(iSyma,jDen)*(ias-1)  ! pointer to C(1,a)
* Replaced by:
                NBA=NBAS(ISYMA)
                NBB=NBAS(ISYMB)
                ISCA=IOFFC(ISYMA)+IAS+NBA*(ISTART(ISYMA,JDEN)-1)

               kchot = NBB*NUMV
     &               + NBB*(jVref+JVEC-2) + ibs
     &               + ipChoT(iSymb,jDen) - 1

*               CALL DAXPY_(nPorb(iSyma,jDen),Scr(kscr),Work(ISMOSA),1,
               CALL DAXPY_(nUse(iSyma,jDen),Scr(kscr),CMO(ISCA),NBA,
     &                    Work(kchot),NBA*NUMV)

             END DO

            ENDIF  ! skipping blocks check

         End Do  ! jRab loop

        ENDIF ! total symmetric vectors check

      END DO


      ELSE   ! iSwap check

         write(6,*)SECNAM//': invalid argument. Iswap= ',Iswap
         irc = 66
         Return

      ENDIF  ! iSwap check


      irc=0

      Call qExit(SECNAM)

      Return
      END

**************************************************************
