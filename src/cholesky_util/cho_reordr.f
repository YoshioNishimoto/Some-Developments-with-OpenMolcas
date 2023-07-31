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
************************************************************************
      SUBROUTINE CHO_REORDR(irc,scr,lscr,jVref,JVEC1,JNUM,NUMV,JSYM,
     &                      IREDC,iSwap,ipChoV,iSkip)
************************************************************
*   Author: F. Aquilante
*
*   Purpose:  SCR(lscr) contains JNUM cholesky vectors
*             starting from JVEC1 and stored in reduced
*             sets. The routine performs a reallocation
*             of these elements in a set of target
*             arrays identified by the pointers ipChoV.
*             In the target arrays, the vectors are
*             stored in full dimension and as a
*             subset of a a given NUMV number of vectors.
*             Each pointer should thereby point to a
*             location where the corresponding Cholesky
*             vector of a given unique symmetry pair
*             of indices has to be stored
*
*   Input:
*       Ivec1 =  first vector to be copied
*       JNum  =  # of vectors to be copied
*
*       NumV  =  total # of vectors in the target arrays
*
*       iSwap :   = 0   L(a,b,J) is returned
*                       (in LT-storage if sym(a)=sym(b))
*                 = 1   L(a,J,b) is returned
*                 = 2   L(a,J,b) is returned
*                       (in SQ-storage if sym(a)=sym(b))
*
*       iSkip(syma)=0 : skip the symmetry block a.
*                    Any vector L(ab) or L(ba) with syma x symb=JSYM
*                    won't be returned in the target array
*
*       IREDC :  reduced set in core at the moment of
*                the call to the routine.
*                Can be set to -1 by the calling routine
*
*********************************************************
      use ChoArr, only: nDimRS, iRS2F
      use ChoSwp, only: InfVec, IndRed
      Implicit Real*8 (a-h,o-z)
      Real*8  Scr(lscr)
      Integer ipChoV(*),iSkip(*)

#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      Integer  cho_isao
      External cho_isao

************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
************************************************************************

**********************************************************
C
C    From Reduced sets to full storage
C    ---------------------------------
C
C    iSwap = 0
C
C     L{a,b,J} ---> L(a,b,J)
C
C    iSwap = 1
C
C     L{a,b,J} ---> L(a,J,b)
C
C    iSwap = 2
C
C     L{a,b,J} ---> L(ab,J)  ! with squaring of the
C                            ! "diagonal" symmetry blocks
C
**********************************************************

      iLoc = 3 ! use scratch location in reduced index arrays


      IF (jSym.eq.1 .and. iSwap.eq.0) Then  ! L(ab),J

         NREAD = 0
         kchov = 0
         kchov2= 0
         DO JVEC=1,JNUM

            JRED = InfVec(JVEC1-1+JVEC,2,jSym)

            IF (JRED .NE. IREDC) THEN ! JRED is not the rs in core
               Call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
               IREDC = JRED
            END IF

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,JRED)

            Do jRab=1,nnBstR(jSym,iLoc)

               kRab = iiBstr(jSym,iLoc) + jRab
               iRab = IndRed(kRab,iLoc)

               iag   = iRS2F(1,iRab)  !global address
               ibg   = iRS2F(2,iRab)

               iSyma = cho_isao(iag)  !symmetry block

               kscr  = kscr + 1

               If(iSkip(iSyma).ne.0)Then

                ias   = iag - ibas(iSyma)  !addr within that symm block
                ibs   = ibg - ibas(iSyma)

                iabf  = iTri(ias,ibs)
                kchov = (JVEC-1)*nBas(iSyma)*(nBas(iSyma)+1)/2 + iabf
     &                + ipChoV(iSyma) - 1

                Work(kchov) = Scr(kscr)

               EndIf

            End Do

         END DO


      ELSEIF (jSym.eq.1 .and. iSwap.eq.1) Then  ! LaJ,b

         NREAD = 0
         kchov = 0
         kchov2= 0
         DO JVEC=1,JNUM

            JRED = InfVec(JVEC1-1+JVEC,2,jSym)

            IF (JRED .NE. IREDC) THEN ! JRED is not the rs in core
               Call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
               IREDC = JRED
            END IF

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,JRED)

            Do jRab=1,nnBstR(jSym,iLoc)

               kRab = iiBstr(jSym,iLoc) + jRab
               iRab = IndRed(kRab,iLoc)

               iag   = iRS2F(1,iRab)  !global address
               ibg   = iRS2F(2,iRab)

               iSyma = cho_isao(iag)  !symmetry block

               kscr  = kscr + 1

               If(iSkip(iSyma).ne.0)Then

                ias   = iag - ibas(iSyma)  !addr within that symm block
                ibs   = ibg - ibas(iSyma)

                kchov1 = nBas(iSyma)*NUMV*(ibs-1)
     &                 + nBas(iSyma)*(jVref+JVEC-2) + ias
     &                 + ipChoV(iSyma) - 1
                kchov2 = nBas(iSyma)*NUMV*(ias-1)
     &                 + nBas(iSyma)*(jVref+JVEC-2) + ibs
     &                 + ipChoV(iSyma) - 1

                Work(kchov1) = Scr(kscr)
                Work(kchov2) = Scr(kscr)

               EndIf

            End Do

         END DO


      ELSEIF (jSym.eq.1 .and. iSwap.eq.2) Then  ! L[ab],J

         NREAD = 0
         kchov = 0
         kchov2= 0
         DO JVEC=1,JNUM

            JRED = InfVec(JVEC1-1+JVEC,2,jSym)

            IF (JRED .NE. IREDC) THEN ! JRED is not the rs in core
               Call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
               IREDC = JRED
            END IF

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,JRED)

            Do jRab=1,nnBstR(jSym,iLoc)

               kRab = iiBstr(jSym,iLoc) + jRab
               iRab = IndRed(kRab,iLoc)

               iag   = iRS2F(1,iRab)  !global address
               ibg   = iRS2F(2,iRab)

               iSyma = cho_isao(iag)  !symmetry block

               kscr  = kscr + 1

               If(iSkip(iSyma).ne.0)Then

                ias   = iag - ibas(iSyma)  !addr within that symm block
                ibs   = ibg - ibas(iSyma)

                kchov1 = nBas(iSyma)*nBas(iSyma)*(JVEC-1)
     &                 + nBas(iSyma)*(ibs-1) + ias
     &                 + ipChoV(iSyma) - 1
                kchov2 = nBas(iSyma)*nBas(iSyma)*(JVEC-1)
     &                 + nBas(iSyma)*(ias-1) + ibs
     &                 + ipChoV(iSyma) - 1

                Work(kchov1) = Scr(kscr)
                Work(kchov2) = Scr(kscr)

               EndIf

            End Do

         END DO


      ELSEIF (jSym.gt.1 .and. iSwap.eq.0) Then

         NREAD = 0
         kchov = 0
         kchov2= 0
         DO JVEC=1,JNUM

            JRED = InfVec(JVEC1-1+JVEC,2,jSym)

            IF (JRED .NE. IREDC) THEN ! JRED is not the rs in core
               Call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
               IREDC = JRED
            END IF

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,JRED)

            Do jRab=1,nnBstR(jSym,iLoc)

               kRab = iiBstr(jSym,iLoc) + jRab
               iRab = IndRed(kRab,iLoc)

               iag   = iRS2F(1,iRab)  !global address
               ibg   = iRS2F(2,iRab)

               iSyma = cho_isao(iag)  !symmetry block
               iSymb = MulD2h(jSym,iSyma) ! sym(a) > sym(b)

               kscr  = kscr + 1

               If(iSkip(iSyma).ne.0)Then

                ias   = iag - ibas(iSyma)  !addr within that symm block
                ibs   = ibg - ibas(iSymb)

                kchov = nBas(iSyma)*nBas(iSymb)*(JVEC-1)
     &                + nBas(iSyma)*(ibs-1) + ias
     &                + ipChoV(iSyma) - 1

                Work(kchov) = Scr(kscr)

               EndIf

            End Do

         END DO


      ELSEIF (jSym.gt.1 .and. iSwap.eq.1) Then

         NREAD = 0
         kchov = 0
         kchov2= 0
         DO JVEC=1,JNUM

            JRED = InfVec(JVEC1-1+JVEC,2,jSym)

            IF (JRED .NE. IREDC) THEN ! JRED is not the rs in core
               Call Cho_X_SetRed(irc,iLoc,JRED) !set indx arrays at iLoc
               IREDC = JRED
            END IF

            kscr = NREAD
            NREAD= NREAD + nDimRS(jSym,JRED)

            Do jRab=1,nnBstR(jSym,iLoc)

               kRab = iiBstr(jSym,iLoc) + jRab
               iRab = IndRed(kRab,iLoc)

               iag   = iRS2F(1,iRab)  !global address
               ibg   = iRS2F(2,iRab)

               iSyma = cho_isao(iag)  !symmetry block
               iSymb = MulD2h(jSym,iSyma) ! sym(a) > sym(b)

               kscr  = kscr + 1

               If(iSkip(iSyma).ne.0)Then

                ias   = iag - ibas(iSyma)  !addr within that symm block
                ibs   = ibg - ibas(iSymb)

                kchov = nBas(iSyma)*NUMV*(ibs-1)
     &                + nBas(iSyma)*(jVref+JVEC-2) + ias
     &                + ipChoV(iSyma) - 1

                Work(kchov) = Scr(kscr)

               EndIf

            End Do

         END DO


      ELSE

       write(6,*)'Wrong parameters combination. JSYM,iSwap= ',JSYM,iSwap
       irc = 66
       Return

      ENDIF


      irc=0

      Return
      END

**************************************************************
