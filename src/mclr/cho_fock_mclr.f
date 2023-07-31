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
* Copyright (C) Mickael G. Delcey                                      *
************************************************************************
      SUBROUTINE CHO_Fock_MCLR(DA,G2,W_JA,W_KA,W_FkA,CVa,W_CMO,nIsh,
     &                         nAsh,LuAChoVec)
************************************************************************
*                                                                      *
*  Author : M. G. Delcey                                               *
*                                                                      *
************************************************************************
      use ChoArr, only: nBasSh, nDimRS
      use ChoSwp, only: InfVec
      use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
      Implicit Real*8 (a-h,o-z)
#include "warnings.h"
      Character(LEN=13), Parameter :: SECNAM = 'CHO_FOCK_MCLR'
      Integer   ipLpq(8,3)
      Integer   LuAChoVec(8)
      Integer   nAsh(8),nIsh(8)
#include "cholesky.fh"
#include "choorb.fh"
#include "real.fh"
#include "stdalloc.fh"
      Type (DSBA_type) CVa, JA(1), KA, Fka, CMO, Scr
      Real*8 DA(*), G2(*), W_JA(*), W_KA(*), W_FkA(*), W_CMO(*)
      Real*8, parameter:: xone=-One, FactCI = -Two, FactXI = Half
      Integer , External :: Cho_LK_MaxVecPerBatch
      Integer, Allocatable:: kOffSh(:,:)
      Real*8, Allocatable::  Fab(:), Lrs(:), LF(:)
      Logical add
************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
************************************************************************
*
      Call Allocate_DT(JA(1),nBas,nBas,nSym,aCase='TRI',Ref=W_JA )
      Call Allocate_DT(KA ,nBas,nBas,nSym,            Ref=W_KA )
      Call Allocate_DT(FkA,nBas,nBas,nSym,            Ref=W_FkA)
      Call Allocate_DT(CMO,nBas,nBas,nSym,            Ref=W_CMO)
*
**    Compute Shell Offsets ( MOs and transformed vectors)
*
      Call mma_allocate(kOffSh,nShell,nSym,Label='kOffSh')
      Do iSyma=1,nSym
         LKsh=0
         Do iaSh=1,nShell    ! kOffSh(iSh,iSym)
            kOffSh(iaSh,iSyma) = LKsh
            LKsh = LKsh + nBasSh(iSyma,iaSh)
         End Do
      End Do
*
*     memory for the Q matrices --- temporary array
      Call Allocate_DT(Scr,nBas,nBas,nSym)
      Scr%A0(:)=Zero
*
      MaxVecPerBatch=Cho_LK_MaxVecPerBatch()
*
**    Start looping!
*
      DO jSym=1,nSym
        NumCV=NumCho(jSym)
        Call GAIGOP_SCAL(NumCV,'max')
        If (NumCV .lt. 1) GOTO 1000
*
        iLoc = 3 ! use scratch location in reduced index arrays
*
**    Estimate memory need
*
        mTvec = 0
        do l=1,nSym
           k=Muld2h(l,JSYM)
           mTvec = mTvec + nAsh(k)*nBas(l)*3
        end do
*
        JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
        JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec

c --- entire red sets range for parallel run
        Call GAIGOP_SCAL(JRED1,'min')
        Call GAIGOP_SCAL(JRED2,'max')

        Do JRED=JRED1,JRED2
          CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)
          If (nVrs.eq.0) GOTO 999
          if (nVrs.lt.0) then
            Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
     &                            ,nVrs
            call Abend
          endif
          Call Cho_X_SetRed(irc,iLoc,JRED)
c         !set index arrays at iLoc
          if(irc.ne.0)then
            Write(6,*)SECNAM//'cho_X_setred non-zero return code.'//
     &                        ' rc= ',irc
            call Abend
          endif
          IREDC=JRED

          nRS = nDimRS(JSYM,JRED)
          If (jSym.eq.1) Then
            Call mma_allocate(Fab,nRS,Label='Fab')
            Fab(:)=Zero
          EndIf

          Call mma_MaxDBLE(LWORK)
          nVec = min(LWORK/(nRS+mTvec+1),min(nVrs,MaxVecPerBatch))
          If (nVec.lt.1) Then
             WRITE(6,*) SECNAM//': Insufficient memory for batch'
             WRITE(6,*) 'LWORK= ',LWORK
             WRITE(6,*) 'min. mem. need= ',nRS+mTvec+1
             WRITE(6,*) 'nRS= ',nRS
             WRITE(6,*) 'mTvec= ',mTvec
             WRITE(6,*) 'jsym= ',jsym
             CALL Quit(_RC_MEMORY_ERROR_)
             nBatch = -9999  ! dummy assignment
          End If
          LREAD = nRS*nVec

          Call mma_allocate(Lrs,LREAD,Label='Lrs')
          CALL mma_allocate(LF,mTvec*nVec,Label='LF')

          nBatch = (nVrs-1)/nVec + 1

          DO iBatch=1,nBatch

            If (iBatch.eq.nBatch) Then
               JNUM = nVrs - nVec*(nBatch-1)
            else
               JNUM = nVec
            endif
************************************************************************
*                                                                      *
*           START WORKING                                              *
*                                                                      *
************************************************************************
*
**          Read Cholesky vector
*
            JVEC = nVec*(iBatch-1) + iVrs
            IVEC2 = JVEC - 1 + JNUM

            CALL CHO_VECRD(Lrs,LREAD,JVEC,IVEC2,JSYM,
     &                     NUMV,IREDC,MUSED)

            If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
               RETURN
            End If
*
*            CALL CWTIME(TCINT1,TWINT1)

C --- Set up the skipping flags and the pointers ipLpq
C --- The memory used before for the full-dimension AO-vectors
C ---     is now re-used to store half and full transformed
C ---     vectors in the active space
C -------------------------------------------------------------
             lChoa=0
             Do i=1,nSym

                k = Muld2h(i,JSYM)

                ipLpq(k,1) = 1 + lChoa   ! Lvb,J
                ipLpq(k,2) = ipLpq(k,1)     ! Lvi,J i general MO index
     &                     + nAsh(k)*nBas(i)*JNUM
                ipLpq(k,3) = ipLpq(k,2)     ! L~vi,J ~ transformed index
     &                       + nAsh(k)*nBas(i)*JNUM

                lChoa= lChoa + nAsh(k)*nBas(i)*3*JNUM

             End Do
*
**  Read half-transformed cho vectors
*
             ioff=0
             Do i=1,nSym
                k = Muld2h(i,JSYM)
                lvec=nAsh(k)*nBas(i)*JNUM
                iAdr=(JVEC-1)*nAsh(k)*nBas(i)+ioff
                call DDAFILE(LuAChoVec(Jsym),2,LF(ipLpq(k,1)),
     &                       lvec,iAdr)
                ioff=ioff+nAsh(k)*nBas(i)*NumCho(jSym)
             End Do
*             CALL CWTIME(TCINT2,TWINT2)
*             tint1(1) = tint1(1) + (TCINT2 - TCINT1)
*             tint1(2) = tint1(2) + (TWINT2 - TWINT1)
C --------------------------------------------------------------------
C --- Active-Active transformation  Lvw,J = sum_b  Lvb,J * C2(w,b)
C --------------------------------------------------------------------
             Do iSymb=1,nSym

               iSymv = MulD2h(JSYM,iSymb)
               NAv = nAsh(iSymv)
               NAw = nAsh(iSymb)

               If(NAv*NBAS(iSymb).ne.0)Then

*                CALL CWTIME(TCINT2,TWINT2)

                 Do JVC=1,JNUM
                  ipLvb = ipLpq(iSymv,1) + NAv*NBAS(iSymb)*(JVC-1)
                  ipLvw = ipLpq(iSymv,2) + NAv*Naw*(JVC-1)
                  CALL DGEMM_('N','T',NAv,Naw,NBAS(iSymb),
     &                       One,LF(ipLvb),NAv,
     &                       CVa%SB(iSymb)%A2,Naw,
     &                      Zero,LF(ipLvw),NAv)
                 End Do
*                 CALL CWTIME(TCINT2,TWINT2)
*                 tint1(1) = tint1(1) + (TCINT2 - TCINT3)
*                 tint1(2) = tint1(2) + (TWINT2 - TWINT3)
C
C
C ************ EVALUATION OF THE ACTIVE FOCK MATRIX *************
*Coulomb term
                 ipVJ = ipLpq(iSymv,3)
                 ipLvtw = ipLpq(iSymv,2)

                 CALL DGEMV_('T',Nav*Naw,JNUM,
     &                  ONE,LF(ipLvtw),Nav*Naw,
     &                  DA,1,ZERO,LF(ipVJ),1)
*
                 CALL DGEMV_('N',nRS,JNUM,
     &                -FactCI,Lrs,nRS,
     &                LF(ipVJ),1,1.0d0,Fab,1)

*                 CALL CWTIME(TCINT2,TWINT2)
*                 tact(1) = tact(1) + (TCINT2 - TCINT3)
*                 tact(2) = tact(2) + (TWINT2 - TWINT3)
C --------------------------------------------------------------------
C --- Formation of the Q matrix Qpx = L~py Lvw Gxyvw
C --------------------------------------------------------------------
                 Do JVC=1,JNUM
*Lxy=Lvw Gxyvw
*MGD probably additional nSym loop
                  ipLvb = ipLpq(iSymv,1) + NAv*NBAS(iSymb)*(JVC-1)
                  ipLvw = ipLpq(iSymv,2) + NAv*Naw*(JVC-1)
                  ipLxy = ipLpq(iSymv,3) + NAv*Naw*(JVC-1)
                  CALL DGEMV_('N',NAv*Naw,NAv*Naw,
     &               ONE,G2,NAv*Naw,
     &               LF(ipLvw),1,ZERO,LF(ipLxy),1)
*Qpx=Lpy Lxy
                  Call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,
     &                       One,LF(ipLvb),NAv,
     &                       LF(ipLxy),Naw,
     &                      ONE,Scr%SB(iSymb)%A2,NBAS(iSymb))
                 End Do
*                 CALL CWTIME(TCINT3,TWINT3)
*                 tQmat(1) = tQmat(1) + (TCINT3 - TCINT2)
*                 tQmat(2) = tQmat(2) + (TWINT3 - TWINT2)
C
C
C ************ EVALUATION OF THE ACTIVE FOCK MATRIX *************
*Exchange term
                 Do JVC=1,JNUM
                   ipLvb = ipLpq(iSymv,1)+ NAv*NBAS(iSymb)*(JVC-1)
                   ipLwb = ipLpq(iSymv,2)+ NAv*NBAS(iSymb)*(JVC-1)
                   Call DGEMM_('T','N',NBAS(iSymb),Nav,Nav,
     &                         ONE,LF(ipLvb),Nav,
     &                         DA,Nav,ZERO,
     &                         LF(ipLwb),NBAS(iSymb))
                 End Do
*                 CALL CWTIME(TCINT2,TWINT2)
*                 tact(1) = tact(1) + (TCINT2 - TCINT3)
*                 tact(2) = tact(2) + (TWINT2 - TWINT3)
                 Do JVC=1,JNUM
                   ipLwb = ipLpq(iSymv,2)+ NAv*NBAS(iSymb)*(JVC-1)
                   Do is=1,NBAS(iSymb)
                    ipLtvb = ipLpq(iSymv,1)+ NAv*NBAS(iSymb)*(JVC-1)
     &                      + Nav*(is-1)
                    CALL DGEMV_('N',NBAS(iSymb),Nav,
     &                   -FactXI,LF(ipLwb),NBAS(iSymb),
     &                   LF(ipLtvb),1,ONE,KA%SB(iSymb)%A2(:,is),1)

                  EndDo
                 End Do
*                 CALL CWTIME(TCINT3,TWINT3)
*                 tact(1) = tact(1) + (TCINT3 - TCINT2)
*                 tact(2) = tact(2) + (TWINT3 - TWINT2)
               End If
             End Do
*
**          All good things come to an end
*
          END DO  ! end batch loop

c --- backtransform fock matrix to full storage
          If(JSYM.eq.1)Then
             add = .True.
             nMat = 1
             Call swap_rs2full(irc,iLoc,nRS,nMat,JSYM,
     &                           JA,Fab,add)
             Call mma_deallocate(Fab)
          EndIf
          Call mma_deallocate(Lrs)
          Call mma_deallocate(LF)
 999      Continue
        End Do  ! loop over red sets
 1000   CONTINUE
      End Do    ! loop over JSYM
************************************************************************
*                                                                      *
*     POST PROCESSING                                                  *
*                                                                      *
************************************************************************
*
**    Accumulate Coulomb and Exchange contributions
*
      Do iSym=1,nSym

         Do iaSh=1,nShell
            ioffa = kOffSh(iaSh,iSym)

            Do ibSh=1,nShell

               ioffb = kOffSh(ibSh,iSym)

               Do ib=1,nBasSh(iSym,ibSh)
                Do ia=1,nBasSh(iSym,iaSh)
*MGD warning with sym

                  iag = ioffa + ia
                  ibg = ioffb + ib
                  iabg= iTri(iag,ibg)

                  FkA%SB(iSym)%A2(iag,ibg)
     &                    = JA(1)%SB(iSym)%A1(iabg)
     &                    + KA%SB(iSym)%A2(iag,ibg)
     &                    + KA%SB(iSym)%A2(ibg,iag)
                End Do
               End Do
            End Do
         End Do
      End Do

      Call Deallocate_DT(JA(1))
      Call Allocate_DT(JA(1) ,nBas,nBas,nSym,           Ref=W_JA )
*
**Transform Fock and Q matrix to MO basis
*
      Do iS=1,nSym
        jS=iS
        If (nBas(iS).ne.0) Then
          Call DGEMM_('T','N',nBas(jS),nBas(iS),nBas(iS),
     &                1.0d0,FkA%SB(iS)%A2,nBas(iS),
     &                      CMO%SB(iS)%A2,nBas(iS),
     &                0.0d0,JA(1)%SB(iS)%A2,nBas(jS))
          FkA%SB(is)%A2(:,:)=Zero
          Call DGEMM_('T','N',nBas(jS),nIsh(jS),nBas(iS),
     &                1.0d0,JA(1)%SB(iS)%A2,nBas(iS),
     &                      CMO%SB(jS)%A2,nBas(jS),
     &                0.0d0,FkA%SB(iS)%A2,nBas(jS))
*         ioff=nIsh(iS)+1
          iOff = 1 + nIsh(iS) * nBas(iS)
          Call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),
     &                 1.0d0,CMO%SB(iS)%A2,nBas(jS),
     &                       Scr%SB(iS)%A2,nBas(jS),
*    &                 0.0d0,FkA%SB(iS)%A2(1:,ioff),nBas(jS))
     &                 0.0d0,FkA%SB(iS)%A1(iOff:),nBas(jS))
        EndIf
      End Do
************************************************************************
*                                                                      *
*     TERMINATING                                                      *
*                                                                      *
************************************************************************
      Call deallocate_DT(Scr)
      Call mma_deallocate(kOffSh)
      Call Deallocate_DT(CMO)
      Call Deallocate_DT(FkA)
      Call Deallocate_DT(KA)
      Call Deallocate_DT(JA(1))

      Return
      END
