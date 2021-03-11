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
      SUBROUTINE CHO_Fock_MCLR(ipDA,ipG2,ipJA,ipKA,ipFkA,
     &                      ipAsh,ipCMO,nIsh,nAsh,LuAChoVec)

************************************************************************
*                                                                      *
*  Author : M. G. Delcey                                               *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "warnings.fh"
      Character*13 SECNAM
      Parameter (SECNAM = 'CHO_FOCK_MCLR')
      Integer   ISTLT(8),ISTSQ(8),ISSQ(8,8),kaOff(8),ipLpq(8,3)
      Integer   ipAsh(*),ipAorb(8,2),LuAChoVec(8)
      Integer   nAsh(8),nIsh(8)
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
      parameter ( N2 = InfVec_N2 )
      parameter (zero = 0.0D0, one = 1.0D0, xone=-1.0D0)
      parameter (FactCI = -2.0D0, FactXI = 0.5D0)
      Character*6 mode
      Integer   Cho_LK_MaxVecPerBatch
      External  Cho_LK_MaxVecPerBatch
************************************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
******
      iTri(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
******
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)
******
      NBASSH(I,J)=IWORK(ip_NBASSH-1+NSYM*(J-1)+I)
******
      kOffSh(i,j) = iWork(ip_kOffSh+nShell*(j-1)+i-1)
************************************************************************
*
      nDen=1
*
**    Compute offsets
*
      ISTLT(1)=0
      ISTSQ(1)=0
      kAOff(1)=0
      nnA=0
      nsBB=nBas(1)**2
      DO ISYM=2,NSYM
        ISTSQ(ISYM)=nsBB
        nsBB = nsBB + nBas(iSym)**2
        NBB=NBAS(ISYM-1)*(NBAS(ISYM-1)+1)/2
        ISTLT(ISYM)=ISTLT(ISYM-1)+NBB
        nnA = nnA + nAsh(iSym-1)
        kAOff(iSym)=nnA
      End Do
      nnA = nnA + nAsh(nSym)
*
      nnBSQ=0
      DO LSYM=1,NSYM
         DO KSYM=LSYM,NSYM
            ISSQ(KSYM,LSYM) = nnBSQ
            ISSQ(LSYM,KSYM) = nnBSQ ! symmetrization
            nnBSQ = nnBSQ + nBas(kSym)*nBas(lSym)
         END DO
      END DO
*
**    Compute Shell Offsets ( MOs and transformed vectors)
*
      Call GetMem('ip_kOffSh','Allo','Inte',ip_kOffSh,nShell*nSym)
      Do iSyma=1,nSym
         LKsh=0
         Do iaSh=1,nShell    ! kOffSh(iSh,iSym)
            iWork(ip_kOffSh+nShell*(iSyma-1)+iaSh-1) = LKsh
            LKsh = LKsh + nBasSh(iSyma,iaSh)
         End Do
      End Do
*
      DO jDen=1,nDen
         ipAorb(1,jDen)= ipAsh(jDen)
         DO ISYM=2,NSYM
            ipAorb(iSym,jDen) = ipAorb(iSym-1,jDen)
     &                        + nAsh(iSym-1)*nBas(iSym-1)
         END DO
      END DO
*     memory for the Q matrices --- temporary array
      Call GetMem('Qmat','ALLO','REAL',ipScr,nsBB*nDen)
      Call Fzero(Work(ipScr),nsBB*nDen)
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
        myJRED1=JRED1 ! first red set present on this node
        myJRED2=JRED2 ! last  red set present on this node

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
            Call GetMem('rsFC','Allo','Real',ipFab,nRS)
            Call Fzero(Work(ipFab),nRS)
          EndIf

          Call GetMem('MaxM','Max','Real',KDUM,LWORK)
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

          Call GetMem('rsL','Allo','Real',ipLrs,LREAD)
          CALL GetMem('FullV','Allo','Real',ipLF,mTvec*nVec)

          nBatch = (nVrs-1)/nVec + 1

          DO iBatch=1,nBatch

            If (iBatch.eq.nBatch) Then
               JNUM = nVrs - nVec*(nBatch-1)
            else
               JNUM = nVec
            endif
**********************************************************************
*                                                                    *
*           START WORKING                                            *
*                                                                    *
**********************************************************************
*
**          Read Cholesky vector
*
            JVEC = nVec*(iBatch-1) + iVrs
            IVEC2 = JVEC - 1 + JNUM

            CALL CHO_VECRD(Work(ipLrs),LREAD,JVEC,IVEC2,JSYM,
     &                     NUMV,IREDC,MUSED)

            If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
               rc=77
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

                ipLpq(k,1) = ipLF + lChoa   ! Lvb,J
                ipLpq(k,2) = ipLpq(k,1)     ! Lvi,J i general MO index
     &                     + nAsh(k)*nBas(i)*JNUM
                ipLpq(k,3) = ipLpq(k,2)     ! L~vi,J ~ transformed index
     &                       + nAsh(k)*nBas(i)*JNUM

                lChoa= lChoa + nAsh(k)*nBas(i)*3*JNUM

             End Do

             iSwap = 0  ! Lvb,J are returned
             kMOs = 1  !
             nMOs = 1  ! Active MOs (1st set)
*
**  Read half-transformed cho vectors
*
             ioff=0
             Do i=1,nSym
                k = Muld2h(i,JSYM)
                lvec=nAsh(k)*nBas(i)*JNUM
                iAdr=(JVEC-1)*nAsh(k)*nBas(i)+ioff
                call DDAFILE(LuAChoVec(Jsym),2,Work(ipLpq(k,1)),
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
     &                       One,Work(ipLvb),NAv,
     &                       Work(ipAOrb(iSymb,1)),Naw,
     &                      Zero,Work(ipLvw),NAv)
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
     &                  ONE,Work(ipLvtw),Nav*Naw,
     &                  Work(ipDA),1,ZERO,Work(ipVJ),1)
*
                 CALL DGEMV_('N',nRS,JNUM,
     &                -FactCI,Work(ipLrs),nRS,
     &                Work(ipVJ),1,1.0d0,Work(ipFab),1)

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
                  ipG    = ipG2
                  CALL DGEMV_('N',NAv*Naw,NAv*Naw,
     &               ONE,Work(ipG),NAv*Naw,
     &               Work(ipLvw),1,ZERO,Work(ipLxy),1)
*Qpx=Lpy Lxy
                  ipQpx=ipScr
                  Call DGEMM_('T','N',NBAS(iSymb),NAw,Nav,
     &                       One,Work(ipLvb),NAv,
     &                       Work(ipLxy),Naw,
     &                      ONE,Work(ipQpx),NBAS(iSymb))
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
     &                         ONE,Work(ipLvb),Nav,
     &                         Work(ipDA),Nav,ZERO,
     &                         Work(ipLwb),NBAS(iSymb))
                 End Do
*                 CALL CWTIME(TCINT2,TWINT2)
*                 tact(1) = tact(1) + (TCINT2 - TCINT3)
*                 tact(2) = tact(2) + (TWINT2 - TWINT3)
                 Do JVC=1,JNUM
                   ipLwb = ipLpq(iSymv,2)+ NAv*NBAS(iSymb)*(JVC-1)
                   Do is=1,NBAS(iSymb)
                    ipLtvb = ipLpq(iSymv,1)+ NAv*NBAS(iSymb)*(JVC-1)
     &                      + Nav*(is-1)
                    ipFock=ipKA+nBas(iSymb)*(is-1)
                    CALL DGEMV_('N',NBAS(iSymb),Nav,
     &                   -FactXI,Work(ipLwb),NBAS(iSymb),
     &                   Work(ipLtvb),1,ONE,Work(ipFock),1)

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
             mode = 'tofull'
             Call play_rassi_sto(irc,iLoc,JSYM,ISTLT,
     &                           ISSQ,ipJA,ipFab,mode)
          Call GetMem('rsFC','Free','Real',ipFab,nRS)
          EndIf
          Call GetMem('rsL','Free','Real',ipLrs,LREAD)
          CALL GetMem('FullV','Free','Real',ipLF,mTvec*nVec)
 999      Continue
        End Do  ! loop over red sets
 1000   CONTINUE
      End Do    ! loop over JSYM
**********************************************************************
*                                                                    *
*     POST PROCESSING                                                *
*                                                                    *
**********************************************************************
*
**    Accumulate Coulomb and Exchange contributions
*
      Do iSym=1,nSym
         ipFAc= ipJA + ISTLT(iSym)
         ipKAc= ipKA  + ISTSQ(iSym)
         ipFA = ipFkA + ISTSQ(iSym)

         Do iaSh=1,nShell
            ioffa = kOffSh(iaSh,iSym)

            Do ibSh=1,nShell

               ioffb = kOffSh(ibSh,iSym)

               Do ib=1,nBasSh(iSym,ibSh)
                Do ia=1,nBasSh(iSym,iaSh)
*MGD warning with sym

                  iag = ioffa + ia
                  ibg = ioffb + ib

                  jFA= ipFac- 1 + iTri(iag,ibg)
                  jKa = ipKac- 1 + nBas(iSym)*(ibg-1) + iag
                  jKa2= ipKac- 1 + nBas(iSym)*(iag-1) + ibg

                  jSA= ipFA - 1 + nBas(iSym)*(ibg-1) + iag
                  Work(jSA)= Work(jFa)+ Work(jKa)+ Work(jKa2)
                End Do
               End Do
            End Do
         End Do
      End Do

*
**Transform Fock and Q matrix to MO basis
*
      Do iS=1,nSym
        jS=iS
        If (nBas(iS).ne.0) Then
          Call DGEMM_('T','N',nBas(jS),nBas(iS),nBas(iS),
     &                1.0d0,Work(ipFkA+ISTSQ(iS)),nBas(iS),
     &                Work(ipCMO+ISTSQ(iS)),nBas(iS),0.0d0,
     &                Work(ipJA+ISTSQ(iS)),nBas(jS))
          call dcopy_(nBas(jS)*nBas(iS),[0.0d0],0,
     &                Work(ipFkA+ISTSQ(iS)),1)
          Call DGEMM_('T','N',nBas(jS),nIsh(jS),nBas(iS),
     &                1.0d0,Work(ipJA+ISTSQ(iS)),
     &                nBas(iS),Work(ipCMO+ISTSQ(jS)),nBas(jS),
     &                0.0d0,Work(ipFkA+ISTSQ(iS)),nBas(jS))
          ioff=nIsh(iS)*nBas(jS)+ISTSQ(iS)
          Call DGEMM_('T','N',nBas(jS),nAsh(iS),nBas(jS),
     &                 1.0d0,Work(ipCMO+ISTSQ(iS)),nBas(jS),
     &                 Work(ipScr+ISTSQ(iS)),nBas(jS),0.0d0,
     &                 Work(ipFkA+ioff),nBas(jS))
        EndIf
      End Do
**********************************************************************
*                                                                    *
*     TERMINATING                                                    *
*                                                                    *
**********************************************************************
      Call GetMem('ip_kOffSh','Free','Inte',ip_kOffSh,nShell*nSym)
      Call GetMem('Qmat','FREE','REAL',ipScr,nsBB*nDen)

      CAll QExit(SECNAM)
      Return
      END
