      Subroutine CASPT2_Res
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

C#include "SysDef.fh"
C
      !! 1) Calculate the derivative of the CASPT2 energy with respect
      !!    to the amplitude.
      !! 2) In the standard CASPT2, solve the CASPT2 equation. In the
      !!    diagonal CASPT2, compute the lambda directly.
C
C     write (*,*) "in CASPT2_res"
      !! Copy the solution vector to the residual space
      Do iCase = 1, 13
C       if (icase.ne.12.and.icase.ne.13) cycle
C       if (icase.ne.10.and.icase.ne.11) cycle
C       if (icase.ne. 8.and.icase.ne. 9) cycle
        Do iSym = 1, nSym
          nIN = nINDEP(iSym,iCase)
          IF(NIN.EQ.0) Cycle
          nAS = nASUP(iSym,iCase)
          nIS = nISUP(iSym,iCase)
C Remember: NIN values in BDIAG, but must read NAS for correct
C positioning.
          Call GETMEM('LBD','ALLO','REAL',LBD,nAS)
          Call GETMEM('LID','ALLO','REAL',LID,nIS)
          iD = iDBMat(iSym,iCase)
          Call dDaFile(LUSBT,2,Work(LBD),nAS,iD)
          Call dDaFile(LUSBT,2,Work(LID),nIS,iD)
C         if (icase.eq.4) then
C           do i = 1, nis
C             write (*,*) "ir = ",i
C             do j = 1, nin
C               write (*,'(i3,3f20.10)') i,work(lbd+j-1),work(lid+i-1),
C    *          1.0d+00/(work(lbd+j-1)+work(lid+i-1))
C             end do
C           end do
C         end if

          Call RHS_ALLO(nIN,nIS,lg_V1)
          Call RHS_ALLO(nIN,nIS,lg_V2)
          !! Read the solution vector
          Call RHS_Read(nIN,nIS,lg_V2,iCase,iSym,iVecX)
          !! Save it in the residual vector space immediately
C         Call RHS_Save(nIN,nIS,lg_V2,iCase,iSym,iVecR)
          !! Read the RHS vector
          Call RHS_Read(nIN,nIS,lg_V1,iCase,iSym,iRHS)
          If (MaxIt.eq.0) Then
            !! Scale the RHS vector appropriately (compute lambda)
            Call CASPT2_ResD(1,nIN,nIS,lg_V1,Work(LBD),Work(LID))
            !! T <- T + lambda
            Call DScal_(nIN*nIS,2.0D+00,Work(lg_V1),1)
C           call dcopy_(nin*nis,0.0d+00,0,work(lg_v2),1)
C           Call DaXpY_(nIN*nIS,2.0D+00,Work(lg_V1),1,Work(lg_V2),1)
            !! Save the modified T in the original T
C           Call RHS_Save(nIN,nIS,lg_V2,iCase,iSym,iVecX)
C           write (*,*) "lambda"
C       do i = 1, 10
C       write (*,*) i,work(lg_v1)
C       end do
            Call RHS_Save(nIN,nIS,lg_V1,iCase,iSym,iVecR)
          Else
          End If
          Call RHS_Free(nIN,nIS,lg_V1)
          Call RHS_Free(nIN,nIS,lg_V2)

          Call GETMEM('LBD','FREE','REAL',LBD,nAS)
          Call GETMEM('LID','FREE','REAL',LID,nIS)
        End Do
      End Do
C
      RETURN
      END
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE CASPT2_ResD(Mode,NIN,NIS,lg_W,DIN,DIS)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION DIN(*),DIS(*)

C Apply the resolvent of the diagonal part of H0 to an RHS array

#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
C-SVC: get the local vertical stripes of the lg_W vector
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)
        IF (iLo.NE.0.AND.jLo.NE.0) THEN
          NROW=iHi-iLo+1
          NCOL=jHi-jLo+1
          CALL GA_Access (lg_W,iLo,iHi,jLo,jHi,mW,LDW)
          CALL CASPT2_ResD2(MODE,NROW,NCOL,DBL_MB(mW),LDW,DIN(iLo),
     &                DIS(jLo),SHIFT,SHIFTI)
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
        CALL GA_Sync()
C       CALL GAdSUM_SCAL(DOVL)
      ELSE
        CALL CASPT2_ResD2(MODE,NIN,NIS,WORK(lg_W),NIN,DIN,DIS,
     &                   SHIFT,SHIFTI)
      END IF
#else
      CALL CASPT2_ResD2(MODE,NIN,NIS,WORK(lg_W),NIN,DIN,DIS,
     &                 SHIFT,SHIFTI)
#endif

      END
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
      SUBROUTINE CASPT2_ResD2(Mode,NROW,NCOL,W,LDW,DIN,DIS,
     &                  SHIFT,SHIFTI)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION W(LDW,*),DIN(*),DIS(*)

      DO J=1,NCOL
        DO I=1,NROW
          If (Mode.eq.1) Then
            DELTA  = SHIFT+DIN(I)+DIS(J)
            DELINV = DELTA/(DELTA**2+SHIFTI**2)
            !! The following SCAL is the actual residual
            SCAL   = 1.0D+00 - (DIN(I)+DIS(J))*DELINV
C           write (*,*) "residue = ", scal
C           if (abs(residue).ge.1.0d-08) write (*,*) "residue = ", scal
            !! Another scaling is required for lambda
            SCAL   =-SCAL*DELINV
          ELse If (Mode.eq.2) Then
            SCAL   =-SHIFTI/(DIN(I)+DIS(J))
          End If
          W(I,J) = SCAL*W(I,J)
        END DO
      END DO
      END
