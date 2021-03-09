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
      SUBROUTINE CHO_CAS_DRV(rc,CMO,DI,FI,DA1,FA,DA2,TraOnly)
      Use Data_Structures, only: CMO_Type, Allocate_CMO,
     &                           Deallocate_CMO
      Implicit real*8 (a-h,o-z)

      Integer   rc
      Real*8    DA1(*),DI(*),DA2(*),FI(*),FA(*),CMO(*)
      Integer   nForb(8),nIorb(8),nAorb(8),nChM(8),nChI(8)
      Integer   ipDSA2(8,8,8),nnA(8,8),ipKLT(2)
      Logical   TraOnly

#include "chotodo.fh"
#include "chopmat.fh"
#include "chlcas.fh"
#include "cholk.fh"

      Character(LEN=11), Parameter:: SECNAM = 'CHO_CAS_DRV'

#include "rasdim.fh"
#include "wadr.fh"
#include "general.fh"
#include "rasscf.fh"
#include "WrkSpc.fh"

      Type (CMO_type) CVa(2), POrb(3)
C ************************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
C  **************************************************

      rc=0

      IF (TraOnly) THEN
c
c --- It only performs the MO transformation of FI and FA
c -------------------------------------------------------
c
*     transform FI from AO to MO basis  (LT-storage)
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        Call GetMem('Scr1','Allo','Real',iTmp1,iBas*iBas)
        Call GetMem('Scr2','Allo','Real',iTmp2,iOrb*iBas)
        Call Square(FI(iOff1),Work(iTmp1),1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Work(iTmp1),
     &               iBas,CMO(iOff2+(iFro*iBas)),max(iBas,iBas),
     &               0.0d0,Work(iTmp2),iBas)
        Call MXMT(Work(iTmp2),iBas,1,
     &            CMO(iOff2+(iFro*iBas)),1,iBas,
     &            FI(iOff3),
     &            iOrb,iBas)
        Call GetMem('Scr2','Free','Real',iTmp2,iOrb*iBas)
        Call GetMem('Scr1','Free','Real',iTmp1,iBas*iBas)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

*     transform FA from AO to MO basis  (LT-storage)
      iOff1 = 1
      iOff2 = 1
      iOff3 = 1
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        Call GetMem('Scr1','Allo','Real',iTmp1,iBas*iBas)
        Call GetMem('Scr2','Allo','Real',iTmp2,iOrb*iBas)
        Call Square(FA(iOff1),Work(iTmp1),1,iBas,iBas)
        Call DGEMM_('N','N',iBas,iOrb,iBas,1.0d0,Work(iTmp1),
     &               iBas,CMO(iOff2+(iFro*iBas)),max(iBas,iBas),
     &               0.0d0,Work(iTmp2),iBas)

        Call MXMT(Work(iTmp2),iBas,1,
     &            CMO(iOff2+(iFro*iBas)),1,iBas,
     &            FA(iOff3),
     &            iOrb,iBas)
        Call GetMem('Scr2','Free','Real',iTmp2,iOrb*iBas)
        Call GetMem('Scr1','Free','Real',iTmp1,iBas*iBas)
        iOff1 = iOff1 + (iBas*iBas+iBas)/2
        iOff2 = iOff2 + iBas*iBas
        iOff3 = iOff3 + (iOrb*iOrb+iOrb)/2
      End Do

c**************************************************************************


      ELSE

c --- It only computes FI and FA  in AO-basis and returns
c --- the active integrals (tw|xy)
c --- If specified in input, the routine also computes
c --- the auxiliary Q-matrix stored as Q(av), where a is an AO index
c --- and v refers to the active orbitals only

      Do iSym=1,nSym
         nForb(iSym) = nFro(iSym)
         nIorb(iSym) = nIsh(iSym)
         nAorb(iSym) = nAsh(iSym)
      End Do


C --- Build the packed densities from the Squared ones
      Call Getmem('DILT','Allo','Real',ipDILT,NTot1)
      Call Getmem('DALT','Allo','Real',ipDALT,NTot1)

      Call Fold(nSym,nBas,DI,Work(ipDILT))
      Call Fold(nSym,nBas,DA1,Work(ipDALT))

      FactXI = -1.0D0

!AMS - should this be set differently for ExFac.ne.1?
!      FactXI = 0-ExFac

      If (Deco) Then

         FactXI = -0.5D0

!AMS - should this be set differently for ExFac.ne.1?
!         FactXI = 0-(ExFac*.5d0)

c --- decompose the Inactive density on request
         CALL GETMEM('choIn','allo','real',ipIna,NTot2)
         CALL GETMEM('ddec','allo','real',ipddec,NTot2)
         call dcopy_(NTot2,DI(1),1,Work(ipddec),1)

         ipInc = ipIna

         Thr = 1.0d-12
         ipd = ipddec
         ipV = ipIna
         incs=0
         Do i=1,nSym
            if((nForb(i)+nIorb(i)).gt.0)then
             CALL CD_InCore(Work(ipd),nBas(i),Work(ipV),nBas(i),
     &                     NumV,Thr,rc)
             If (rc.ne.0) Then
              write(6,*)SECNAM//': ill-defined dens decomp for Inact'
              write(6,*) 'rc value produced = ', rc
              Call abend()
             EndIf
             nChI(i) = NumV
             if ( NumV .ne. nIsh(i)+nForb(i) ) then
               write(6,*)'Warning! The number of occupied from the deco'
     &                 //'mposition of the Inactive density matrix is ',
     &                   numV,' in symm. ',i
               write(6,*)'Expected value = ',nIsh(i)+nForb(i)
               incs=incs+1
               Ymax=0.0d0
               do ja=1,nBas(i)
                  jaa=ipd-1+nBas(i)*(ja-1)+ja
                  Ymax=Max(Ymax,Work(jaa))
               end do
               write(6,*)'Max diagonal of the density in symm. ',i,' is'
     &                   //' equal to ',Ymax
             endif
            else
             nChI(i) = 0
            endif
            ipd = ipd + nBas(i)**2
            ipV = ipV + nBas(i)**2
         End Do

         If (incs.gt.0 .and. DoLocK) then
            dmpk_old = dmpk
            dmpk = 1.0d-2*dmpk
            write(6,*)'LK-damping decreased from ',dmpk_old,' to ',dmpk
         EndIf

         CALL GETMEM('ddec','free','real',ipddec,NTot2)

c --- to get the right input arguments for CHO_FCAS_AO and CHO_FMCSCF
         If(.not.DoLocK)Then
           Do i=1,nSym
              nForb(i) = 0
              nIorb(i) = nChI(i)
           End Do
         EndIf


      Else

        ipInc = ip_of_Work(CMO(1))

        Do i=1,nSym
           nChI(i) = nForb(i)+nIorb(i)
        End Do

      EndIf

C --- Reordering of the MOs coefficients to fit cholesky needs
      If(.not.DoLocK)Then

        Call Allocate_CMO(POrb(1),nChI,nBas,nSym)
        Call Allocate_CMO(POrb(3),nAOrb,nBas,nSym)

        nOcs=0
        ioff1=0
        Do iSym=1,nSym

           do ikk=1,nChI(iSym)
              ioff2=ioff1+nBas(iSym)*(ikk-1)
              POrb(1)%SB(iSym)%A(ikk,:) =
     &           Work(ipInc+ioff2 : ipInc+ioff2 -1 + nBas(iSym))
           end do

           ioff2=ioff1+nBas(iSym)*(nForb(iSym)+nIorb(iSym))
           do ikk=1,nAorb(iSym)
              POrb(3)%SB(iSym)%A(ikk,:) =
     &             CMO( ioff2+nBas(iSym)*(ikk-1) + 1 :
     &                  ioff2+nBas(iSym)*(ikk-1) + nBas(iSym))
           end do
           ioff1=ioff1+nBas(iSym)**2
           nOcs = nOcs + nAorb(iSym)**2

        End Do

      Else

C *** Only the active orbitals MO coeff need reordering
           Call Allocate_CMO(CVa(1),nAorb,nBas,nSym)

           ioff1 = 0
           Do iSym=1,nSym
            ioff2 = ioff1 + nBas(iSym)*(nForb(iSym)+nIorb(iSym))
            do ikk=1,nAorb(iSym)
               ioff = ioff2+nBas(iSym)*(ikk-1)
               CVa(1)%SB(iSym)%A(ikk,:) =
     &           CMO(ioff +  1 : ioff + nBas(iSym))
            end do
            ioff1 = ioff1 + nBas(iSym)**2
           End Do

      EndIf

C --- Optional Section for Q-matrix evaluation
C --- Reorder the 2-el density matrix to fit cholesky needs
      If (DoQmat.and.ALGO.ne.1) Then

         Call set_nnA(nSym,nAorb,nnA)

         nPmat=0   ! P[vw],xy
         Do iSymXY=1,nSym
            Do iSymy=1,nSym
               iSymx=MulD2h(iSymXY,iSymy)
               if (iSymx.le.iSymy) then
               Do iSymw=1,nSym
                  iSymv=MulD2h(iSymXY,iSymw)
                  nPmat = nPmat
     &                  + nAorb(iSymv)*nAorb(iSymw)*nnA(iSymx,iSymy)
               End Do
               endif
            End do
         End Do

         Call Getmem('P-mat','Allo','Real',ipPmat,nPmat)
         Call Getmem('P-scal','Allo','Real',ipPL,NACPR2)
         Call CHO_Pmat(DA2,Work(ipPmat))

         Call Fzero(Work(ipPmat),nPmat)

c         ipDA2 = ip_of_Work(DA2(1))
         ipDA2 = ipPL

         Call Reord_Pmat(ipDA2,ipPmat,ipDSA2)

         Call Getmem('P-scal','Free','Real',ipPL,NACPR2)

      EndIf



      If (DoActive) Then
C ---  Decompose the active density  -----------------------------

#ifdef _DEBUGPRINT_
       koff=0
       do i=1,nSym
          CALL CD_TESTER(rc,Work(ipDALT+koff),nBas(i),.true.)
          write(6,*) 'DALT for sym=', i
          CALL TRIPRT('DALT',' ',Work(ipDALT+koff),nBas(i))
          koff = koff + nBas(i)*(nBas(i)+1)/2
       end do
#endif

        Call Allocate_CMO(CVa(2),nBas,nBas,nSym)
        CALL GETMEM('ddec','allo','real',ipddec,NTot2)
        call dcopy_(NTot2,DA1(1),1,Work(ipddec),1)

        Thr = 1.0d-12
        ipd = ipddec
        Do i=1,nSym
           if(nAorb(i).gt.0)then
! NOTE(Giovanni): CD will proceed with approx. decompos for QMC
!                 This will avoid warnings for negative-definit
             call CD_InCore(Work(ipd),nBas(i),
     &                      CVa(2)%SB(i)%A,nBas(i),
     &                      NumV,Thr,rc)
             If (rc.ne.0) Then
                write(6,*)SECNAM//': ill-defined dens decomp for active'
                write(6,*) 'rc value produced = ', rc
                Call abend()
             EndIf
             nChM(i) = NumV
           else
             nChM(i) = 0
           endif
           ipd = ipd + nBas(i)**2
        End Do

        CALL GETMEM('ddec','free','real',ipddec,NTot2)

      Else

        ! Dummy allocation
        Call Allocate_CMO(CVa(2),[1],[1],1)
        nChM(:) = 0

      EndIf

      If (.not.DoLocK .and. DoActive) Then

c --- reorder "Cholesky MOs" to Cva storage

        Call Allocate_CMO(POrb(2),nChM,nBas,nSym)
        Do iSym=1,nSym
           If (nBas(iSym)*nChM(iSym).ne.0) Then
               do ikk=1,nChM(iSym)
                  POrb(2)%SB(iSym)%A(ikk,:) =
     &               CVa(2)%SB(iSym)%A(:,ikk)
               end do
           EndIf
        End Do

      Else

        Call Allocate_CMO(POrb(2),[1],[1],1)

      EndIf
C ----------------------------------------------------------------

      Call Fzero(FI(1),nTot1) ! LT-storage
      Call Fzero(FA(1),nTot1) ! LT-storage

      ipFI = ip_of_Work(FI(1))
      ipFA = ip_of_Work(FA(1))


      IF (ALGO.eq.1 .and. .not. DoLocK) THEN

         ipInt = lpwxy   ! (PU|VX) integrals are computed
         ipCM = ip_of_work(CMO(1))  ! MOs coeff. in C(a,p) storage

         CALL CHO_FMCSCF(rc,ipFA,ipFI,nForb,nIorb,nAorb,FactXI,
     &                   ipDILT,ipDALT,DoActive,POrb,nChM,ipInt,ExFac)


      ELSEIF (ALGO.eq.1 .and. DoLocK) THEN

         ipInt = lpwxy   ! (PU|VX) integrals are computed
         ipCM = ip_of_work(CMO(1))  ! MOs coeff. in C(a,p) storage

         Call Getmem('KILT','Allo','Real',ipKLT(1),NTot1)
         Call Fzero(Work(ipKLT(1)),NTot1)

         If (DoActive) Then
           Call Getmem('KALT','Allo','Real',ipKLT(2),NTot1)
           Call Fzero(Work(ipKLT(2)),NTot1)
         EndIf

         CALL CHO_LK_CASSCF(ipDILT,ipDALT,ipFI,ipFA,ipKLT,ipInc,ipInt,
     &                      FactXI,nChI,nAorb,nChM,CVa,DoActive,
     &                      nScreen,dmpK,abs(CBLBM),ExFac)

         If(DoActive) Call Getmem('KALT','Free','Real',ipKLT(2),NTot1)
         Call Getmem('KILT','Free','Real',ipKLT(1),NTot1)

      ELSE

         write(6,*)SECNAM//': wrong input parameter. ALGO= ',ALGO
         rc=55
         Return

      ENDIF

      If (Allocated(POrb(3)%A0)) Call Deallocate_CMO(POrb(3))
      If (Allocated(POrb(2)%A0)) Call Deallocate_CMO(POrb(2))
      If (Allocated(POrb(1)%A0)) Call Deallocate_CMO(POrb(1))
      If (Allocated(CVa(1)%A0)) Call Deallocate_CMO(CVa(1))
      If (Allocated(CVa(2)%A0)) Call Deallocate_CMO(CVa(2))

      If (DoQmat.and.ALGO.ne.1) Then
         Call Getmem('P-mat','Free','Real',ipPmat,NPmat)
      EndIf

      If (Deco) CALL GETMEM('choIn','free','real',ipIna,NTot2)

      Call Getmem('DALT','Free','Real',ipDALT,NTot1)
      Call Getmem('DILT','Free','Real',ipDILT,NTot1)



      ENDIF


      Return
      END
