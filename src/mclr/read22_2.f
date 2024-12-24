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
      SubRoutine Read22_2(MO1,Fock,Q,FockI,FockA,Temp2,Scr,Temp3)
********************************************************************
*                                                                  *
*   Constructs         everything                                  *
*                                                                  *
*                                                                  *
*                                                                  *
*   Output:MO     :MO integrals                                    *
*          Fock   :Fock matrix (one index transformed integrals)   *
*          MOtilde:MO (one index transformed integrals)            *
*                                                                  *
********************************************************************
      use Arrays, only: W_CMO=>CMO, W_CMO_Inv=>CMO_Inv, Int1, G1t, G2t
      use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
      Implicit Real*8(a-h,o-z)
#include "real.fh"
#include "Pointers.fh"
#include "standard_iounits.fh"
#include "Input.fh"
#include "stdalloc.fh"
#include "Files_mclr.fh"
      Real*8 Fock(nDens2),FockI(nDens2),FockA(nDens2),
     &       Temp2(nDens2),Temp3(ndens2),Q(nDens2),
     &       MO1(*), Scr(*)
      Logical Fake_CMO2,DoAct
      Real*8, Allocatable:: G2x(:)
      Type (DSBA_Type) CVa(2), DLT(1), DI, DA, Kappa, JI(1), KI, JA, KA,
     &                 FkI, FkA, QVec, CMO, CMO_Inv
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
      call dcopy_(ndens2,[0.0d0],0,focki,1)
      call dcopy_(ndens2,[0.0d0],0,focka,1)
      If(TwoStep.and.(StepType.eq.'RUN2')) Then
        nm=0
        Do iS=1,nSym
           nm=nAsh(is)+nm
        End Do
        nAtri=nm*(nm+1)/2
        nAtri=nAtri*(nAtri+1)/2
        Call dcopy_(ndens2,[0.0d0],0,fock,1)
        Call dcopy_(ndens2,[0.0d0],0,   Q,1)
        Call dcopy_(nAtri ,[0.0d0],0, MO1,1)
        Call ddafile(LuQDAT,2,FockA,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,2,FockI,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,2,Fock ,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,2,Q    ,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,2,MO1  ,nAtri ,iaddressQDAT)
        goto 101
      End If

*
      nas=0
      Do is=1,nSym
       nAS=nAS+nAsh(is)
      end do
*                                                                      *
************************************************************************
*                                                                      *
      Select Case (NewCho)
*                                                                      *
************************************************************************
*                                                                      *
       Case (.False.)   ! Cho-MO algorithm
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
         Do jS=1,iS
            ijS=iEOr(iS-1,jS-1)+1
            Do kS=1,nSym
               Do lS=1,ks
                  If (nOrb(iS)*nOrb(jS)*nOrb(ks)*nOrb(lS).eq.0)
     &                Go To 100
                  If (iEOr(kS-1,lS-1).ne.ijS-1) Goto 100
*                                                                      *
************************************************************************
*                                                                      *
                  ijB1=0
                  Do iB=1,nB(iS)
                     nnB=nB(jS)
                     If (iS.eq.jS) nnB=iB
                     Do jB=1,nnB
*
*                       Read a symmetry block
*
                        Call COUL(kS,lS,iS,jS,iB,jB,Temp2,Scr)
                        ijB1=ijB1+1
*                                                                      *
************************************************************************
*                                                                      *
*                       Add to unrotated inactive fock matrix
*
*                       Fkl sum(i)     2(ii|kl) -(ki|li)
*
*
*                       Coulomb term: F  =2(ii|kl)
*                                      kl
*
                        If (iS.eq.jS.and.iB.eq.jB .and.
     &                     (iB.le.nIsh(iS)))
     &                     Call DaXpY_(nOrb(kS)*nOrb(lS),Two,
     &                                Temp2,1,Focki(ipCM(kS)),1)

*                                                                      *
************************************************************************
*                                                                      *
*                       Add to unrotated active fock matrix
*
*                       Coulomb term: F  =2(ij|kl)d   i<j
*                                      kl          ij
*
                        If (iMethod.eq.iCASSCF) Then
*
                           If (iS.eq.jS) Then
                              If (((iB.gt.nIsh(is)).and.(nAsh(iS).ne.0))
     &                                .and.
     &                            ((jB.gt.nIsh(js)).and.(nAsh(jS).ne.0))
     &                           ) Then
                                 ipD=iTri(jB-nIsh(jS)+nA(jS),
     &                               iB-nIsh(is)+nA(iS))
                                 Fact=Two
                                 If (iB.eq.jB) Fact=one
                                 Call DaXpY_(nOrb(kS)*nOrb(lS),
     &                                      Fact*G1t(ipD),
     &                                      Temp2,1,FockA(ipCM(kS)),1)
                              End If
                           End If
*
                        End If
*                                                                      *
************************************************************************
*                                                                      *
                     End Do  ! jB
                  End Do     ! iB
*                                                                      *
************************************************************************
*                                                                      *
 100              Continue
               End Do          ! lS
            End Do             ! kS
         End Do                ! jS
      End Do                   ! iS
*                                                                      *
************************************************************************
*                                                                      *
*    Construct Q matrix: Q = sum(jkl)(pj|kl)d
*                         pi                 ijkl

      If (iMethod.eq.iCASSCF) Then

         Call CreQ2(Q,G2t,1,Temp2,Scr,nDens2)
*
*        Sort out MO (ij|kl)
*
         Do iS=1,nSym
            Do jS=1,iS
               Do kS=1,iS
                  lS=iEOR(iEOr(iS-1,jS-1),kS-1)+1
                  If (lS.gt.kS.or.(iS.eq.kS.and.lS.gt.jS)) Goto 123
*
                  Do iB=1,nAsh(iS)
                     iib=ib+nA(iS)
                     nnb=nAsh(jS)
                     If (iS.eq.jS) nnb=ib
                     Do jB=1,nnB
                        jjb=jb+nA(jS)
*
                        Call Coul(kS,lS,iS,jS,iB+nIsh(iS),jB+nIsh(jS),
     &                            Temp2,Scr)
*
                        nnK=nAsh(kS)
                        If (iS.eq.kS)  nnK=iB
                        Do kB=1,nnk
                           kkb=kb+nA(kS)
                           nnL=nAsh(lS)
                           If (kS.eq.lS)  nnL=kB
                           If (iib.eq.kkb) nnL=jB
                           Do lB=1,nnL
                              llb=lb+nA(lS)
*
                              ip2 = (lB+nIsh(lS)-1)*nBas(kS)
     &                            + kB+nIsh(kS)
*
                              ip1=iTri(iTri(iib,jjb),iTri(kkb,llb))
                              MO1(ip1) = Temp2(ip2)
                           End Do
                        End Do
                     End Do
                  End Do
*
 123              Continue
               End Do
            End Do
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Do iS=1,nSym
         kS=iS
         Do js=1,nSym
            lS=jS
            If (iEor(iEor(is-1,js-1),iEor(ks-1,ls-1)).ne.0) Cycle
            If (nOrb(iS)*nOrb(jS)*nOrb(ks)*nOrb(lS).eq.0)   Cycle
            Do LB=1,nB(LS)
               Do JB=1,nB(JS)
*                                                                      *
************************************************************************
*                                                                      *
                  Call EXCH(is,js,ks,ls,jb,lb,Temp2,Scr)
*                                                                      *
************************************************************************
*                                                                      *
*                 Add to unrotated inactive fock matrix
*
*                 Exchange term: F  =-(ij|kj)
*                                 ik
*
                  If (jS.eq.lS.and.jB.eq.lB.and.
     &                (jB.le.nIsh(jS)))
     &               Call DaXpY_(nOrb(iS)*nOrb(kS),-One,
     &                          Temp2,1,Focki(ipCM(iS)),1)

*                                                                      *
************************************************************************
*                                                                      *
*                 Add to unrotated active fock matrix
*
*                 Exchange term: F  =-1/2(ij|kl)d
*                                 ik             jl
*
                  If (iMethod.eq.iCASSCF) Then
                     If (jS.eq.lS) Then
                        If (((jB.gt.nIsh(js)).and.(nAsh(jS).ne.0)).and.
     &                      ((lB.gt.nIsh(ls)).and.(nAsh(lS).ne.0))
     &                     ) Then
                          ipD=iTri(lB-nIsh(lS)+nA(lS),
     &                             jB-nIsh(js)+nA(jS))
                          Call DaXpY_(nOrb(iS)*nOrb(kS),-half*G1t(ipD),
     &                               Temp2,1,FockA(ipCM(iS)),1)
                        End If
                     End If
                  End If
*                                                                      *
************************************************************************
*                                                                      *
               End Do
            End Do
         End Do ! jS
      End Do ! iS
*                                                                      *
************************************************************************
*                                                                      *
      Case (.TRUE.)   ! Cho-Fock Algorithm
*                                                                      *
************************************************************************
*                                                                      *
        Fake_CMO2=.true.
        DoAct=.true.
*
**      Construct inactive density matrix
*
        call dcopy_(nDens2,[0.0d0],0,temp2,1)
        Do is=1,nSym
          Do iB=1,nIsh(is)
            ip=ipCM(iS)+(ib-1)*nOrb(is)+ib-1
            Temp2(ip)=2.0d0
          End Do
        End Do
*
**      Transform to AO basis
*
        Do iS=1,nSym
           If (nIsh(iS).ne.0) Then
              jS=iS
              Call DGEMM_('T','T',nIsh(jS),nOrb(iS),nIsh(iS),
     &                    1.0d0,Temp2(ipCM(iS)),nOrb(iS),
     &                    W_CMO(ipCM(is)),nOrb(iS),
     &                    0.0d0,Temp3(ipMat(jS,iS)),nOrb(jS))
              Call DGEMM_('T','T',nOrb(jS),nOrb(jS),nIsh(iS),
     &                    1.0d0,Temp3(ipMat(jS,iS)),nOrb(iS),
     &                    W_CMO(ipCM(js)),nOrb(jS),
     &                    0.0d0,Temp2(ipCM(iS)),nOrb(jS))
           EndIf
        End Do
*
        Call Allocate_DT(DLT(1),nBas,nBas,nSym,aCase='TRI')
        call Fold_Mat(nSym,nOrb,Temp2,DLT(1)%A0)
*
**      Form active CMO and density
*
        nAct=0
        If (iMethod.eq.iCASSCF) Then
          na2=0
          nG2=0
          Do iSym=1,nSym
            na2=na2+nAsh(iSym)**2
            nAct=nAct+nAsh(iSym)
            nAG2=0
            Do jSym=1,nSym
              kSym=iEOr(jsym-1,isym-1)+1
              nAG2=nAg2+nAsh(jSym)*nAsh(kSym)
            End Do
            nG2=nG2+nAG2**2
          End Do
          Call Allocate_DT(CVa(1),nAsh,nOrb,nSym)
          CVa(1)%A0(:)=0.0D0
          Call Allocate_DT(CVa(2),nAsh,nOrb,nSym)
          CVa(2)%A0(:)=0.0D0
          Call Allocate_DT(DA,nAsh,nAsh,nSym)
*
          ioff=0
          Do iSym=1,nSym
            ioff2 = ioff + nOrb(iSym)*nIsh(iSym)
            do ikk=1,nAsh(iSym)
               ioff3=ioff2+nOrb(iSym)*(ikk-1)
               CVa(1)%SB(iSym)%A2(ikk,:) =
     &            W_CMO(ioff3+1:ioff3+nOrb(iSym))
               ik=ikk+nA(iSym)
               Do ill=1,ikk-1
                 il=ill+nA(iSym)
                 ikl=ik*(ik-1)/2+il
                 DA%SB(iSym)%A2(ill,ikk)=G1t(ikl)
                 DA%SB(iSym)%A2(ikk,ill)=G1t(ikl)
               End Do
               ikl=ik*(ik-1)/2+ik
               DA%SB(iSym)%A2(ikk,ikk)=G1t(ikl)
            End Do
            ioff=ioff+nOrb(iSym)**2
          End Do
          DA%A0(:) = Half * DA%A0(:)
*
**      Expand 2-body density matrix
*
          Call mma_allocate(G2x,nG2,Label='G2x')
          ipGx=0
          Do ijS=1,nSym
            Do iS=1,nSym
              jS=iEOR(is-1,ijS-1)+1
              Do kS=1,nSym
                lS=iEOR(kS-1,ijS-1)+1
                Do kAsh=1,nAsh(ks)
                  Do lAsh=1,nAsh(ls)
                    ikl=itri(lAsh+nA(lS),kAsh+nA(kS))
                    Do iAsh=1,nAsh(is)
                      Do jAsh=1,nAsh(js)
                        iij =itri(iAsh+nA(is),jAsh+nA(jS))
                        ipGx=ipGx+1
                        G2x(ipGx)=G2t(itri(iij,ikl))
                      End Do
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
        Else
          na2=1
          nG2=1
          Call Allocate_DT(CVa(1),[1],[1],1) ! dummy allocation
          Call Allocate_DT(CVa(2),[1],[1],1)
          Call Allocate_DT(DA,nAsh,nAsh,nSym)
          Call mma_allocate(G2x,nG2,Label='G2x')
        EndIf
*
**      Let's go
*
        Call Allocate_DT(JA,nBas,nBas,nSym)
        JA%A0(:)=Zero
        Call Allocate_DT(KA,nBas,nBas,nSym)
        KA%A0(:)=Zero
*
        call dcopy_(nDens2,[0.0d0],0,Q,1)
*
        Call Allocate_DT(DI,nBas,nBas,nSym,Ref=Temp2)
        Call Allocate_DT(JI(1),nBas,nBas,nSym,aCase='TRI',Ref=Temp3)
        JI(1)%A0(:)=Zero
        Call Allocate_DT(KI,nBas,nBas,nSym,Ref=Scr)
        KI%A0(:)=Zero
        Call Allocate_DT(FkI,nBas,nBas,nSym,Ref=FockI)
        FkI%A0(:)=Zero
        Call Allocate_DT(FkA,nBas,nBas,nSym,Ref=FockA)
        FkA%A0(:)=Zero
        Call Allocate_DT(QVec,nBas,nAsh,nSym,Ref=Q)
        Call Allocate_DT(CMO,nBas,nBas,nSym,Ref=W_CMO)
        Call Allocate_DT(CMO_Inv,nBas,nBas,nSym,Ref=W_CMO_Inv)
        istore=1 ! Ask to store the half-transformed vectors

        CALL CHO_LK_MCLR(DLT,DI,DA,G2x,Kappa,JI,KI,JA,KA,FkI,FkA,
     &                   MO1,QVec,CVa,CMO,CMO_inv,
     &                   nIsh,nAsh,doAct,Fake_CMO2,
     &                   LuAChoVec,LuIChoVec,istore)

        nAtri=nAct*(nAct+1)/2
        nAtri=nAtri*(nAtri+1)/2
        Call DScal_(nAtri,0.25D0,MO1,1)
        FkI%A0(:) = -Half * FkI%A0(:)
*
        Call Deallocate_DT(CMO_Inv)
        Call Deallocate_DT(CMO)
        Call Deallocate_DT(QVec)
        Call Deallocate_DT(FkA)
        Call Deallocate_DT(FkI)
        Call Deallocate_DT(KI)
        Call Deallocate_DT(JI(1))
        Call Deallocate_DT(DI)
        Call Deallocate_DT(JA)
        Call Deallocate_DT(KA)
        Call deallocate_DT(DLT(1))
        Call mma_deallocate(G2x)
        Call Deallocate_DT(CVa(2))
        Call Deallocate_DT(CVa(1))
        Call deallocate_DT(DA)
*
        Call GADSum(FockI,nDens2)
        Call GADSum(FockA,nDens2)
        Call GADSum(    Q,nDens2)
        Call GADSum(  MO1,nAtri)
*                                                                      *
************************************************************************
*                                                                      *
      End Select
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Do iSym = 1, nSym
        Write (6,*) 'iSym=',iSym
        Call RecPrt('FockI',' ',FockI(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
        Call RecPrt('FockA',' ',FockA(ipCM(iSym)),nOrb(iSym),nIsh(iSym))
        Call RecPrt('Q',' ',Q(ipMatba(iSym,iSym)),nOrb(iSym),nAsh(iSym))
      End Do
      nAtri=nas*(nas+1)/2
      nAtri=nAtri*(nAtri+1)/2
      Call RecPrt('MO1',' ',MO1,1,nAtri)
#endif
      Call DaXpY_(ndens2,One,Int1,1,FockI,1)
      call dcopy_(ndens2,[0.0d0],0,Fock,1)
*
      Do iS=1,nSym
         If (nOrb(iS).eq.0) Go To 300
*
         If (nIsh(iS).gt.0)
     &      Call DYaX(nOrb(iS)*nIsh(is),2.0d0,
     &                FockI(ipCM(iS)),1,
     &                Fock (ipCM(iS)),1)
         If (iMethod.eq.iCASSCF) Then
            If (nIsh(iS).gt.0)
     &         Call DaXpY_(nOrb(iS)*nIsh(is),2.0d0,
     &                    FockA(ipCM(iS)),1,
     &                    Fock (ipCM(iS)),1)
            If (nAsh(iS).gt.0)
     &         Call DYaX(nOrb(iS)*nAsh(is),1.0d0,
     &                   Q(ipMatba(iS,is)),1,
     &                   Fock(ipCM(iS)+nIsh(is)*nOrb(is)),1)
            Do iAsh=1,nAsh(is)
               ipi=ipCM(iS)+nOrb(is)*(nIsh(is)+iAsh-1)
               Do jAsh=1,nAsh(is)
                  ipj=ipCM(iS)+nOrb(is)*(nIsh(is)+jAsh-1)
                  ni=nA(is)+iAsh
                  nj=nA(is)+jAsh
                  ipD=iTri(ni,nj)
                  call daxpy_(nOrb(is),G1t(ipD),
     &                       FockI(ipi),1,
     &                       Fock (ipj),1)
               End Do
            End Do
         End If
*
 300     Continue
      End Do
 101  Continue
      renergy=0.0d0
      rcora=0.0d0
      Do iS=1,nSym
      Do iB=1,nAsh(is)+nIsh(is)
      rEnergy=rEnergy+Fock(ipCM(is)+nOrb(iS)*(ib-1)+ib-1)
      End Do
      End Do
      rcorei=0.0d0
      rcorea=0.0d0
      rcor=0.0d0
      Do iS=1,nSym
       iptmp=ipCM(iS)
       Do iB=1,nIsh(is)
       rcorei=rcorei+2.0d0*Int1(iptmp)
       rcor=rcor+2.0d0*Focki(iptmp)
       iptmp=iptmp+nOrb(iS)+1
       End Do

       Do iB=1,nAsh(iS)
        Do jB=1,nAsh(iS)
         iiB=nA(iS)+ib
         ijB=nA(iS)+jb
         iij=iTri(iib,ijb)
         iiB=nIsh(iS)+ib
         ijB=nIsh(iS)+jb
         rcorea=rcorea+G1t(iij)*Int1(ipCM(is)-1+nOrb(is)*(iib-1)+ijB)

         rcora=rcora+G1t(iij)*Focki(ipCM(is)+nOrb(is)*(iib-1)+ijB-1)
        End Do
       End Do
      End Do
      rin_ene=0.5d0*(rcor+rcorei)
      rcore=rCorei+rcoreA
      If (debug) Then
         Write(6,*) 'Checking energy',0.5d0*renergy+potnuc+half*rcore
         Write(6,*) 'Checking energy',0.5d0*renergy,potnuc,half*rcore
         write(6,*)
      End if

      If(TwoStep.and.(StepType.eq.'RUN1')) Then
        iaddressQDAT=0
        nm=0
        Do iS=1,nSym
           nm=nAsh(is)+nm
        End Do
        nAtri=nm*(nm+1)/2
        nAtri=nAtri*(nAtri+1)/2
        Call ddafile(LuQDAT,1,FockA,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,1,FockI,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,1,Fock ,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,1,Q    ,nDens2,iaddressQDAT)
        Call ddafile(LuQDAT,1,MO1  ,nAtri ,iaddressQDAT)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
