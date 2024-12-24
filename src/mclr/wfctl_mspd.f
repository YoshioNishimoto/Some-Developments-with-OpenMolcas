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
* Copyright (C) 2021, Jie J. Bao                                       *
************************************************************************
      SubRoutine WfCtl_MSPD(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,
     &                    iRHSDisp,converged,iPL)
************************************************************************
*                                                                      *
*                                                                      *
*     called from: MCLR                                                *
*                                                                      *
*                                                                      *
************************************************************************
      use Exp, only: Exp_Close
      use ipPage, only: W
      use cmslag, only: ResQaaLag2
      Implicit Real*8 (a-h,o-z)
*
#include "stdalloc.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"

      Logical CI
#include "crun_mclr.fh"
      Character*8   Fmt2
      Integer iKapDisp(nDisp),isigDisp(nDisp)
      Integer iRHSDisp(nDisp)
      Integer iCIDisp(nDisp),iCIsigDisp(nDisp)
      Integer opOut
      Logical lPrint,converged(8)
      Real*8 rchc(mxroot)
      Real*8 rdum(1)
      Real*8, Allocatable:: Kappa(:), dKappa(:), Sigma(:),
     &                      Temp3(:), Temp4(:),
     &                      Sc1(:), Sc2(:), Fancy(:)
      Integer LURot,IsFreeUnit,JRoot
      External IsFreeUnit
      CHARACTER(len=16)::VecName

      interface
        subroutine RHS_NAC(Fock,SLag)
          Real*8 Fock(*)
          real*8, optional :: SLag(*)
        end subroutine
      end interface

*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      Call StatusLine(' MCLR:',
     &                ' Computing Lagrangian multipliers for MS-PDFT')
*

      lPaper=132
      lLine =120
      left=(lPaper-lLine)/2
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*

      iDis=0
*
      fail=.false.
      Do i=1,8
       Converged(i)=.true.
      end do
*MGD I think this is nice when printed...
      lprint=.true.
      reco=-One
      Lu_50=50
      If (SAVE) CALL DANAME(Lu_50,'RESIDUALS')
      If (SAVE) Then
         Write (6,*) 'WfCtl_MSPD: SAVE option not implemented'
         Call Abend()
      End If
      If (iAnd(kprint,2).eq.2) lprint=.true.
      isym=1
      nconf1=ncsf(State_Sym)

      CI=.false.
      If (iMethod.eq.2.and.nconf1.gt.0) CI=.true.

*          Initiate CSF <-> SD
           Call InCSFSD(iEor(iSym-1,State_Sym-1)+1,
     &                  State_sym,.false.)
*
*
*          Calculate length of the density, Fock and Kappa matrix etc
*          notice that this matrices are not necessarily symmetric.
*          Store pointers.
*
*          Input:
*                 iSym: Symmetry of perturbation
*
*          Output: Commonblocks (Pointers.fh)
*
      nConf3=nint(Max(xispsm(State_SYM,1),xispsm(State_SYM,1)))

      Call Setup_MCLR(iSym)

*
*     Determine if we should page CI vectors
*                                [2]
*     Calculate the diagonal of E    and store in core/disc
*
      Call mma_allocate(FANCY,nroots**3,Label='FANCY')
      Call CIDia_SA(State_Sym,rCHC,Fancy)

      irc=ipOut(ipdia)
*
*     Allocate disk/memory space
*
*
*     This areas should be addressed through ipIn
*     ipOut will page them out to disk and free the memory area
*     if page mode is used
*
*     opOut will release the memory area without update the disk
*
      ipS1 =ipGet(nconf3*nroots)
      ipS2 =ipGet(nconf3*nroots)
      ipST =ipGet(nconf3*nroots)
      ipCIT=ipGet(nconf1*nroots)
      ipCId=ipGet(nconf1*nroots)
*
      npre2=npre(isym)
      ipPre2=ipGet(npre2)
      irc=ipIn(ipPre2)
      If (TwoStep.and.(StepType.eq.'RUN2')) Then
         ! fetch data from LuQDAT and skip the call to "Prec"
         Call ddafile(LuQDAT,2,W(ipPre2)%Vec,npre2,iaddressQDAT)
      Else
         Call Prec(W(ipPre2)%Vec,isym)
         irc=ipOut(ippre2)
         If (TwoStep.and.(StepType.eq.'RUN1')) Then
            ! save the computed data in "Prec" to LuQDAT and skip the
            ! following part of this function
            irc=ipIn(ipPre2)
            Call ddafile(LuQDAT,1,W(ipPre2)%Vec,npre2,iaddressQDAT)
            Go To 193
         End If
      End If
*
*
*     OK START WORKING
*
*     idisp=1
      jspin=0
*
*     Allocate areas for scratch and state variables
*
      Call mma_allocate(Kappa,nDens2+6,Label='Kappa')
      Call mma_allocate(dKappa,nDens2+6,Label='dKappa')
      Call mma_allocate(Sigma,nDens2+6,Label='Sigma')
      Call mma_allocate(Temp3,nDens2+6,Label='Temp3')
      Call mma_allocate(Temp4,nDens2+6,Label='Temp4')
      Call mma_allocate(Sc1,nDens2+6,Label='Sc1')
      Call mma_allocate(Sc2,nDens2+6,Label='Sc2')
*
      do iDisp=1,nDisp
         Kappa(1:nDens2)=Zero
         dKappa(1:nDens2)=Zero
         Sigma(1:nDens2)=Zero
*
*-----------------------------------------------------------------------------
*
*     Calculate RHS for the perturbation
*
*-----------------------------------------------------------------------------
*
*     (T1,T2,T3,T4,T5,T6,T7,Kappa1,CI1)
*
      If (debug) Then
         If (isNAC) Then
            Write(6,*)'States: ',NACstates(1),NACstates(2)
         Else
            Write(6,*)'State: ',irlxroot
         EndIf
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      If (isNAC) Then
        Call RHS_NAC(Temp4)
      Else
       LURot=233
       LURot=IsFreeUnit(LURot)
       Call Molcas_Open(LURot,'ROT_VEC')
       DO JRoot=1,nRoots
        read(LURot,*) VecName
       END DO
       read(LURot,*) VecName
       CLOSE(LURot)
       if(VecName.eq.'CMS-PDFT') THEN
        Call RHS_CMS(Temp4,W(ipST)%Vec)
        CALL DMinvCI_SA(ipST,W(ipS2)%Vec,rdum(1),isym,Fancy)
       else
        write(6,'(6X,A)')'Error: Lagrangian Not Implemented for MS-PDFT'
        write(6,'(6X,A)')'       Other Than CMS-PDFT'
        CALL xQuit(96)
       end if
      End If


*
      irc=opOut(ipci)
*
      If (lprint) Write(6,*)
     &      '       Iteration       Delta           Res(kappa)     '//
     &      '  Res(CI)          DeltaK           DeltaC'
      iLen=nDensC
      iRHSDisp(iDisp)=iDis
      Call Compress(Temp4,Sigma,iSym)
      r1=ddot_(nDensc,Sigma,1,Sigma,1)
      If(debug)Write(6,*) 'Hi how about r1',r1
      Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
*
      irc=ipIn(ipCIT)
*      irc=ipIn(ipST)
      irc=ipIn(ipCID)
      call dcopy_(nConf1*nroots,[Zero],0,W(ipCIT)%Vec,1)
*     already initialized in rhs_mspdft
*      call dcopy_(nConf1*nroots,[Zero],0,W(ipST)%Vec,1)
*      call dcopy_(nConf1*nroots,[Zero],0,W(ipCID)%Vec,1)
*      call dcopy_(nConf1*nroots,W(ipST)%Vec,1,W(ipCID)%Vec,1)
      irc=ipOut(ipCIT)
      Call DSCAL_(nDensC,-One,Sigma,1)
*
      irc=ipIn(ipPre2)
      Call DMInvKap(W(ipPre2)%Vec,Sigma,nDens2+6,
     &              Kappa,nDens2+6,Temp3,nDens2+6,
     &              isym,iter)

      irc=opOut(ippre2)
      r2=ddot_(ndensc,Kappa,1,Kappa,1)
      If(debug)Write(6,*) 'In that case I think that r2 should be:',r2
      If (r2.gt.r1) Write(6,*) 'Warning ',
     &    ' perturbation number ',idisp,' might diverge'
*
      call dcopy_(ndensC,Kappa,1,dKappa,1)
*
*      In MS-PDFT deltaC is no longer zero initially
*      deltaC=Zero
*      Use following two lines instead
      CALL DMinvCI_SA(ipST,W(ipS2)%Vec,rdum(1),isym,Fancy)
      call dcopy_(nConf1*nroots,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
      deltaC=ddot_(nConf1*nroots,W(ipST)%Vec,1,W(ipS2)%Vec,1)
      irc=ipOut(ipcid)
      deltaK=ddot_(nDensC,Kappa,1,Sigma,1)
      Kappa(1:nDens)=Zero
      delta=deltac+deltaK
      delta0=delta
      iter=1
      If (delta.eq.Zero) Goto 300
*-----------------------------------------------------------------------
*
200   Continue
         Call TimesE2MSPDFT(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

*
*-----------------------------------------------------------------------
*
*                   delta
*        rAlpha=------------
*               dKappa:dSigma
*
*-----------------------------------------------------------------------
*
         rAlphaK=Zero
         rAlphaK=ddot_(nDensC,Temp4,1,dKappa,1)
         rAlphaC=Zero
         irc=ipIn(ipS1)
         irc=ipIn(ipCId)
         rAlphaC=ddot_(nConf1*nroots,W(ipS1)%Vec,1,W(ipCId)%Vec,1)
         rAlpha=delta/(rAlphaK+rAlphaC)
*
*-------------------------------------------------------------------*
*
*        Kappa=Kappa+rAlpha*dKappa
         Call DaxPy_(nDensC,ralpha,dKappa,1,Kappa,1)
*        Sigma=Sigma-rAlpha*dSigma       Sigma=RHS-Akappa
         Call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
         resk=sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
         resci=Zero
         irc=ipIn(ipCIT)
         Call DaXpY_(nConf1*nroots,ralpha,W(ipCId)%Vec,1,W(ipCIT)%Vec,1)
         irc=ipOut(ipcit)
*        ipST =ipST -rAlpha*ipS1         ipST=RHS-A*ipCIT
         irc=ipIn(ipS1)
         irc=ipIn(ipST)
         Call DaXpY_(nConf1*nroots,-ralpha,W(ipS1)%Vec,1,W(ipST)%Vec,1)
         irc=opOut(ipS1)
         irc=ipIn(ipST)
         resci=sqrt(ddot_(nconf1*nroots,W(ipST)%Vec,1,W(ipST)%Vec,1))

*-------------------------------------------------------------------*
*
*        Precondition......
*           -1
*        S=M  Sigma
*
         irc=opOut(ipcid)

         irc=ipIn(ipS2)
         Call DMinvCI_SA(ipST,W(ipS2)%Vec,rdum(1),isym,Fancy)
         irc=opOut(ipci)
         irc=opOut(ipdia)

         irc=ipIn(ipPre2)
         Call DMInvKap(W(ipPre2)%Vec,Sigma,nDens2+6,
     &                 Sc2,nDens2+6,Sc1,nDens2+6,
     &                 iSym,iter)
         irc=opOut(ippre2)
*
*-------------------------------------------------------------------*
*             s:Sigma (k+1)     s:Sigma (k+1)
*        Beta=-------        =  -------------
*              delta  (k)        s:Sigma (k)
*
*        delta=s:sigma
*
*        dKappa=s+Beta*dKappa
*
         deltaC=ddot_(nConf1*nroots,W(ipST)%Vec,1,W(ipS2)%Vec,1)
         irc=ipOut(ipST)
*
         deltaK=ddot_(nDensC,Sigma,1,Sc2,1)
         If (.not.CI) Then
            rBeta=deltaK/delta
            delta=deltaK
            Call DScal_(nDensC,rBeta,dKappa,1)
            Call DaXpY_(nDensC,One,Sc2,1,dKappa,1)
         Else
            rbeta=(deltac+deltaK)/delta
            delta=deltac+deltaK
            irc=ipIn(ipCID)
            Call DScal_(nConf1*nroots,rBeta,W(ipCID)%Vec,1)
            Call DScal_(nDensC,rBeta,dKappa,1)
            Call DaXpY_(nConf1*nroots,One,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
            Call DaXpY_(nDensC,One,Sc2,1,dKappa,1)
            irc=opOut(ipS2)
            irc=ipOut(ipCID)
         End If

*    ######  #    #  #####        #####    ####    ####
*    #       ##   #  #    #       #    #  #    #  #    #
*    #####   # #  #  #    #       #    #  #       #
*    #       #  # #  #    #       #####   #       #  ###
*    #       #   ##  #    #       #       #    #  #    #
*    ######  #    #  #####        #        ####    ####
*
*-------------------------------------------------------------------*
*
*
         res=Zero ! dummy initialize
         If (iBreak.eq.1) Then
            If (abs(delta).lt.abs(Epsilon**2*delta0)) Goto 300
         Else If (iBreak.eq.2) Then
            res=sqrt(resk**2+resci**2+ResQaaLag2)
            if (doDMRG) res=sqrt(resk**2)
            If (res.lt.abs(epsilon)) Goto 300
         Else
            If (abs(delta).lt.abs(Epsilon**2*delta0).and.
     &          res.lt.abs(epsilon))  Goto 300
         End If
         If (iter.ge.niter) goto 210
         If (lprint)
     &   Write(6,Fmt2//'I7,4X,ES17.9,ES17.9,ES17.9,ES17.9,ES17.9)')
     &          iter,delta/delta0,resk,resci,deltac,deltak
         iter=iter+1
*
         Goto 200
*
************************************************************************
*
 210     Continue
         Write(6,Fmt2//'A,I4,A)')
     &         'No convergence for perturbation no: ',
     &          idisp,'. Increase Iter.'
         converged(isym)=.false.
         fail=.true.
         Goto 310
 300  Continue
      If (iPL.ge.2) Then
          Write(6,Fmt2//'A,I4,A,I4,A)')
     &          'Perturbation no: ',idisp,' converged in ',
     &          iter-1,' steps.'
      End If
      irc=ipnout(-1)
*
 310  Continue
      If (iPL.ge.2) Write(6,*)
      iLen=ndensC
      iKapDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Kappa,iLen,iDis)
      iSigDisp(iDisp)=iDis
      Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
      ilen=nconf1*nroots
      iCIDisp(iDisp)=iDis
*
      irc=ipin(ipCIT)
      Call dDaFile(LuTemp,1,W(ipCIT)%Vec,iLen,iDis)
*
**MGD This last call seems unused, so I comment it
*
*      Call TimesE2(Kappa,ipCIT,1,reco,jspin,ipS2,
*     &             Temp4,ipS2)
      iCISigDisp(iDisp)=iDis
      irc=ipin(ipST)
      Call dDaFile(LuTemp,1,W(ipST)%Vec,iLen,iDis)
      end do ! iDisp
*
      Call mma_deallocate(Temp4)
      Call mma_deallocate(Temp3)
      Call mma_deallocate(Sigma)
      Call mma_deallocate(dKappa)
      Call mma_deallocate(Kappa)
      Call mma_deallocate(Sc2)
      Call mma_deallocate(Sc1)
*
*     Free all memory and remove from disk all data
*     related to this symmetry
*
193   Continue
      Call mma_deallocate(Fancy)

      irc=ipclose(ipdia)
      If (.not.CI) irc=ipclose(ipPre2)
*
      Call Exp_Close()

      If (debug) Then
      Write(6,*)  '****************************************'//
     &            '****************************************'
      Write(6,*)
      End If
      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      End

      Subroutine TimesE2MSPDFT(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,
     & ipCiOut)
      use ipPage, only: w
      Implicit Real*8(a-h,o-z)
#include "stdalloc.fh"
#include "Pointers.fh"
#include "dmrginfo_mclr.fh"
#include "real.fh"
#include "Input.fh"
      Integer opOut
      Real*8 Kap(*),KapOut(*)
      Real*8 rdum(1)
      Real*8, Allocatable:: Temp3(:), Temp4(:),
     &                      Sc1(:), Sc2(:), Sc3(:), RMOAA(:),MSHam(:)
      INTEGER kRoot,lRoot
      Real*8 ECOff
*
      Call mma_allocate(RMOAA,n2Dens,Label='RMOAA')
      Call mma_allocate(Sc1,nDens2,Label='Sc1')
      Call mma_allocate(Sc2,nDens2,Label='Sc2')
      Call mma_allocate(Sc3,nDens2,Label='Sc3')
      Call mma_allocate(Temp3,nDens2,Label='Temp3')
      Call mma_allocate(Temp4,nDens2,Label='Temp4')
*
      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
      Call Uncompress(Kap,Sc1,isym)

! Integral derivative !yma
      Call RInt_generic(SC1,rmoaa,rdum,
     &                 Sc2,
     &                 Temp3,Temp4,Sc3,
     &                 isym,reco,jspin)

      Call Kap_CI(Temp4,nDens2,rmoaa,n2Dens,ipCIOUT)
      Call Ci_Ci(ipcid,ipS2)
      Call CI_KAP(ipCid,Sc1,Sc3,isym)

      Call DZaXpY(nDens,One,Sc2,1,Sc3,1,Sc1,1)
*
      Call Compress(Sc1,KapOut,isym)   ! ds
*     Call RecPrt('Ex',' ',KapOut,ndensC,1)
*
      irc=ipin(ipS2)
      irc=ipin(ipCIOUT)
      Call DaXpY_(nConf1*nroots,One,
     &               W(ipS2)%Vec,1,
     &               W(ipCIOUT)%Vec,1)
      irc=opOut(ipCId)


******Adding Hkl contribution
      Call mma_allocate(MSHam,nRoots**2)
      CALL CMSRdMat(MSHam,nRoots,nRoots,'ROT_HAM',7)
      DO kRoot=1,nRoots
       Do lRoot=1,nRoots
        IF (kRoot.eq.lRoot) Cycle
        ECOff=-MSHam(nRoots*(kRoot-1)+lRoot)*Weight(1)*2.0d0
        Call DaXpY_(nConf1,ECOff,
     &              W(ipCId)%Vec((lRoot-1)*nConf1+1),1,
     &            W(ipCIOut)%Vec((kRoot-1)*nConf1+1),1)
       End Do
      END DO
      Call mma_deallocate(MSHam)
******End of adding Hkl contribution


*
      Call mma_deallocate(Temp4)
      Call mma_deallocate(Temp3)
      Call mma_deallocate(Sc3)
      Call mma_deallocate(Sc2)
      Call mma_deallocate(Sc1)
      Call mma_deallocate(rmoaa)

      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
*
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      End
