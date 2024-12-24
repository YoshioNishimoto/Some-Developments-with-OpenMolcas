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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      SubRoutine WfCtl_SA(iKapDisp,iSigDisp,iCIDisp,iCIsigDisp,
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
      use pcm_grad, only: do_RF, PCM_grad_CLag, iStpPCM, def_solv,
     *                    PCM_mod_ERASSCF, ISRot, PCM_grad_PT2
      use ISRotation, only: CGS,DMInvISR,InvSCF,ISR,ISR_final,ISR_init,
     *                      ISR_projection,ISR_RHS,InvEne
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
     &                      Sc1(:), Sc2(:), Fancy(:),
     &                      SLag(:), wrk(:)
      Real*8, Allocatable :: R0(:),Pvec(:),Qvec(:),Uvec(:)
*
      interface
        subroutine RHS_NAC(Fock,SLag_pt2,ipS2)
          Real*8 Fock(*)
          real*8, optional :: SLag_pt2(*)
          integer, optional :: ipS2
        end subroutine
      end interface

      interface
        subroutine rhs_sa(Fock,SLag_pt2,ipS2)
          Real*8 Fock(*)
          real*8, optional :: SLag_pt2(*)
          integer, optional :: ipS2
        end subroutine
      end interface

*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
      Call StatusLine(' MCLR:',
     &                ' Computing Lagrangian multipliers for SA-CASSCF')
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
         Write (6,*) 'WfCtl_SA: SAVE option not implemented'
         Call Abend()
      End If
      If (iAnd(kprint,2).eq.2) lprint=.true.
      isym=1
      nconf1=ncsf(State_Sym)
      !! iStpPCM has been set somewhere, but just as a mnemonic device
      iStpPCM = 1

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
      !! Check whether we need to consider internal state rotations
      !! iteratively or non-interatively
      InvEne = .true.
      if (.not.InvSCF) InvEne = .false.
      if (PT2 .and. nRoots > 0) InvEne = .false.
      if (do_RF .and. def_solv /= 3) InvEne = .false.
C       InvSCF=.false.
C     if (do_RF) InvEne = .false.
      !! add unequal state-averaging
*
      !! initialize CGS and some for InvEne
      call ISR_Init(iPL,nRoots,ncsf,State_Sym,
     *              ERASSCF)
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
      if (CGS) then
        ipPvec=ipGet(nconf1*nroots)
        ipQvec=ipGet(nconf1*nroots)
        ipUvec=ipGet(nconf1*nroots)
        ipR0  =ipGet(nconf1*nroots)
      end if
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
      if (CGS) then
        Call mma_allocate(Pvec,nDens2+6,Label='Pvec')
        Call mma_allocate(Qvec,nDens2+6,Label='Qvec')
        Call mma_allocate(Uvec,nDens2+6,Label='Uvec')
        Call mma_allocate(R0,nDens2+6,Label='R0')
      end if
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

      Call mma_allocate(SLag,nRoots**2,Label='SLag')
      SLag(1:nRoots**2) = Zero
      Call DCopy_(nconf1*nRoots,[Zero],0,W(ipST)%Vec,1)
      if (PT2) then
        Call RHS_PT2(Kappa,W(ipST)%Vec,SLag)
        write (6,*) "SLAG from PT2"
        call sqprt(slag,nroots)
        !! SLag given by the CASPT2 routine should be treated as the
        !! initial rotation of the reference states.
        !! If we need to treat parameters as internal rotations, they
        !! will be computed later using configuration parametrs.
C       if (.not.InvSCF) then
C         call dcopy_(nRoots**2,SLag,1,ISR%p,1)
C         SLag = Zero
C       end if
        if (do_RF) call PCM_grad_PT2()
      end if
*
      If (do_RF) Then
        !! SA-CASSCF/PCM is not stationary wrt CI coefficients
        !! Use ipS2 for the moment (some intermediate vectors are SD)
C       if (PT2) SLag(1:nRoots**2) = Zero
        Call PCM_grad_CLag(1,ipCI,ipS2,SLag)
        if (InvSCF) call daxpy_(nRoots**2,One,ISR%Rvec,1,SLag,1)
        Call DCopy_(nconf1*nRoots,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
        Call DCopy_(nconf1*nRoots,[Zero],0,W(ipS2)%Vec,1)
C       if (.not.NonInv) call daxpy_(nRoots**2,+One,SLag,1,ISR%p,1)
      End If
C     if (.not.NonInv) call dcopy_(nRoots**2,ISR%p,1,SLag,1)
*
      If (isNAC) Then
        Call RHS_NAC(Temp4,SLag,ipS2)
      Else
        Call RHS_SA(Temp4,SLag,ipS2)
      End If
      if (allocated(SLag)) Call mma_deallocate(SLag)
*
      !! add CI derivative contributions (if evaluated in rhs_sa/nac)
      !! due to the non-iterative internal state rotations
C     if (do_RF .and. (.not.InvEne .and. InvSCF)) then
      If (do_RF .and. ((PT2 .and. nRoots > 1)
     *            .or. (.not.InvEne .and. InvSCF))) then
        Call Daxpy_(nconf1*nRoots,One,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
C       Call DScal_(nconf1*nRoots,Two,W(ipCID)%Vec,1)
C       if (def_solv /= 3) then
C        Call daxpy_(nconf1*nRoots,One,W(ipS2)%Vec,1,W(ipCID)%Vec,1)
C      Call DCopy_(nconf1*nRoots,W(ipST)%Vec,1,W(ipCID)%Vec,1)
C         Call DCopy_(nconf1*nRoots,[Zero],0,W(ipST)%Vec,1)
C       end if
      End If
C       write (6,*) "Temp4"
C       call sqprt(temp4,nbas(1))
C       write (6,*) "ipST initial"
C       do i = 1, nRoots
C         write (6,*) "iroot = ", i
C         do j = 1, nconf1
C          write (6,'(i4,2f20.10)') j,w(ipST)%Vec(j+nconf1*(i-1)),
C    *       w(ipCID)%Vec(j+nconf1*(i-1))
C         end do
C       end do
      irc=ipOut(ipST)
*
      If (PT2) Then
        Call DaXpY_(nDens2,1.0D+00,Kappa,1,Temp4,1)
        Kappa(1:nDens2)=Zero
      End If
      irc=opOut(ipci)
*
      If (lprint) Write(6,*)
     &      '       Iteration       Delta           Res(kappa)     '//
     &      '  Res(CI)          DeltaK           DeltaC'
      iLen=nDensC
      iRHSDisp(iDisp)=iDis
      Call Compress(Temp4,Sigma,iSym)
      r1=ddot_(nDensc,Sigma,1,Sigma,1)
      If (PT2)
     *  R1 = R1 + DDot_(nConf1*nRoots,W(ipST)%Vec,1,W(ipST)%Vec,1)
      If (do_RF)
     *  R1 = R1 + DDot_(nConf1*nRoots,W(ipCID)%Vec,1,W(ipCID)%Vec,1)
      If(debug)Write(6,*) 'Hi how about r1',r1
      Call dDaFile(LuTemp,1,Sigma,iLen,iDis)
*
      If (PT2.or.do_RF) then
        Call DSCAL_(nConf1*nRoots,-One,W(ipST)%Vec,1)
        If (CI) Then
          !! The order of CSF coefficients in CASPT2 and MCLR is somehow
          !! different, so the CI lagrangian computed in CASPT2 must be
          !! reordered so that it can be used here.
          If (PT2) Then
            Call mma_allocate(wrk, nConf1, Label='wrk')
            Do iR = 1, nRoots
              Call DCopy_(nConf1,W(ipST)%Vec(1+nConf1*(iR-1):nConf1*iR),
     *                    1,wrk,1)
              Call GugaCtl_MCLR(wrk,1)
              Call DCopy_(nConf1,wrk,1,
     *                    W(ipST)%Vec(1+nConf1*(iR-1):nConf1*iR),1)
            End Do
            Call mma_deallocate(wrk)
          End If
C             Call ISR_RHS(W(ipCI)%Vec,W(ipST)%Vec)
C             isr%rvec = 0.0d+00
C             call dcopy_(nconf1*nroots,[zero],0,W(ipST)%vec,1)
*
C           write (*,*) "first irs_rhs"
          ! Call ISR_RHS(W(ipCI)%Vec,W(ipST)%Vec)
          ! Call ISR_Projection(W(ipCI)%Vec,W(ipST)%Vec)
          ! isr%rvec = isr%rvec * 2.0d+00
          ! do i = 1, nroots
          !   isr%rvec(i,i) = isr%rvec(i,i)*0.5d+00
          ! end do
          If (do_RF) Then
            !! Yes, it should be multiplied by two
            Call Daxpy_(nConf1*nRoots,-One,W(ipCID)%Vec,1,W(ipST)%Vec,1)
            call dcopy_(nConf1*nroots,[Zero],0,W(ipCID)%Vec,1)
          End If
*
          !! Now, compute the RHS of the internal state rotation
          if (.not.InvSCF) then
            write (6,*) "second irs_rhs"
            Call ISR_RHS(W(ipCI)%Vec,W(ipST)%Vec)
C           write (6,*) "isr%p to be added"
C           call sqprt(isr%p,nroots)
            if (PT2) ISR%Rvec = ISR%Rvec + ISR%p*0.5d+00
            write (6,*) "initial state rotation + PT2 rotation"
            call sqprt(isr%rvec,nroots)
          end if
          Call ISR_Projection(W(ipCI)%Vec,W(ipST)%Vec)
*
          !! precondition (z0 = M^{-1}*r0)
          if (CGS) then
            Call dcopy_(nConf1*nRoots,W(ipST)%Vec,1,W(ipS2)%Vec,1)
          else
            Call DMinvCI_sa(ipST,W(ipS2)%Vec,rdum(1),isym,fancy)
          end if
          irc=opOut(ipci)
          irc=opOut(ipdia)
          !! z0 <= p0
          Call DCopy_(nConf1*nRoots,W(ipS2)%Vec,1,
     *                              W(ipCId)%Vec,1)
          irc=ipIn(ipCIT)
          call dcopy_(nConf1*nroots,[Zero],0,W(ipCIT)%Vec,1)
        End If
C       write (6,*) "ipST initial"
C       do i = 1, nRoots
C         write (6,*) "iroot = ", i
C         do j = 1, nconf1
C          write (6,'(i4,2f20.10)') j,w(ipST)%Vec(j+nconf1*(i-1)),
C    *       w(ipCID)%Vec(j+nconf1*(i-1))
C         end do
C       end do
      Else
        irc=ipIn(ipCIT)
        irc=ipIn(ipST)
        irc=ipIn(ipCID)
        call dcopy_(nConf1*nroots,[Zero],0,W(ipCIT)%Vec,1)
        call dcopy_(nConf1*nroots,[Zero],0,W(ipST)%Vec,1)
        call dcopy_(nConf1*nroots,[Zero],0,W(ipCID)%Vec,1)
      End If
      irc=ipOut(ipCIT)
*
      if (.not.InvSCF.and..not.CGS) then
        ISR%Rvec = +ISR%Rvec
        Call DMInvISR(ISR%Rvec,ISR%prec)
        ISR%p = ISR%prec
C
C       ISRot(:,:,3) = -ISRot(:,:,1)
C       Call DMInvISR(ISRot(:,:,3),ISRot(:,:,4))
C       ISRot(:,:,2) = ISRot(:,:,4)
        write (6,*) "isrot 4 (0)"
C       call sqprt(isrot(:,:,4),nroots)
        call sqprt(ISR%prec,nroots)
        write (6,*) "isrot 3 (0)"
C       call sqprt(isrot(:,:,3),nroots)
        call sqprt(ISR%Rvec,nroots)
      end if
*
      Call DSCAL_(nDensC,-One,Sigma,1)
*
      if (CGS) then
        call dcopy_(nDensC,Sigma,1,Kappa,1)
      else
        irc=ipIn(ipPre2)
        Call DMInvKap(W(ipPre2)%Vec,Sigma,nDens2+6,
     &                Kappa,nDens2+6,Temp3,nDens2+6,
     &                isym,iter)
        irc=opOut(ippre2)
      end if
*
      r2=ddot_(ndensc,Kappa,1,Kappa,1)
      If (PT2.or.do_RF)
     *  R2 = R2 + DDot_(nConf1*nRoots,W(ipS2)%Vec,1,W(ipS2)%Vec,1)
      If(debug)Write(6,*) 'In that case I think that r2 should be:',r2
      If (r2.gt.r1) Write(6,*) 'Warning ',
     &    ' perturbation number ',idisp,' might diverge'
*
      call dcopy_(ndensC,Kappa,1,dKappa,1)
*
      if (CGS) then
        ! R0
        call dcopy_(nDensC,Kappa,1,R0,1)
        call dcopy_(nConf1*nRoots,W(ipST)%Vec,1,W(ipR0)%Vec,1)
        if (.not.InvSCF) ISR%R0 = +ISR%Rvec
C       if (.not.InvSCF) ISRot(:,:,9) = -ISRot(:,:,1)

        ! u_k
        call dcopy_(nDensC,Kappa,1,Uvec,1)
        call dcopy_(nConf1*nRoots,W(ipST)%Vec,1,W(ipUvec)%Vec,1)
        if (.not.InvSCF) ISR%Uvec = +ISR%Rvec
C       if (.not.InvSCF) ISRot(:,:,8) = -ISRot(:,:,1)

        ! p_k
        call dcopy_(nDensC,Uvec,1,dKappa,1)
        Call DCopy_(nDensC,dKappa,1,Pvec,1)
        call dcopy_(nConf1*nRoots,W(ipUvec)%Vec,1,W(ipCId)%Vec,1)
        Call DCopy_(nconf1*nRoots,W(ipCId)%Vec,1,W(ipPvec)%Vec,1)
        if (.not.InvSCF) then
          ISR%Pvec = ISR%Uvec
          ISR%Rvec = ISR%R0
          ISR%Xvec = 0.0d+00
        end if
C       if (.not.InvSCF) ISRot(:,:,6) = ISRot(:,:,8)
C       if (.not.InvSCF) ISRot(:,:,3) = ISRot(:,:,9)
C       if (.not.InvSCF) ISRot(:,:,5) = 0.0d+00

        deltaC = ddot_(nConf1*nroots,W(ipR0)%Vec,1,W(ipR0)%Vec,1)
        if (.not.InvSCF) deltaC = deltaC
C    *    + ddot_(nRoots**2,ISRot(:,:,9),1,ISRot(:,:,9),1)
     *    + ddot_(nRoots**2,ISR%R0,1,ISR%R0,1)
        deltaK = ddot_(nDensC,R0,1,R0,1)
        Kappa(1:nDens)=Zero
        delta  = deltaC + deltaK
      else
        deltaC=Zero
        If (PT2.or.do_RF)
     &    deltaC=ddot_(nConf1*nroots,W(ipST)%Vec,1,W(ipS2)%Vec,1)
        if (.not.InvSCF) deltaC = deltaC
C    *    + ddot_(nRoots**2,ISRot(:,:,3),1,ISRot(:,:,4),1)
     *    + ddot_(nRoots**2,ISR%Rvec,1,ISR%prec,1)
        irc=ipOut(ipcid)
        deltaK=ddot_(nDensC,Kappa,1,Sigma,1)
        Kappa(1:nDens)=Zero
        delta=deltac+deltaK
      end if
      delta0=delta
      iter=1
      If (delta.eq.Zero) Goto 300
      iStpPCM = 2
*-----------------------------------------------------------------------
*
200   Continue
*
         if (CGS) then
           !! (Preconditioned) Conjugate Gradient Squared (CGS) method
           !! instead of the default PCG
           !! CGS is preferred, if internal rotations are needed (e.g.,
           !! PCM)
           !! The implementation is based on Algorithm 4 in
           !! "Preconditioned Algorithm of the CGS Method Focusing on
           !! Its Deriving Process" (Japanese literature)

           !! precondition p (K-1*p)
           irc=ipIn(ipS2)
           Call DMinvCI_SA(ipPvec,W(ipCId)%Vec,rdum(1),isym,Fancy)
           irc=opOut(ipci)
           irc=opOut(ipdia)

           irc=ipIn(ipPre2)
           Call DMInvKap(W(ipPre2)%Vec,Pvec,nDens2+6,
     &                   dKappa,nDens2+6,Sc1,nDens2+6,
     &                   iSym,iter)
           irc=opOut(ippre2)

C          if (.not.InvSCF) Call DMInvISR(ISRot(:,:,6),ISRot(:,:,2))
           if (.not.InvSCF) Call DMInvISR(ISR%Pvec,ISR%p)

           !! A*K-1*p
           Call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

           !! compute alpha
           rAlphaK=ddot_(nDensC,R0,1,Temp4,1)
           rAlphaC=ddot_(nConf1*nroots,W(ipR0)%Vec,1,W(ipS1)%Vec,1)
           if (.not.InvSCF) rAlphaC = rAlphaC
C    *       + ddot_(nRoots**2,ISRot(:,:,9),1,ISRot(:,:,1),1)
     *       + ddot_(nRoots**2,ISR%R0,1,ISR%Ap,1)
           rAlpha=delta/(rAlphaK+rAlphaC)

           ! q_k = u_k - alpha*(A*K-1*p)
           call dcopy_(nDensC,Uvec,1,Qvec,1)
           call daxpy_(nDensC,-ralpha,Temp4,1,Qvec,1)
           call dcopy_(nConf1*nRoots,W(ipUvec)%Vec,1,W(ipQvec)%Vec,1)
           call daxpy_(nConf1*nRoots,-ralpha,W(ipS1)%Vec,1,
     *                                       W(ipQvec)%Vec,1)
C          if (.not.InvSCF) then
C            ISRot(:,:,7) = ISRot(:,:,8) - ralpha*ISRot(:,:,1)
C          end if
           if (.not.InvSCF) ISR%Qvec = ISR%Uvec - ralpha*ISR%Ap

           !! precondition u + q
           call daxpy_(nDensC,+One,Qvec,1,Uvec,1)
           call daxpy_(nConf1*nRoots,+One,W(ipQvec)%Vec,1,
     *                                    W(ipUvec)%Vec,1)
C         if (.not.InvSCF)
C    *      call daxpy_(nRoots**2,+One,ISRot(:,:,7),1,ISRot(:,:,8),1)
           if (.not.InvSCF) ISR%Uvec = ISR%Uvec + ISR%Qvec

           irc=ipIn(ipS2)
           Call DMinvCI_SA(ipUvec,W(ipCId)%Vec,rdum(1),isym,Fancy)
           irc=opOut(ipci)
           irc=opOut(ipdia)

           irc=ipIn(ipPre2)
           Call DMInvKap(W(ipPre2)%Vec,Uvec,nDens2+6,
     &                   dKappa,nDens2+6,Sc1,nDens2+6,
     &                   iSym,iter)
           irc=opOut(ippre2)

C          if (.not.InvSCF) Call DMInvISR(ISRot(:,:,8),ISRot(:,:,2))
           if (.not.InvSCF) Call DMInvISR(ISR%Uvec,ISR%p)

           !! x_k+1 = x_k + alpha*K-1*(u+q)
           call daxpy_(nDensC,+ralpha,dKappa,1,Kappa,1)
           Call DaXpY_(nConf1*nroots,ralpha,W(ipCId)%Vec,1,
     *                                      W(ipCIT)%Vec,1)
C          if (.not.InvSCF) ISRot(:,:,5) = ISRot(:,:,5) + ralpha*ISRot(:,:,2)
           if (.not.InvSCF) ISR%Xvec = ISR%Xvec + ralpha*ISR%p

           !! A*K-1*(u+q)
           Call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

           !! r_k+1 = r_k - alpha*(A*K-1*(u+q))
           Call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
           resk=sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
           Call DaXpY_(nConf1*nroots,-ralpha,W(ipS1)%Vec,1,
     *                                       W(ipST)%Vec,1)
           resci=sqrt(ddot_(nconf1*nroots,W(ipST)%Vec,1,W(ipST)%Vec,1))
           if (.not.InvSCF) then
C            ISRot(:,:,3) = ISRot(:,:,3) - ralpha*ISRot(:,:,1)
C            resci = resci
C    *         + sqrt(ddot_(nRoots**2,ISRot(:,:,3),1,ISRot(:,:,3),1))
             ISR%Rvec = ISR%Rvec - ralpha*ISR%Ap
             resci = resci
     *         + sqrt(ddot_(nRoots**2,ISR%Rvec,1,ISR%Rvec,1))
           end if

           !! compute beta
           deltaK=ddot_(nDensC,R0,1,Sigma,1)
           deltaC=ddot_(nConf1*nroots,W(ipR0)%Vec,1,W(ipST)%Vec,1)
           if (.not.InvSCF) then
             deltaC = deltaC
C    *         + ddot_(nRoots**2,ISRot(:,:,9),1,ISRot(:,:,3),1)
     *         + ddot_(nRoots**2,ISR%R0,1,ISR%Rvec,1)
           end if
           rbeta=(deltaC+deltaK)/delta
           delta=deltac+deltaK

           !! u_k = r_k + beta*q
           call dcopy_(nDensC,Sigma,1,Uvec,1)
           call daxpy_(nDensC,+rbeta,Qvec,1,Uvec,1)
           call dcopy_(nConf1*nRoots,W(ipST)%Vec,1,W(ipUvec)%Vec,1)
           Call DaXpY_(nConf1*nroots,+rbeta,W(ipQvec)%Vec,1,
     *                                      W(ipUvec)%Vec,1)
C          if (.not.InvSCF) then
C            ISRot(:,:,8) = ISRot(:,:,3) + rbeta*ISRot(:,:,7)
C          end if
           if (.not.InvSCF) ISR%Uvec = ISR%Rvec + rbeta*ISR%Qvec

           !! p = u + beta*q + beta*beta*p
           call dscal_(nDensC,rbeta**2,Pvec,1)
           call daxpy_(nDensC,+One,Uvec,1,Pvec,1)
           call daxpy_(nDensC,+rbeta,Qvec,1,Pvec,1)
           call dscal_(nConf1*nRoots,rbeta**2,W(ipPvec)%Vec,1)
           Call DaXpY_(nConf1*nroots,+One,W(ipUvec)%Vec,1,
     *                                      W(ipPvec)%Vec,1)
           Call DaXpY_(nConf1*nroots,+rbeta,W(ipQvec)%Vec,1,
     *                                      W(ipPvec)%Vec,1)
           if (.not.InvSCF) then
C            ISRot(:,:,6) = ISRot(:,:,8) + rbeta*ISRot(:,:,7)
C    *                    + rbeta*rbeta*ISRot(:,:,6)
             ISR%Pvec = ISR%Uvec + rbeta*ISR%Qvec + rbeta*rbeta*ISR%Pvec
           end if
           !! end of CGS
         else
           !! PCG starts here

         Call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)
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
         if (.not.InvSCF) rAlphaC=rAlphaC
C    *     + ddot_(nRoots**2,ISRot(:,:,1),1,ISRot(:,:,2),1)
     *     + ddot_(nRoots**2,ISR%Ap,1,ISR%p,1)
         rAlpha=delta/(rAlphaK+rAlphaC)
         if (.not.InvSCF) then
         write (6,*) "isrot 1 and 2 after TimesE2"
C        call sqprt(isrot(:,:,1),nroots)
C        call sqprt(isrot(:,:,2),nroots)
         call sqprt(isr%Ap,nroots)
         call sqprt(isr%p,nroots)
         end if
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
         if (.not.InvSCF) then
         write (6,*) "ralpha = ",ralpha
C          ISRot(:,:,5) = ISRot(:,:,5) + ralpha*ISRot(:,:,2)
C          ISRot(:,:,3) = ISRot(:,:,3) - ralpha*ISRot(:,:,1)
C          ress=sqrt(ddot_(nRoots**2,ISRot(:,:,3),1,ISRot(:,:,3),1))
           ISR%Xvec = ISR%Xvec + ralpha*ISR%p
           ISR%Rvec = ISR%Rvec - ralpha*ISR%Ap
           ress=sqrt(ddot_(nRoots**2,ISR%Rvec,1,ISR%Rvec,1))
           write (6,*) "ress = ", ress
           resci = resci + ress
        write (6,*) "isrot 1 (1)"
C       call sqprt(isrot(:,:,1),nroots)
        call sqprt(isr%Ap,nroots)
        write (6,*) "isrot 3 (1)"
C       call sqprt(isrot(:,:,3),nroots)
        call sqprt(isr%Rvec,nroots)
         end if

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

C        if (.not.InvSCF) Call DMInvISR(ISRot(:,:,3),ISRot(:,:,4))
         if (.not.InvSCF) Call DMInvISR(ISR%Rvec,ISR%prec)
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
         if (.not.InvSCF) then
           deltaC = deltaC
C    *       + ddot_(nRoots**2,ISRot(:,:,3),1,ISRot(:,:,4),1)
     *       + ddot_(nRoots**2,ISR%Rvec,1,ISR%prec,1)
         end if
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
            if (.not.InvSCF) then
C             Call DScal_(nRoots**2,rBeta,ISRot(:,:,2),1)
C             Call DaXpY_(nRoots**2,One,ISRot(:,:,4),1,ISRot(:,:,2),1)
              ISR%p = ISR%prec + rbeta*ISR%p
        write (6,*) "isrot 4 (3)"
C       call sqprt(isrot(:,:,4),nroots)
        call sqprt(isr%prec,nroots)
        write (6,*) "isrot 2 (3)"
C       call sqprt(isrot(:,:,2),nroots)
        call sqprt(isr%p,nroots)
            end if
         End If
         end if ! CGS vs PCG

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
            res=sqrt(resk**2+resci**2)
            if (doDMRG) res=sqrt(resk**2)
            If (res.lt.abs(epsilon)) Goto 300
         Else
            If (abs(delta).lt.abs(Epsilon**2*delta0).and.
     &          res.lt.abs(epsilon))  Goto 300
         End If
         If (iter.ge.niter) goto 210
         If (lprint)
     &   Write(6,Fmt2//'I7,4X,ES17.9,ES17.9,ES17.9,ES17.9,ES17.9)')
     &          iter,delta/delta0,resk,resci,deltak,deltac
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
*
      !! Save the total internal rotations so that CIDens_sa can use
      !! the parameters to generate correct density
C     If (do_RF .and. .not.InvSCF) ISRot(:,:,2) = +ISRot(:,:,5)
      If (.not.InvSCF) ISR%p = ISR%Xvec
      end do ! iDisp
       Call TimesE2(Kappa,ipCIT,1,reco,jspin,ipS2,
     *              Temp4,ipS1)
C      write (6,*) "temp4"
C     do i = 1, ndensc
C     write (6,*) i,temp4(i)
C     end do
C      write (6,*) "ipS2"
C     do i = 1, nconf1*nroots
C     write (6,*) i,W(ips1)%Vec(i)
C     end do
C     if (.not.InvSCF) then
C     write (6,*) "isrot solution"
C     call sqprt(isr%xvec,nroots)
C     write (6,*) "isrot"
C     call sqprt(isr%ap,nroots)
C     end if
*
      Call mma_deallocate(Temp4)
      Call mma_deallocate(Temp3)
      Call mma_deallocate(Sigma)
      Call mma_deallocate(dKappa)
      Call mma_deallocate(Kappa)
      Call mma_deallocate(Sc2)
      Call mma_deallocate(Sc1)
      if (CGS) then
        Call mma_deallocate(Pvec)
        Call mma_deallocate(Qvec)
        Call mma_deallocate(Uvec)
        Call mma_deallocate(R0)
      end if
*
      call ISR_final()
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
      iStpPCM = 3
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

      Subroutine TimesE2(Kap,ipCId,isym,reco,jspin,ipS2,KapOut,ipCiOut)
      use ipPage, only: w
      use PCM_grad, only: do_RF,PCM_grad_TimesE2,W_SOLV
      use ISRotation, only: CGS,InvSCF,ISR,ISR_TimesE2
      use Arrays, only: CMO,g1t
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
     &                      Sc1(:), Sc2(:), Sc3(:), RMOAA(:)
*
      Call mma_allocate(RMOAA,n2Dens,Label='RMOAA')
      Call mma_allocate(Sc1,nDens2,Label='Sc1')
      Sc1(:)=Zero
      Call mma_allocate(Sc2,nDens2,Label='Sc2')
      Call mma_allocate(Sc3,nDens2,Label='Sc3')
      Call mma_allocate(Temp3,nDens2,Label='Temp3')
      Call mma_allocate(Temp4,nDens2,Label='Temp4')
      If (ActRot) Then
      Call GetMem('SCR5','ALLO','REAL',ipTemp5,ntash*ntash)
      End If
*
      if(doDMRG)then ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
      Call Uncompress(Kap,Sc1,isym)

! Integral derivative !yma
      Call RInt_generic(SC1,rmoaa,rdum(1),
     &                 Sc2,
     &                 Temp3,Temp4,Sc3,
     &                 isym,reco,jspin)

      Call Kap_CI(Temp4,nDens2,rmoaa,n2Dens,ipCIOUT)
      Call Ci_Ci(ipcid,ipS2)
      Call CI_KAP(ipCid,Sc1,Sc3,isym)
      if (do_RF) then
        Call DaXpY_(nConf1*nroots,One,
     &              W(ipS2)%Vec,1,W(ipCIOUT)%Vec,1)
        Call Uncompress(Kap,Sc1,isym)
        if (CGS .or. .not.InvSCF) ISR%Ap = 0.0d+00
        call PCM_grad_TimesE2(1,isym,Sc1,Sc3,ipS2) ! ipS2 is overwritten
C         call PCM_grad_projection(1,W(ipCI)%Vec,W(ipS2)%Vec,Temp4,
C    *                             Temp4,Temp4)
C         write (6,*) "ipCIOUT from PCM cont."
C         do i = 1, nRoots
C           write (6,*) "iroot = ", i
C           do j = 1, nconf1
C            write (6,'(i4,f20.10)') j,w(ipS2)%Vec(j+nconf1*(i-1))
C           end do
C         end do
*
        !! evaluate the diagonal S-S block if SCF is not invariant
        !! Also, scale the CI derivative contributions
        Call ISR_TimesE2(W(ipCI)%Vec,W(ipS2)%Vec,W_SOLV)
      end if

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
*
      !if (do_RF) then
      !  !! Temp4 is used as a working array
      !  if (.not.InvSCF) then
C     !    write (6,*) "ipCIOUT before projection"
C     !    do i = 1, nRoots
C     !      write (6,*) "iroot = ", i
C     !      do j = 1, nconf1
C     !       write (6,'(i4,f20.10)') j,w(ipCIOUT)%Vec(j+nconf1*(i-1))
C     !      end do
C     !    end do
      !    Call
      !  ! call PCM_grad_projection(1,W(ipCI)%Vec,W(ipCIOUT)%Vec,Temp4,
     *!  !                          Temp4,Temp4)
      !    call PCM_grad_projection(0,W(ipCI)%Vec,W(ipCIOUT)%Vec,Temp4,
     *!                             Temp4,Temp4)
C     !    write (6,*) "ipCIOUT after projection"
C     !    do i = 1, nRoots
C     !      write (6,*) "iroot = ", i
C     !      do j = 1, nconf1
C     !       write (6,'(i4,f20.10)') j,w(ipCIOUT)%Vec(j+nconf1*(i-1))
C     !      end do
C     !    end do
      !  else
      !    call PCM_grad_projection(0,W(ipCI)%Vec,W(ipCIOUT)%Vec,Temp4,
     *!                             Temp4,Temp4)
      !  end if
      !else
      !! Project out the internal rotation contribution
      do iR = 1, nroots
        do jR = 1, nroots
          ovl = ddot_(nconf1,W(ipCIOUT)%Vec(1+(iR-1)*nconf1),1,
     *                       W(ipCI)%Vec(1+(jR-1)*nconf1),1)
          call daxpy_(nconf1,-ovl,W(ipCI)%Vec(1+(jR-1)*nconf1),1,
     *                            W(ipCIOUT)%Vec(1+(iR-1)*nconf1),1)
        end do
      end do
*
      Call mma_deallocate(Temp4)
      Call mma_deallocate(Temp3)
      Call mma_deallocate(Sc3)
      Call mma_deallocate(Sc2)
      Call mma_deallocate(Sc1)
      Call mma_deallocate(rmoaa)
      If (ActRot) Then
      Call GetMem('SCR5','FREE','REAL',ipTemp5,ntash*ntash)
      End If

      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
*
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      End
