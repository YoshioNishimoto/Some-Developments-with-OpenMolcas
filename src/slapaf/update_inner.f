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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      Subroutine Update_inner(kIter,Beta,Beta_Disp,Step_Trunc,nWndw,
     &                        mIter,Kriging_Hessian,qBeta,iOpt_RS,
     &                        First_MicroIteration,Iter,qBeta_Disp)
************************************************************************
*     Object: to update coordinates                                    *
*                                                                      *
*    Input:                                                            *
*      kIter          : iteration counter                              *
*      Beta           : damping factor step length                     *
*      Beta_Disp      : damping factor variance                        *
*                                                                      *
*    OutPut:                                                           *
*      Step_Trunc     : character label to denote truncation type      *
*                                                                      *
*                                                                      *
*     Author: Roland Lindh                                             *
*             2000                                                     *
************************************************************************
      use Slapaf_info, only: GNrm, Lambda, Energy, qInt, dqInt, Shift,
     &                       BMx, Degen, nStab, Smmtrc, Lbl
      use Slapaf_Parameters, only: iRow_c, iInt, nFix, iOptH,
     &                             HrmFrq_Show, Curvilinear, FindTS,
     &                             nBVec, nDimBC, iOptC, iNeg,
     &                             TSConstraints, GNrm_Threshold, Mode,
     &                             GrdMax, StpMax, nLambda
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Integer iDum(1)
      Logical Found, Kriging_Hessian, First_MicroIteration
      Character Step_Trunc, File1*8, File2*8, Step_Trunc_
      Real*8, Allocatable:: Hessian(:,:), Wess(:,:), AMat(:),
     &                      RHS(:), ErrVec(:,:), EMtrx(:,:)
      Integer, Allocatable:: Index(:), iFlip(:)
      Real*8, Allocatable:: R(:,:), dRdq(:,:,:), CInt(:), CInt0(:),
     &                      T(:,:), d2L(:,:,:), BM(:,:), dBM(:,:,:),
     &                      Energy_L(:), dEdx(:,:), x(:,:), du(:),
     &                      dEdq(:,:), dx(:,:), dy(:), BVec(:,:),
     &                      Scr1(:), Scr2(:), Value(:), Value0(:),
     &                      Mult(:), dBVec(:), Tmp(:)
      Character(LEN=8), Allocatable:: Lbl_Tmp(:)
*                                                                      *
************************************************************************
*                                                                      *
* This option should careful debugged more.
*                                              RL December 2020
#define _NEW_CODE_
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _NEW_CODE_
      Real*8, Allocatable:: QC(:,:,:)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      nQQ = SIZE(qInt,1)
      nsAtom=SIZE(Degen,2)
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Lu=6
      Write (Lu,*)'Update_inner:iOpt_RS,Beta,Beta_Disp=',
     &                     iOpt_RS,Beta,Beta_Disp
      Call RecPrt('Update_inner: qInt',' ',qInt,nQQ,kIter)
      Call RecPrt('Update_inner: Shift',' ',Shift,nQQ,kIter-1)
      Call RecPrt('Update_inner: GNrm',' ',GNrm,kIter,1)
      Call RecPrt('Update_inner: Energy',' ',Energy,kIter,1)
      n1=3*nsAtom
      Call RecPrt('Update_inner: dQ/dx(BMx)',' ',BMx,n1,nQQ)
#endif
*
      GrdMax=Zero
      StpMax=Zero
*
      Step_Trunc='N'
      nA = (Max(nQQ,kIter)+1)**2
      Call mma_Allocate(AMat,nA,Label='AMat')
      Call mma_Allocate(Index,kIter,Label='Index')
      Call mma_Allocate(ErrVec,nQQ,(kIter+1),Label='ErrVec')
      Call mma_Allocate(EMtrx,kIter+1,kIter+1,Label='EMtrx')
      Call mma_Allocate(RHS,kIter+1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Pick up the Hessian in internal coordinates and then update it
*     according to some Hessian update method (BFGS, MSP, etc.)
*
*
      If (Kriging_Hessian) Then
         iOptH_ = iOr(8,iAnd(iOptH,32))
      Else
         iOptH_ = iOptH
      End If
      Call mma_Allocate(Hessian,nQQ,nQQ,Label='Hessian')
      Call Get_dArray('Hss_Q',Hessian,nQQ**2)
*
*     Perform the Hessian update, in case of GEK it essentially will
*     modify the Hessian if it is needed to guide 2nd order
*     optimization towards a minimum or a TS.
*
#ifdef _DEBUGPRINT_
      Write (Lu,*)
      Write (Lu,*)
      Write (Lu,*) ' *** Updating the molecular Hessian ***'
      Write (Lu,*)
      jPrint=99
#else
      jPrint=5
#endif
      Call Update_H(nWndw,Hessian,nQQ,
     &              mIter,iOptC,
     &              Shift(1,kIter-mIter+1),dqInt(1,kIter-mIter+1),
     &              iOptH_,jPrint,GNrm(kIter),
     &              nsAtom,.True.,
     &              First_MicroIteration)
*
*     Call RecPrt('Update_inner: Hessian',' ',Hessian,nQQ,nQQ)
*     Write (6,*) 'After corrections'
*     Call DiagMtrx(Hessian,nQQ,jNeg)
*
*     Save the number of internal coordinates on the runfile.
*
      Call Put_iScalar('No of Internal coordinates',nQQ)
*
*     Save the force constant matrix on the runfile.
*
      Call Put_dArray('Hess',Hessian,nQQ**2)
*
*     Optional frequency analysis
*
      If (HrmFrq_Show) Call GF_on_the_fly(iDo_DipM)

*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---- If this is a constrained optimization with the goal of finding a
*     transition state and we have one negative eigenvalue of the
*     Hessian then shift over to Mode Following RF to find the TS.
*
*     If TSConstraints were given, merge them with global constraints if
*     they exist, unless we are in TS regime.
*     If TS constraints were not given, remove global constraints when in
*     TS regime.
*
      Call qpg_darray('TanVec',Found,nRP)
      If (FindTS.and.First_MicroIteration) Then
         File1='UDC'
         File2='TSC'
         If (.not.TSConstraints) File2=''
         If (iNeg(1).ge.1) Then
            If ((GNrm(kIter).le.GNrm_Threshold).or.Found) Then
*              Change to MFRF optimization.
               Mask=1+2+4+8+16+32+64+256+512+1024+2048+8192
               iOptC=iAnd(iOptC,Mask)
               Call Put_lScalar('TS Search',.True.)
               If (TSConstraints) Then
                  File2=''
               Else
                  File1=''
               End If
            End If
         End If
         Call Merge_Constraints(File1,File2,'UDC',nLambda,iRow_c)
         Call Fix_UDC(iRow_c,nLambda,nsAtom,nStab,.False.)
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If ( nLambda.eq.0) Then
         M=3*nsAtom
         N=nQQ
         NRHS=1
         Call mma_Allocate(Tmp,M,Label='Tmp:M')
*
*        Iterate to avoid too large displacement in cartesians.
*        If the maximum displacement is more than 2*Beta, reduce the step.
*
*        Initial setup to ensure fCart=1.0 at first iteration.
         rInter=Beta
         fCart=Ten
         rCart=fCart*rInter
         nLoop=0
         Do While (rCart.ge.Two*Beta)
            nLoop=nLoop+1
            If (nLoop.gt.10) Exit
            If (rCart.gt.rInter) Then
              fCart=fCart*rInter/rCart
            Else
              fCart=fCart*0.9D0
            End If
            nA = (Max(nQQ,kIter)+1)**2
*                                                                      *
************************************************************************
*                                                                      *
            gBeta=One
            xBeta=One
            gg_last=Beta
            dxdx_last=Beta
            Sf=Sqrt(Two)
            If (iOpt_RS.eq.0) Then
               kStart=Max(1,kIter-4)
               iEnd=kIter
            Else
               kStart=Max(1,Iter-4)
               iEnd=Iter
            End If
            Thr=1.0D-6
            Do iIter = kStart, iEnd
*
               If (iIter.ne.kIter) Then
                  dxdx=
     &             Sqrt(DDot_(nQQ,Shift(1,iIter),1,Shift(1,iIter),1))
*
                  If (dxdx.lt.0.75D0*dxdx_last.and.
     &                dxdx.lt.(Beta-Thr)) Then
*                    Increase trust radius
C                    xBeta=xBeta*Sf
                     xBeta=Min(Two,xBeta*Sf)
                  Else If (dxdx.gt.1.25D0*dxdx_last.or.
     &                     dxdx.ge.(Beta+Thr)) Then
*                    Reduce trust radius
                     xBeta=Max(One/Five,xBeta/Sf)
                  End If
                  dxdx_last=dxdx
               End If
*
               gg=Sqrt(DDot_(nQQ-nLambda,dqInt(:,iIter),1,
     &                                      dqInt(:,iIter),1))
               If (gg.lt.0.75D0*gg_last.and.gg.lt.(Beta-Thr)) Then
*                 Increase trust radius
C                 gBeta=gBeta*Sf
                  gBeta=Min(Two,gBeta*Sf)
               Else If (gg.gt.1.25D0*gg_last.or.gg.ge.(Beta+Thr)) Then
*                 Reduce trust radius
                  gBeta=Max(One/Five,gBeta/Sf)
               End If
               gg_last=gg
*
            End Do
            tBeta= Max(Beta*Min(xBeta,gBeta),Beta/Ten)
C           Write (6,*) 'tBeta=',tBeta
*                                                                      *
************************************************************************
*                                                                      *
*---------- Compute updated geometry in Internal coordinates
*
            fact=One
            qBeta=fCart*tBeta
            Thr_RS=1.0D-7
            Do
               Step_Trunc_=Step_Trunc
               Call Newq(qInt,nQQ,kIter,Shift,Hessian,dqInt,
     &                   ErrVec,EMtrx,RHS,
     &                   AMat,nA,
     &                   qBeta,nFix,Index,
     &                   Energy,Step_Trunc_,Thr_RS)
               If (Step_Trunc.eq.'N') Step_Trunc=' '
               If (iOpt_RS.eq.0) Then
                  If (Step_Trunc_.eq.'N') Step_Trunc_=' '
                  Step_Trunc=Step_Trunc_
                  Exit
               End If
               If (Step_Trunc//Step_Trunc_.eq.' *') Step_Trunc='.'
*
               qInt(:,kIter+1)=qInt(:,kIter)+Shift(:,kIter)
               Call Dispersion_Kriging_Layer(qInt(1,kIter+1),Disp,
     &                                       nQQ)
#ifdef _DEBUGPRINT_
               Write (6,*) 'Disp,Beta_Disp=',Disp,Beta_Disp
#endif
               fact=Half*fact
               qBeta=Half*qBeta
               If (One-disp/Beta_Disp.gt.1.0D-3) Exit
               If ((fact.lt.1.0D-5) .or. (disp.lt.Beta_Disp)) Exit
               Step_Trunc='*'
            End Do
#ifdef _DEBUGPRINT_
               Write (6,*) 'Step_Trunc=',Step_Trunc
#endif
*
            Call MxLbls(nQQ,dqInt(1,kIter),Shift(1,kIter),Lbl)
*
*---------- Set shift vector to zero for frozen internal coordinates.
*
            If (nFix.gt.0)
     &         Call dcopy_(nFix,[Zero],0,Shift(iInt+1,kIter),1)
*
*           Rough conversion to Cartesians
*
            Call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     Shift(1,kIter),Tmp)
            rInter=Sqrt(dDot_(N,Shift(1,kIter),1,Shift(1,kIter),1))
            rCart=Zero
            Do i=1,nsAtom
               rCart=Max(rCart,
     &                   Sqrt(dDot_(3,Tmp(i),nsAtom,Tmp(i),nsAtom)))
            End Do
         End Do
         Call mma_deallocate(Tmp)
         Call Put_dScalar('Max error',Zero)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*        Optimization with constraints using Lagrangian technique.
*
*        Allocate memory for the new arrays
*
         Call mma_allocate(R,nLambda,kIter,Label='R')
         Call mma_allocate(dRdq,nQQ,nLambda,kIter,Label='dRdq')
         dRdq(:,:,:)=Zero
         Call mma_allocate(T,nQQ,nQQ,Label='T')
         T(:,:)=Zero
         Call mma_allocate(dy,nLambda,Label='dy')
         Call mma_allocate(dx,(nQQ-nLambda),kIter,Label='dx')
         dx(:,:)=Zero
         Call mma_allocate(dEdq,nQQ,kIter,Label='dEdQ')
         Call mma_allocate(du,nQQ,Label='du')
         Call mma_allocate(X,(nQQ-nLambda),(kIter+1),Label='X')
         X(:,:)=Zero
         call mma_allocate(dEdx,nQQ-nLambda,kIter,Label='dEdx')
         Call mma_allocate(Wess,nQQ-nLambda,nQQ-nLambda,
     &                     Label='Wess')
         Wess(:,:)=Zero
         Call mma_allocate(Energy_L,kIter,Label='Energy_L')
*
         Energy_L(:)=Energy(1:kIter)
*
         If (nQQ+nLambda.gt.SIZE(Lbl)) Then
            Call mma_allocate(Lbl_tmp,nQQ+nLambda,Label='Lbl_tmp')
            Lbl_tmp(:)=''
            Lbl_tmp(1:SIZE(Lbl))=Lbl(:)
            Call mma_deallocate(Lbl)
            Call mma_allocate(Lbl,nQQ+nLambda,Label='Lbl')
            Lbl(:)=Lbl_tmp(:)
            call mma_deallocate(Lbl_tmp)
         End If
*
         nBVec=iRow_c-nLambda-1
*
         n1=3*nsAtom
         Call mma_allocate(BM,n1,nLambda,Label='BM')
         Call mma_allocate(dBM,n1,n1,nLambda,Label='dBM')
*
         Call mma_allocate(BVec,3*nsAtom,nBVec,Label='BVec')
         Call mma_allocate(Value,nBVec,Label='Value')
         Call mma_allocate(Value0,nBVec,Label='Value0')
         Call mma_allocate(CInt,nLambda,Label='CInt')
         Call mma_allocate(CInt0,nLambda,Label='CInt0')
         Call mma_allocate(Mult,nBVec**2,Label='Mult')
         Call mma_allocate(iFlip,nBVec,Label='iFlip')
         Call mma_allocate(dBVec,nBVec*(3*nsAtom)**2,Label='dBVec')
*                                                                      *
************************************************************************
*                                                                      *
*        Compute the constraints
*                                                                      *
************************************************************************
*                                                                      *
         Do lIter = 1, kIter
*                                                                      *
************************************************************************
*                                                                      *
            Call DefInt2(BVec,dBVec,nBVec,BM,nLambda,nsAtom,
     &                   iRow_c,Value,cInt,cInt0,Lbl(nQQ+1),
     &                   (lIter.eq.kIter).and.First_MicroIteration,
     &                   Mult,dBM,Value0,lIter,iFlip)
*
*           Assemble r
*
            R(:,lIter)=cInt(:)-cInt0(:)
*
*           Assemble dC/dQ: Solve  B dC/dQ = dC/dx
*
*           where B = dQ/dx
*
            dRdq(:,:,lIter)=Zero
#ifdef _DEBUGPRINT_
            Write (Lu,*) 'Update_inner: lIter=',lIter
            Call RecPrt('Update_inner: dQ/dx(BMx)',' ',BMx,n1,nQQ)
            Call RecPrt('Update_inner: dC/dx(BM)',' ',BM,n1,nLambda)
#endif
*
            M=n1
            N=nQQ
            NRHS=nLambda
*
*           Temporary fix of the dC/dx vector which always
*           is propted up with the full degeneracy factor.
*
            If (.NOT.Curvilinear) Then
               Do iLambda=1,nLambda
                  Do i = 1, n1
                     iAtom = (i+2)/3
                     ixyz = i - (iAtom-1)*3
                     BM(i,iLambda)=BM(i,iLambda)/Degen(ixyz,iAtom)
                  End Do
               End Do
            End If
            If (lIter.eq.kIter) Then
               LudRdX=30
               Call DaName(LudRdX,'dRdX')
               iAd=0
               iDum(1)=nLambda
               Call iDaFile(LudRdX,1,iDum,1,iAd)
               iDum(1)=n1
               Call iDaFile(LudRdX,1,iDum,1,iAd)
               Call dDaFile(LudRdX,1,BM,nLambda*n1,iAd)
               Call DaClos(LudRdX)
            End If
            Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     BM,dRdq(1,1,lIter))
#ifdef _DEBUGPRINT_
            Call RecPrt('Update_inner: dRdq(1,1,lIter)',' ',
     &                   dRdq(1,1,lIter),nQQ,nLambda)
#endif
*                                                                      *
************************************************************************
*                                                                      *
         End Do     ! lIter
*                                                                      *
************************************************************************
*                                                                      *
         Call mma_deallocate(dBVec)
         Call mma_deallocate(iFlip)
         Call mma_deallocate(Mult)
         Call mma_deallocate(cInt0)
         Call mma_deallocate(cInt)
         Call mma_deallocate(Value0)
         Call mma_deallocate(Value)
         Call mma_deallocate(BVec)
*
#ifdef _DEBUGPRINT_
         Call RecPrt('Update_inner: R',' ',R,nLambda,kIter)
         Call RecPrt('Update_inner: dRdq',' ',dRdq,nQQ*nLambda,
     &               kIter)
         Do i = 1, nLambda
            Write (6,*) ' iLambda=',i
            Call RecPrt('Update_inner: d2C/dx2',' ',dBM(1,1,i),n1,n1)
         End Do
         If (Curvilinear) Call dBPrint(nQQ,nDimBC)
#endif
         Call mma_deallocate(BM)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Assemble d^2C/dQ^2
*
*        The expression is written as
*
*        dQ/dx * d^2C/dQ^2 * (dQ/dx)^T = d^2C/dx^2 - Sum (d^2Q/dx^2 * dC/dQ)
*
*        Dimensions
*        dQ/dx     :  3*nsAtom x nQQ
*        d^2C/dx^2 :  3*nsatom x 3*nsatom
*        d^2C/dQ^2 :  nQQ   x nQQ
*
#ifdef _NEW_CODE_
*        Compute Sum (d^2Q/dx^2 * dC/dQ)
*
         Call mma_allocate(QC,nDimBC,nDimBC,nLambda,Label='QC')
         QC(:,:,:)=Zero
         If (Curvilinear) Call dBMult(dRdq(1,1,kIter),
     &                                QC,nQQ,nDimBC,nLambda)
#ifdef _DEBUGPRINT_
         Write (Lu,*) 'Update_inner: kIter=',kIter
         Call RecPrt('dRdq(1,1,kIter)',' ',dRdq(1,1,kIter),nQQ,1)
         Do iLambda=1,nLambda
            Write (6,*) 'Update_inner: iLambda=',iLambda
            Call RecPrt('Update_inner: QC',' ',QC(1,1,iLambda),
     &                   nDimBC,nDimBC)
         End Do
#endif
*
*        Subtract the term from d^2C/dx^2
*
         Do k = 1, nLambda
            j = 0
            Do jx = 1, n1
               jAtom = (jx+2)/3
               jxyz = jx - (jAtom-1)*3
               If (Smmtrc(jxyz,jAtom)) Then
                  j = j + 1
*
                  i=0
                  Do ix = 1, n1
                     iAtom = (ix+2)/3
                     ixyz = ix - (iAtom-1)*3
                     If (Smmtrc(ixyz,iAtom)) Then
                        i = i + 1
                        dBM(ix,jx,k) = dBM(ix,jx,k) - QC(i,j,k)
                     End If
                  End Do
               End If
            End Do
         End Do
         Call mma_deallocate(QC)
#endif
*
         Call mma_allocate(Scr2,nQQ*n1*nLambda,Label='Scr2')
         Call mma_allocate(Scr1,nQQ*n1*nLambda,Label='Scr1')
#ifdef _DEBUGPRINT_
         Call RecPrt('Update_inner: d^2C/dx^2(dBM)',' ',dBM,n1**2,
     &               nLambda)
#endif
*
*        Temporary fix of the d^2C/dx^2 vector which always
*        is propted up with the full degeneracy factor.
*
         If (.NOT.Curvilinear) Then
            Do k = 1, nLambda
               Do j = 1, n1
                  jAtom = (j+2)/3
                  jxyz = j - (jAtom-1)*3
                  Do i = 1, n1
                     iAtom = (i+2)/3
                     ixyz = i - (iAtom-1)*3
                     dBM(i,j,k)=dBM(i,j,k)/(Degen(ixyz,iAtom)
     &                                     *Degen(jxyz,jAtom))
                  End Do
               End Do
            End Do
         End If
#ifdef _DEBUGPRINT_
         Call RecPrt('Update_inner: d^2C/dx^2(dBM)',' ',dBM,n1**2,
     &               nLambda)
#endif
*
*        Solve first for y in
*                dQ/dx y = d^2C/dx^2 - Sum (d^2Q/dx^2 * dC/dQ)
*
*        where y   = d^2C/dQ^2 * (dQ/dx)^T
*
*        and   y^T = dQ/dx * d^2C/dQ^2
*
         M=3*nsAtom
         N=nQQ
         NRHS=n1*nLambda
         Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,dBM,Scr1)
*
*        Generate y^T in Scr2
*
         Call TRNSPS(nQQ,n1*nLambda,Scr1,Scr2)
#ifdef _DEBUGPRINT_
         Call RecPrt('d^2C/dQ^2 * (dQ/dx)^T',' ',Scr1,nQQ,n1*nLambda)
         Call RecPrt('dQ/dx * d^2C/dQ^2',' ',Scr2,n1*nLambda,nQQ)
#endif
         Call mma_deallocate(dBM)
*
*        Followed by solve for x in B x = y^T for which B = dQ/dx
*
         NRHS=nLambda*nQQ
         Call Eq_Solver('N',M,N,NRHS,BMx,Curvilinear,Degen,Scr2,Scr1)
         Call mma_allocate(d2L,nQQ,nQQ,nLambda,Label='d2L')
*
         Call TRNSPS(nQQ*nLambda,nQQ,Scr1,d2L)
#ifdef _DEBUGPRINT_
         Call RecPrt('Scr1',' ',Scr1,nQQ*nLambda,nQQ)
         Do i = 1, nLambda
            Write (6,*) ' iLambda=',i
            Call RecPrt('Update_inner: d2L',' ',d2L(:,:,i),
     &                  nQQ,nQQ)
         End Do
#endif
         Call mma_deallocate(Scr2)
         Call mma_deallocate(Scr1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*------- Compute updated geometry in Internal coordinates
*
         Mode_Save=Mode
         Mode=0
         M=3*nsAtom
         N=nQQ
         NRHS=1
         Call mma_Allocate(Tmp,M,Label='Tmp:M')
*
*        Iterate to avoid too large displacement in cartesians.
*        If the maximum displacement is more than 2*Beta, reduce the step.
*
*        Initial setup to ensure fCart=1.0 at first iteration.
         Thr_RS=1.0D-7
         rInter=Beta
         fCart=Ten
         rCart=fCart*rInter
         nLoop=0
         Do While (rCart.ge.Two*Beta)
            nLoop=nLoop+1
            If (nLoop.gt.100) Exit
            If (rCart.gt.rInter) Then
              fCart=fCart*rInter/rCart
            Else
              fCart=fCart*0.9D0
            End If
            qBeta=fCart*Beta
            Call Con_Opt(R,dRdq,T,dqInt,Lambda,qInt,Shift,dy,dx,
     &                dEdq,du,x,dEdx,Wess,GNrm(kIter),
     &                nWndw,Hessian,nQQ,kIter,
     &                iOptH_,jPrint,Energy_L,nLambda,
     &                ErrVec,EMtrx,RHS,
     &                AMat,nA,qBeta,qBeta_Disp,nFix,
     &                Index,Step_Trunc,Lbl,
     &                d2L,nsAtom,
     &                iOpt_RS,Thr_RS,iter,
     &                First_Microiteration)
            If (iOpt_RS.eq.1) Exit
*
*           Rough conversion to Cartesians
*
            Call Eq_Solver('T',M,N,NRHS,BMx,Curvilinear,Degen,
     &                     Shift(1,kIter),Tmp)
            rInter=Sqrt(dDot_(N,Shift(1,kIter),1,Shift(1,kIter),1))
            rCart=Zero
            Do i=1,nsAtom
               rCart=Max(rCart,
     &                   Sqrt(dDot_(3,Tmp(i),nsAtom,Tmp(i),nsAtom)))
            End Do
         End Do
         Mode=Mode_Save
         Call mma_deallocate(Tmp)
*
#ifdef _DEBUGPRINT_
         Write (Lu,*)
         Write (Lu,*) '********************************************'
         Write (Lu,*) '* Lagrange multipliers for the constraints *'
         Write (Lu,*) '********************************************'
         Write (Lu,'(1X,A,2X,ES13.6)')
     &       (Lbl(nQQ+iInt),-One*Lambda(iInt,mIter),
     &        iInt=1,nLambda)
         Write (Lu,*)
#endif
*
         Call mma_deallocate(d2L)
         Call mma_deallocate(Energy_L)
         Call mma_deallocate(Wess)
         Call mma_deallocate(dEdx)
         Call mma_deallocate(x)
         Call mma_deallocate(du)
         Call mma_deallocate(dEdq)
         Call mma_deallocate(dx)
         Call mma_deallocate(dy)
         Call mma_deallocate(T)
         Call mma_deallocate(dRdq)
         Call mma_deallocate(R)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_Deallocate(Hessian)
      Call mma_Deallocate(RHS)
      Call mma_Deallocate(EMtrx)
      Call mma_Deallocate(ErrVec)
      Call mma_Deallocate(Index)
      Call mma_Deallocate(AMat)
#ifdef _WARNING_WORKAROUND_
      If (Smmtrc(1,1)) i=nDimBC
#endif
*
      Return
      End
