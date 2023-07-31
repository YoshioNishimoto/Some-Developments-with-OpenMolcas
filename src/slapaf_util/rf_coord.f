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
*#define _DEBUGPRINT_
      Subroutine RF_Coord(
     &                 nq,nsAtom,iIter,nIter,Cx,
     &                 Process,Value,
     &                 nB,qLbl,iRef,fconst,
     &                 rMult,LuIC,Indq,
     &                 Proc_dB,mB_Tot,mdB_Tot,
     &                 BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,nqB)
      use Symmetry_Info, only: nIrrep, iOper, VarR, VarT
      use Slapaf_Info, only: nStab, iCoSet, dMass
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "print.fh"
      Real*8 Cx(3,nsAtom,nIter), fconst(nB), Value(nB,nIter), rMult(nB),
     &       Trans(3), RotVec(3), RotMat(3,3),
     &       BM(nB_Tot), dBM(ndB_Tot)
      Integer   nqB(nB),
     &          Indq(3,nB), iBM(nB_Tot), idBM(2,ndB_Tot)
      Logical Process, PSPrint, Proc_dB, Invariant
      Character*3 TR_type(6)
      Character*14 Label, qLbl(nB)
#include "ddvdt_RF.fh"
      Real*8, Dimension(:), Allocatable :: xMass
      Real*8, Dimension(:,:), Allocatable :: currXYZ, Ref123, Grad,
     &                                       dRVdxyz, Hess
      Real*8, Dimension(:,:,:), Allocatable :: d2RV
      Integer, Dimension(:), Allocatable :: Ind, iDCR
      Data TR_type/'Tx ','Ty ','Tz ','Ryz','Rzx','Rxy'/
*
      iRout=151
      iPrint=nPrint(iRout)
#ifdef _DEBUGPRINT_
      iPrint=99
#endif
*
      If (.Not.VarR.and..Not.VarT) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
      nqRF=0
      PSPrint=.False.
      If (iPrint.ge.99) PSPrint=.True.
      If (PSPrint) Write (6,*) ' Enter RF_Coords.'
*
*---- Find nCent and allocate
*
      nCent=0
      Do iAtom = 1, nsAtom
         nCent=nCent+nIrrep/nStab(iAtom)
      End Do
      mB = nCent*3
      Call mma_allocate(currXYZ,3,nCent,label='currXYZ')
      Call mma_allocate(Ref123,3,nCent,label='Ref123')
      Call mma_allocate(Grad,3,nCent,label='Grad')
      Call mma_allocate(dRVdxyz,3,3*nCent,label='dRVdxyz')
      Call mma_allocate(xMass,nCent,label='xMass')
      Call mma_allocate(Ind,nCent,label='Ind')
      Call mma_allocate(iDCR,nCent,label='iDCR')
      Call mma_allocate(Hess,mB,mB,label='Hess')
*
*---- Find index of RF center (origin), etc
*
      iCent=0
      Do iAtom = 1, nsAtom
*
         Do i = 0, nIrrep/nStab(iAtom)-1
            iCent = iCent + 1
            Call OA(iCoSet(i,iAtom),Cx(1:3,iAtom,iIter),
     &              CurrXYZ(1:3,iCent))
            Call OA(iCoSet(i,iAtom),Cx(1:3,iAtom,iRef),
     &              Ref123(1:3,iCent))

*
            Ind(iCent)=iAtom
            iDCR(iCent)=iCoSet(i,iAtom)
         End Do
      End Do
*
C     Fact=One
C     If (.Not.VarR) Fact=2.0D-2
*
*     Write (6,*) 'nCent=',nCent
*     Write (6,*) (Ind(iCent),iCent=1,nCent)
*
      TMass = Zero
      Do iCent = 1, nCent
         iAtom = Ind(iCent)
         xMass(iCent) = dMass(iAtom)
         TMass = TMass + dMass(iAtom)
      End Do
*---- Loop over cartesian components
*
      Do ixyz = 1, 3
*
         Invariant=.False.
         iTest=2**(ixyz-1)
         Do iSym = 0, nIrrep-1
            If (iOper(iSym).eq.iTest) Invariant=.True.
         End Do
         If (Invariant) Go To 199
*
*------- Compute total mass and center of mass of the molecule, the center
*        of the RF cavity is at origin with infinite mass. Hence, the latter
*        is ignored!
*
         COM_xyz = Zero
         Do iCent = 1, nCent
            COM_xyz = COM_xyz + currXYZ(ixyz,iCent)*xMass(iCent)
         End Do
         COM_xyz=COM_xyz/TMass
*
         If (.Not.VarT) Go To 199
*
         iDeg=1
         Deg=Sqrt(DBLE(iDeg))
*
         nq = nq + 1
         If (.Not.Process) mB_Tot = mB_Tot + mB
         If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
         nqRF = nqRF + 1
         Write (LuIC,'(A,I2.2,2A)')
     &              'TR',nqRF,' = ',TR_type(ixyz)
         Label=' '
         Write (Label,'(A,I2.2)') 'TR',nqRF
*
         Val = COM_xyz
*
*------- Compute the gradient
*
         call dcopy_(mB,[Zero],0,Grad,1)
         Do iCent = 1, nCent
            iAtom=Ind(iCent)
*           Write (6,*) 'iAtom,iCOM=',iAtom,iCOM
            Grad(ixyz,iCent) = dMass(iAtom)/TMass
         End Do
#ifdef _DEBUGPRINT_
         Call RecPrt('Grad (Trans)',' ',Grad,3,nCent)
#endif
*
*------- Second derivative is trivially zero!
*
         Call FZero(Hess,mB**2)
         If (Process) Then
*
            Indq(1,nq) = -2**(ixyz)/2
            Indq(2,nq) = 0
            Indq(3,nq) = 0
*
C           fconst(nq)=Sqrt(Fact*Trans_Const)
            fconst(nq)=Sqrt(Trans_Const)
            rMult(nq)=Deg
*
            Value(nq,iIter)=Val
            qLbl(nq)=Label
*
*--------   Project the gradient vector
*
            Call ProjSym(nCent,Ind,currXYZ,
     &                   iDCR,Grad,
     &                   Hess,mB_Tot,mdB_Tot,
     &                   BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                   Proc_dB,nqB,nB,nq,rMult(nq))
*
         End If
*
 199     Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
C     Write (6,*) 'VarR=',VarR
      If (.Not.VarR) Go To 98
*
*     A la Malmqvist
*
      nOrder=2
      nMass=nCent
      Call FZero(Trans,3)
      Call FZero(RotVec,3)
      Call mma_allocate(d2RV,3,3*nCent,3*nCent,label='d2RV')
#ifdef _DEBUGPRINT_
      Call RecPrt('xMass',' ',xMass,1,nMass)
#endif
      Call RotDer(nMass,xMass,currXYZ,ref123,trans,RotAng,
     &            RotVec,RotMat,nOrder,dRVdXYZ,d2RV)
#ifdef _DEBUGPRINT_
      Call RecPrt('RotVec',' ',RotVec,1,3)
      Call RecPrt('RotMat',' ',RotMat,3,3)
      Call RecPrt('dRVdXYZ',' ',dRVdXYZ,3,3*nMass)
#endif
*
      Do ixyz = 1, 3
*
         Invariant=.False.
         If (ixyz.eq.1) Then
            iTest=6
         Else If (ixyz.eq.2) Then
            iTest=5
         Else
            iTest=3
         End If
         Do iSym = 0, nIrrep-1
            If (iOper(iSym).eq.iTest) Invariant=.True.
         End Do
         If (Invariant) Go To 299
*
         jxyz = ixyz+1
         If (jxyz.gt.3) jxyz=1
         kxyz = jxyz+1
         If (kxyz.gt.3) kxyz=1
         iDeg=1
         Deg=Sqrt(DBLE(iDeg))
*
         nq = nq + 1
         If (.Not.Process) mB_Tot = mB_Tot + mB
         If (.Not.Proc_dB) mdB_Tot = mdB_Tot + mB**2
         nqRF = nqRF + 1
         Write (LuIC,'(A,I2.2,2A)')
     &              'TR',nqRF,' = ',TR_type(ixyz+3)
         Label=' '
         Write (Label,'(A,I2.2)') 'TR',nqRF
*
         Val = RotVec(ixyz)
*
*------- Compute the gradient
*
         call dcopy_(mB,[Zero],0,Grad,1)
         call dcopy_(mB,dRVdXYZ(ixyz,1),3,Grad,1)
#ifdef _DEBUGPRINT_
         Call RecPrt('Grad (Rot)',' ',Grad,3,nCent)
#endif
*
*------- Second derivative
*
         Call FZero(Hess,mB**2)
         If (Proc_dB) Call DCopy_(mB**2,d2RV(ixyz,1,1),3,Hess,1)
*
         If (Process) Then
*
            Indq(1,nq) = -( 2**(jxyz)/2  + 2**(kxyz)/2)
            Indq(2,nq) = 0
            Indq(3,nq) = 0
*
            fconst(nq)=Sqrt(Rot_Const)
            rMult(nq)=Deg
*
            Value(nq,iIter)=Val
            qLbl(nq)=Label
*
*--------   Project the gradient vector
*
            Call ProjSym(nCent,Ind,currXYZ,
     &                   iDCR,Grad,
     &                   Hess,mB_Tot,mdB_Tot,
     &                   BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,
     &                   Proc_dB,nqB,nB,nq,rMult(nq))
*
         End If
*
 299     Continue
      End Do
      Call mma_deallocate(d2RV)
C     Write (6,*) 'nqRF=',nqRF
*                                                                      *
************************************************************************
*                                                                      *
 98   Continue
      Call mma_deallocate(currXYZ)
      Call mma_deallocate(Ref123)
      Call mma_deallocate(Grad)
      Call mma_deallocate(dRVdxyz)
      Call mma_deallocate(xMass)
      Call mma_deallocate(Ind)
      Call mma_deallocate(iDCR)
      Call mma_deallocate(Hess)
 99   Continue
      Return
      End
