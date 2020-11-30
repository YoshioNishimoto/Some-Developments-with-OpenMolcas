!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019, Marjan Khamesian                                 *
!               2019, Roland Lindh                                     *
!***********************************************************************
      SubRoutine TWLInt(Alpha,nAlpha,Beta,nBeta,Zeta,ZInv,rKappa,P,    &
                       Final,nZeta,nIC,nComp,la,lb,A,RB,nHer,          &
                       Array,nArr,kVector,nOrdOp,lOper,iChO,           &
                       iStabM,nStabM,                                  &
                       PtChrg,nGrid,iAddPot)
!***********************************************************************
!                                                                      *
! Object: kernel routine for the computation of integrals for the      *
!         interaction between matter and light with orbital angular    *
!         momentum.                                                    *
!***********************************************************************
      use Her_RW
      use Real_Spherical
      Implicit None
      External NrOpr
#include "stdalloc.fh"
!
!     External Arrays and integers
!
      Integer nZeta, la, lb, nIC, nAlpha, nBeta, nArr, nComp, nOrdOp,  &
              nStabM, lc, ld
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),         &
             Zeta(nZeta), Alpha(nAlpha), Beta(nBeta), Gamma(nGamma)    &
             ZInv(nZeta), rKappa(nZeta),                               &
             P(nZeta,3), A(3), RB(3),                                  &
             Array(nZeta*nArr), kvector(3)
      Integer iStabM(0:nStabM-1), iDCRT(0:7),                          &
              iStabO(0:7), lOper(nComp), iChO(nComp)
      Integer nHer, nGrid, iAddPot, NrOpr
      Real*8  PtChrg
      Logical ABeq(3)
!
!     Local Arrays and integers
!
      Integer iOper(0:7)
      Integer iAlpha, iBeta, ixyz
      Integer ipA, ipAOff, ipAxyz, ipB, ipBOff, ipBxyz, ipQxyz, ipRes, &
              ipVxyz, nip, icomp, lDCRT, llOper, LmbdT, nDCRT,  &
              nIrrep, nOp, nStabO, ipP, lAng, ipScr, iOff
      Real*8 Zero, Half, One, Two, Three, Four, Rxy, Fi1, Rxyz, Fi2, Fi3
      Real*8, Allocatable:: TransM(:,:)
      Real*8 kVector_local(3)
      Real*8 A_Local(3), RB_Local(3)

      Integer :: i, l
      Real*8 :: dmm
!
!     since we do not include info.fh -- it is a f77 style line contiuation file --
!     we have to have a dirty work around here.
!
      Call Peek_iScalar('nSym',nIrrep)
      iOper(:)=0
      Call Peek_iOper(iOper,nIrrep)
#define _DEBUG_
!
      Zero =0.0D0
      Half =0.5D0
      One  =1.0D0
      Two  =2.0D0
      Three=3.0D0
      Four =4.0D0
      lAng= 1             ! Temporary set value
      lAng= 0             ! Temporary set value
      nip = 1
!
      If (lAng.ne.0 .and. nOrdOp.eq.0) Then
         Write (6,*) "lAng.ne.0 .and. nOrdOp.eq.0"
         Call Abend()
      End If
!
!     We need here a transformation of the Cartesian coordinates in which
!     the k-vector coincides with the z-vector direction.
!     In particular, we will transform P to the new coordinate system.
!
#ifdef _TEMP_SKIP_
      If (lAng.eq.0) Then
         kVector_Local(1)=kVector(1)
         kVector_Local(2)=kVector(2)
         kVector_Local(3)=kVector(3)
         Go To 114
      End If
#endif
!
!     Use Euler angles (z-x-z). We skip the last rotation.
!     Rotate around the z-axis so that the x-component becomes zero.
      Rxy=Sqrt(kvector(1)**2 + kvector(2)**2)
      Fi1= ATAN2(kvector(1),kvector(2))
      kVector_Local(1)=Rxy
      kVector_Local(2)=Zero
      kVector_Local(3)=kVector(3)
      Write(6,*) 'kVector=',kVector
      Write(6,*) 'kVector_local=',kVector_Local
!
!     Rotate around the x-axis so that the y-component becomes zero.
      Rxyz=Sqrt(kVector_Local(2)**2 + kVector_Local(3)**2)
      Fi2 = ATAN2(kVector_Local(2),kVector_Local(3))
      kVector_Local(1)=Zero
      kVector_Local(2)=Zero
      kVector_Local(3)=Rxyz
      Write(6,*) 'kVector_local=',kVector_Local
!
      Fi3=0.0D0
!
      Call mma_Allocate(TransM,3,3,Label="TransM")
      TransM(1,1)= Cos(Fi1)
      TransM(2,1)= Sin(Fi1)
      TransM(3,1)= Zero
      TransM(1,2)=-Sin(Fi1)*Cos(Fi2)
      TransM(2,2)= Cos(Fi1)*Cos(Fi2)
      TransM(3,2)=          Sin(Fi2)
      TransM(1,3)= Sin(Fi1)*Sin(Fi2)
      TransM(2,3)=-Cos(Fi1)*Sin(Fi2)
      TransM(3,3)=          Cos(Fi2)
#ifdef _DEBUG_
      Write (6,*) 'Fi1,Fi2,Fi3=',Fi1,Fi2,Fi3
      Call RecPrt('TransM',' ',TransM,3,3)
      Call RecPrt('A',' ',A,3,1)
      Call RecPrt('RB',' ',RB,3,1)
      Call RecPrt('P',' ',P,nZeta,3)
#endif
!
!     Transform A, RB, and P
!
      ipP=nip
      nip=nip+nZeta*3
!
      Call DGEMM_('N','N',3,1,3,                                        &
                   1.0D0,TransM,3,                                      &
                         A,3,                                           &
                   0.0D0,A_Local,3)
      Call DGEMM_('N','N',3,1,3,                                        &
                   1.0D0,TransM,3,                                      &
                         RB,3,                                          &
                   0.0D0,RB_Local,3)
      Call DGEMM_('N','T',nZeta,3,3,                                    &
                   1.0D0,P,nZeta,                                       &
                         TransM,3,                                      &
                   0.0D0,Array(ipP),nZeta)
      Call mma_deallocate(TransM)
#ifdef _DEBUG_
      Call RecPrt('A_local',' ',A_local,3,1)
      Call RecPrt('RB_local',' ',RB_local,3,1)
      Call RecPrt('P_local',' ',Array(ipP),nZeta,3)
#endif
!
 114  Continue
#ifdef _TEMP_SKIP_
      If (lAng.eq.0) Then
         Call TWLInt_Internal(Array,A,RB,P)
      Else
         Call TWLInt_Internal(Array,A_Local,RB_Local,Array(ipP))
      End If
#else
      Call TWLInt_Internal(Array,A_Local,RB_Local,Array(ipP))
#endif
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!====================================================================
!
#ifdef _TEMP_SKIP_
      If (lAng.eq.0) Go To 263
#endif
!     Now when all Cartesian components have been computed we
!     transform back to the coordinate system of the molecule.
!
!====================================================================
!     Backtransform the integrals to the original coordinate system
!
      ipScr=nip
      nip = nip - nZeta*((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)
!
!     1) transform the integrals to spherical harmonics (with the contaminant).
!        ij,a,b  ->  B,ij,a
!
      Call DGEMM_('T','T',                                                 &
                  (lb+1)*(lb+2)/2,nZeta*((la+1)*(la+2)/2),(lb+1)*(lb+2)/2, &
                  1.0D0,RSph(ipSph(lb)),(lb+1)*(lb+2)/2,                   &
                        Array(ipRes),nZeta*((la+1)*(la+2)/2),              &
                  0.0D0,Array(ipScr),(lb+1)*(lb+2)/2 )
!
!        B,ij,a  ->  A,B,ij
!
      Call DGEMM_('T','T',                                                 &
                  (la+1)*(la+2)/2,((lb+1)*(lb+2)/2)*nZeta,(la+1)*(la+2)/2, &
                  1.0D0,RSph(ipSph(la)),(la+1)*(la+2)/2,                   &
                        Array(ipScr),((lb+1)*(lb+2)/2)*nZeta,              &
                  0.0D0,Array(ipRes),(la+1)*(la+2)/2 )
!
!     2) backtransform the spherical harmonics to the original coordinate system
!
!        A(new),B,ij -> B,ij,A(orig)  (one set of sphericals at the time)
!
!
!     First generate the matrix which transforms between spherical harmonics
!     in different coordinate system. Here the two coordinates systems are
!     related throught theee Euler angles, rotating z-x-z. We do the rotations
!     in the opposite order as we rotated above.
!
      Call mma_Allocate(TransM,(la+1)*(la+2)/2,(la+1)*(la+2)/2,Label="TransM")
      TransM(:,:)=0.0D0
!     Note the order, this is the same order as the transformation matrix
!     Cartesian to Spherical has!
      iOff=1
      Do i = la,0,-2
!==============================================================================
!        Generate the blocks of the transformation matrix and put them into
!        TransM
         Call dmm_transform(i,-Fi3,-Fi2,-Fi1,TransM(iOff,iOff),            &
                            (la+1)*(la+2)/2)
         iOff = iOff + 2*i + 1
!
!==============================================================================
      End Do
!
!     Time to do the actual transformations.
!
!     We like to do Trans(A(old),A(new)) I(A(new),(B,ij)) to I(A(old),(B,ij))
!
!     However, prepare for the next transformation we like to the order of
!     result to be I((B,ij),A(old)). Thus we will have to change the order
!     of the matrices and transpose them.
!
      Call DGEMM_('T','T',                                                 &
                 ((lb+1)*(lb+2)/2)*nZeta,(la+1)*(la+2)/2,(la+1)*(la+2)/2,  &
                 1.0D0,Array(ipRes),((la+1)*(la+2)/2),                     &
                       TransM,(la+1)*(la+2)/2,                             &
                 0.0D0,Array(ipScr),((lb+1)*(lb+2)/2)*nZeta)
      Call mma_deallocate(TransM)
!
!        Repeat the transformation on the second index, B(new)
!
!        B(new),(ij,A(orig)) -> (ij,A(orig)),B(orig)
!
      Call mma_Allocate(TransM,(lb+1)*(lb+2)/2,(lb+1)*(lb+2)/2,Label="TransM")
      TransM(:,:)=0.0D0
      iOff=1
      Do i = lb,0,-2
!==============================================================================
!        Generate the blocks of the transformation matrix and put them into
!        TransM
         Call dmm_transform(i,-Fi3,-Fi2,-Fi1,TransM(iOff,iOff),            &
                            (lb+1)*(lb+2)/2)
         iOff = iOff + 2*i +1
!
!==============================================================================
      End Do
      Call DGEMM_('T','T',                                                 &
                 nZeta*((la+1)*(la+2)/2),(lb+1)*(lb+2)/2,(lb+1)*(lb+2)/2,  &
                  1.0D0,Array(ipScr),nZeta*((la+1)*(la+2)/2),              &
                        TransM,((lb+1)*(lb+2)/2),                          &
                 0.0D0,Array(ipRes),nZeta*((la+1)*(la+2)/2))

      Call mma_deallocate(TransM)
!
!     3) transform the spherical harmonics to Cartesians.
!
!       Find inverse for RSph(ipSph(ld))
!
        Call mma_Allocate(TransM,(lb+1)*(lb+2)/2,(lb+1)*(lb+2)/2,Label='TransM')
        Call DCopy_(((lb+1)*(lb+2)/2)**2,RSph(ipSph(lb)),1,TransM,1)
        Call MatInvert(TransM,(lb+1)*(lb+2)/2)
!====================================================================
!       (ij,A),B -> b,(ij,A)
        Call DGEMM_('N','T',                                                &
                   (lb+1)*(lb+2)/2,nZeta*((la+1)*(la+2)/2),(lb+1)*(lb+2)/2, &
                   1.0D0,TransM,((lb+1)*(lb+2)/2),                          &
                         Array(ipRes),nZeta*((la+1)*(la+2)/2),              &
                   0.0D0,Array(ipScr),(lb+1)*(lb+2)/2)
        Call mma_deallocate(TransM)

!       Find inverse for RSph(ipSph(lc))
!
        Call mma_Allocate(TransM,(la+1)*(la+2)/2,(la+1)*(la+2)/2,Label='TransM')
        Call DCopy_(((la+1)*(la+2)/2)**2,RSph(ipSph(la)),1,TransM,1)
        Call MatInvert(TransM,(la+1)*(la+2)/2)
!....
!       (b,ij),A -> a,(b,ij)
        Call DGEMM_('T','T',                                                &
                   (la+1)*(la+2)/2,((lb+1)*(lb+2)/2)*nZeta,(la+1)*(la+2)/2, &
                   1.0D0,TransM,(la+1)*(la+2)/2,                            &
                         Array(ipScr),((lb+1)*(lb+2)/2)*nZeta,              &
                   0.0D0,Array(ipRes),(la+1)*(la+2)/2 )
        Call mma_deallocate(TransM)
!
!     4) transpose to the correct order
!        (a,b),ij -> ij,(a,b)  make sure that it is in Array(ipRes)
!
        Call DCopy_(((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)*nZeta,Array(ipRes),1, &
                                                              Array(ipScr),1)
        Call Trnsps(((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2),nZeta,Array(ipScr),Array(ipRes))
!
!=====================================================================
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
 263  Continue
!
      llOper=lOper(1)
      Do iComp = 2, nComp
         llOper = iOr(llOper,lOper(iComp))
      End Do
      Call SOS(iStabO,nStabO,llOper)
      Call DCR(LmbdT,iOper,nIrrep,iStabM,nStabM,iStabO,nStabO,          &
               iDCRT,nDCRT)
!
      Do lDCRT = 0, nDCRT-1
!
!--------Accumulate contributions
!
         nOp = NrOpr(iDCRT(lDCRT),iOper,nIrrep)
         Call SymAdO(Array(ipRes),nZeta,la,lb,nComp,Final,nIC,          &
                     nOp,lOper,iChO,One)
!
      End Do
!
      Return
!
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
!**********************************************************************
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!**********************************************************************
! This is to allow type punning without an explicit interface
    Contains
!======================================================================
    SubRoutine TWLInt_Internal(Array,A,RB,P)
      Use Iso_C_Binding
      Real*8, Target :: Array(*)
      Real*8 A(3), RB(3), P(nZeta,3)
      Complex*16, Pointer :: zAxyz(:),zBxyz(:),zQxyz(:),zVxyz(:)
      Real*8 Alpha_, Beta_
!
      ABeq(1) = A(1).eq.RB(1)
      ABeq(2) = A(2).eq.RB(2)
      ABeq(3) = A(3).eq.RB(3)
!
      ipAxyz = nip
      nip = nip + nZeta*3*nHer*(la+1+nOrdOp) * 2
      ipBxyz = nip
      nip = nip + nZeta*3*nHer*(lb+1+nOrdOp) * 2
      ipQxyz = nip
      nip = nip + nZeta*3*(la+1+nOrdOp)*(lb+1+nOrdOp) * 2
      If (nOrdOp.eq.1) Then
         ipVxyz = nip
         nip = nip + nZeta*6*(la+1)*(lb+1) * 2
         ipA = nip
         nip = nip + nZeta
         ipB = nip
         nip = nip + nZeta
         ipRes = nip
         nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      Else
         ipVxyz = nip
         ipA = nip
         ipB = nip
         ipRes = nip
         nip = nip + nZeta*nElem(la)*nElem(lb)*nComp
      End If
      If (nip-1.gt.nArr*nZeta) Then
         Call WarningMessage(2,'TWLINT: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in TWLINT'
         Call Abend()
      End If
!
#ifdef _DEBUG_
      Call RecPrt(' In TWLINT: A',' ',A,1,3)
      Call RecPrt(' In TWLINT: RB',' ',RB,1,3)
      Call RecPrt(' In TWLINT: KVector',' ',kvector_Local,1,3)
      Call RecPrt(' In TWLINT: P',' ',P,nZeta,3)
      Write (6,*) ' In TWLINT: la,lb=',la,lb
#endif
!
!**********************************************************************
!    Computation of intermediate integrals
!
!    Compute all Cartesian components with the OAM free code. Then
!    compute the x and y components in the OAM formalism replacing the
!    OAM-free x and y components.
!
     Final(:,:,:,:)=Zero
!
!    Compute the Cartesian values of the basis functions angular part
!    Note that these arrays are complex.
!
     Call C_F_Pointer(C_Loc(Array(ipAxyz)),zAxyz,[nZeta*3*nHer*(la+nOrdOp+1)])
     Call CCrtCmp(Zeta,P,nZeta,A,zAxyz,la+nOrdOp,HerR(iHerR(nHer)),nHer,ABeq,kvector_Local)
     Call C_F_Pointer(C_Loc(Array(ipBxyz)),zBxyz,[nZeta*3*nHer*(lb+nOrdOp+1)])
     Call CCrtCmp(Zeta,P,nZeta,RB,zBxyz,lb+nOrdOp,HerR(iHerR(nHer)),nHer,ABeq,kvector_Local)
     Nullify(zAxyz,zBxyz)
!
!    Compute the Cartesian components for the multipole moment
!    integrals. The integrals are factorized into components.
!
     Call C_F_Pointer(C_Loc(Array(ipAxyz)),zAxyz,[nZeta*3*nHer*(la+nOrdOp+1)])
     Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,[nZeta*3*(la+nOrdOp+1)*(lb+nOrdOp+1)])
     Call C_F_Pointer(C_Loc(Array(ipBxyz)),zBxyz,[nZeta*3*nHer*(lb+nOrdOp+1)])
     Call CAssmbl(zQxyz,                                             &
                  zAxyz,la+nOrdOp,                                   &
                  zBxyz,lb+nOrdOp,                                   &
                  nZeta,HerW(iHerW(nHer)),nHer)
     Nullify(zAxyz,zBxyz,zQxyz)
!
!    Compute the cartesian components for the velocity integrals.
!    The velocity components are linear combinations of overlap
!    components.
!
     If (nOrdOp.eq.1) Then
        ipAOff = ipA
        Do iBeta = 1, nBeta
           call dcopy_(nAlpha,Alpha,1,Array(ipAOff),1)
           ipAOff = ipAOff + nAlpha
        End Do
!
        ipBOff = ipB
        Do iAlpha = 1, nAlpha
           call dcopy_(nBeta,Beta,1,Array(ipBOff),nAlpha)
           ipBOff = ipBOff + 1
        End Do
!
        Call C_F_Pointer(C_Loc(Array(ipVxyz)),zVxyz,[nZeta*3*(la+1)*(lb+1)])
        Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,[nZeta*3*(la+1)*(lb+1)])
        Call CVelInt(zVxyz,zQxyz,la,lb,Array(ipA),Array(ipB),nZeta)
        Nullify(zVxyz,zQxyz)
     End If
!
!***********************************************************************
#ifdef _TEMP_SKIP_
      If (lAng.ne.0) Then
#endif
!
!     Now compute the OAM x and y component.
!
         Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,[nZeta*3*(la+1)*(lb+1)])
         If (nOrdOp.eq.1) Then
            Call C_F_Pointer(C_Loc(Array(ipVxyz)),zVxyz,[nZeta*3*(la+1)*(lb+1)])
         Else
            zVxyz=zQxyz
         End If
         call OAM_xy(zVxyz,zQxyz,nZeta,la,lb,nOrdOp,Array(ipAOff),Array(ipBOff),       &
                     Array(ipRes),nComp)
         Nullify(zVxyz,zQxyz)
!
!
!    Combine the cartesian components to the full one electron integral.
!
#ifdef _TEMP_SKIP_
      Else ! OAM-free case

         If (nOrdOp.eq.1) Then
            Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,[nZeta*3*(la+1)*(lb+1)])
            Call C_F_Pointer(C_Loc(Array(ipVxyz)),zVxyz,[nZeta*3*(la+1)*(lb+1)*2])
            Call CCmbnVe(zQxyz,nZeta,la,lb,Zeta,rKappa,Array(ipRes),nComp,zVxyz,kvector_Local)
            Nullify(zQxyz,zVxyz)
         Else
            Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,[nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)])
            Call CCmbnMP(zQxyz,nZeta,la,lb,nOrdOp,Zeta,rKappa,Array(ipRes),nComp,kvector_Local)
            Nullify(zQxyz)
         End If
!
!************************************************************************
      End If
#endif
!
!     At this point the integrals are store in Array starting at the position
!     ipRes. The size of the block is nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2
!
    End SubRoutine TWLInt_Internal
!===========================================================================
!
    Integer Function nElem(ixyz)
      Integer ixyz
      nElem = (ixyz+1)*(ixyz+2)/2
    End Function nElem
!
!===========================================================================
!    Subroutine to calculate the xyz-integral
!---------------------------------------------------------------------------    
!    SUPPORTING INFORMATION of "J. Chem. Theory Comput. 2019, 15, 4180−4186"

! -- Performing the rotational averaging, by rotating the molecule
!    by the angles (α, β, γ).    
! -- Using the the rotation matrix in the ZXZ convention,
! -- With the Euler angles α ∈ [0, 2π], β ∈ [0, π], and γ ∈ [0, 2π].
 
    SubRoutine OAM_xyz(Alpha, Beta, Gamma,a,b,c,xp,yp,zp,expo,m,n,k0,w0)

      implicit none
      Real*8 :: Alpha, Beta, Gamma  
      Integer :: a,b,c,m,n                                             
      Real*8  :: xp, yp, zp, expo                         
      Real*8  :: k0, w0  
      Complex*16 :: gau_her  
      Complex*16,parameter :: cone=(1.d0,0.d0), eye=(0.d0,1.d0), &
                              czero=(0.d0,0.d0) !cone=(1.d0,0.d0)=1+0.i  
      Real*8, parameter :: pi=3.141592653589793d0
      Integer :: r,s,t,r1,s1,t1
      Integer :: d,f,g,h,d1,f1,g1,d2,f2         
      Real*8  :: cx,cy,cz,bx,by,bz,cxx,cyy,czz,cxy,cyz,czx
      Complex*16 :: AA, BB, CC, DD, FF
!                                                                                            
      Complex*16 :: ctmp1, ctmp2, ctmp3, ctmp4, ctmp5, ctmp6, &
                    ctmp7, ctmp8, ctmp9, ctmp10, ctmp11
      Complex*16 :: ctmp12, ctmp13, ctmp14, ctmp15, ctmp16,   &
                    ctmp17, ctmp18, ctmp19, ctmp20
!                                                                                            
      Real*8, external :: fact  ! Function to calculate factorial 
!                                                                                            
!- Following are parameters in "S29 - S45":
!      
    bx=2.d0**0.5d0/w0*(-dcos(Beta)*dcos(Gamma)*dsin(Alpha)-dcos(Alpha)*dsin(Gamma))
    by=2.d0**0.5d0/w0*(dcos(Alpha)*dcos(Beta)*dcos(Gamma)-dsin(Alpha)*dsin(Gamma))
    bz=2.d0**0.5d0/w0*(dcos(Gamma)*dsin(Beta))
!
    cx=2.d0**0.5d0/w0*(dcos(Alpha)*dcos(Gamma)-dcos(Beta)*dsin(Alpha)*dsin(Gamma))
    cy=2.d0**0.5d0/w0*(dcos(Gamma)*dsin(Alpha)+dcos(Alpha)*dcos(Beta)*dsin(Gamma))
    cz=2.d0**0.5d0/w0*dsin(Beta)*dsin(Gamma)
! 
    cxx=(dcos(Alpha)**2.d0+dcos(Beta)**2.d0*dsin(Alpha)**2.d0)/w0**2.d0
    cyy=(dcos(Alpha)**2.d0*dcos(Beta)**2.d0+dsin(Alpha)**2.d0)/w0**2.d0
    czz=(dsin(Beta)**2.d0)/w0**2.d0
!
    cxy=2.d0*dsin(Alpha)*dcos(Alpha)*dsin(Beta)**2.d0/w0**2.d0
    cyz=2.d0*dcos(Alpha)*dsin(Beta)*dcos(Beta)/w0**2.d0
    czx=-2.d0*dcos(Beta)*dsin(Alpha)*dsin(Beta)/w0**2.d0
!                                                                                            
    AA=expo*cone+cone*cyy-cone*cxy**2.d0/4.d0/(cone*expo+cone*cxx)                              ! A 
!
    BB=cone*2.d0*expo*yp+eye*dcos(Alpha)*dsin(Beta)*k0 &
         &-cone*(cone*2.d0*expo*xp-eye*dsin(Alpha)*dsin(Beta)*k0)*cxy/2.d0/(cone*expo+cone*cxx) ! B
!
    CC=cone*cxy*czx/2.d0/(cone*expo+cone*cxx)-cone*cyz                                          ! C
!                                                                                            
    DD=cone*expo+cone*czz-cone*czx**2.d0/4.d0/(cone*expo+cone*cxx)-cone*CC**2.d0/4.d0/AA        ! D
!
    FF=cone*2.d0*expo*zp-cone*(cone*2.d0*expo*xp-eye*dsin(Alpha)*dsin(Beta)*k0)*czx/2.d0/(cone*e\ 
    xpo+cone*cxx) &  
         &+cone*BB*CC/2.d0/AA - eye*dcos(Beta)*k0                                               ! F
!                                              
!                                                                                                                                                                                       
    gau_her=czero       
    ctmp1=czero
    ctmp2=czero
    ctmp3=czero
    ctmp4=czero
    ctmp5=czero
!
!   The exponential terms in equation (S28)     
    ctmp6=-cone*expo*xp**2.d0 &
         &+(cone*2.d0*expo*xp-eye*dsin(Alpha)*dsin(Beta)*k0)**2.d0/4.d0/(cone*expo+cone*cxx) &
         &-cone*expo*yp**2.d0 &
         &+BB**2.d0/4.d0/AA &
         &-cone*expo*zp**2.d0 &
         &+FF**2.d0/4.d0/DD
!
    ctmp7=cdexp(ctmp6)
!                                                                                            
    If ((m .lt. 0) .or. (n .lt. 0)) Then
       gau_her=czero
    Else
!       
!      Sum over r, s, t
       Do r=0,m/2           ! 'm' comes from the Hermite polynomials H_m(x) used in (S27)
!                              m = 0, 1, ...
          Do s=0, m-2*r
             Do t=0,m-2*r
                ctmp1=czero                                                                
                If ((s+t) .le. (m-2*r)) Then
!
!                  ctmp1: The first fraction in equation (S28)
                   ctmp1=cone*(-1.d0)**dble(r)*fact(m)*cx**dble(s)*cy**dble(t)*cz**dble(m-2*r-s-t) &
                        &*2.d0**dble(m-2*r)/fact(r)/fact(s)/fact(t)/fact(m-2*r-s-t)
!
!                  Sum over r1(r'), s1(s'), t1(t') 
                   Do r1=0,n/2                       ! 'n' comes from the Hermite polynomials H_n(x) used in (S27)
!                                                       n = 0, 1, ...
                      Do s1=0,n-2*r1
                         Do t1=0,n-2*r1
                            ctmp2=czero
                            If ((s1+t1) .le. (n-2*r1)) Then
!                               
!                              ctmp2: The 2nd fraction from the first term in (S28)
                               ctmp2=cone*(-1.d0)**dble(r1)*fact(n)*bx**dble(s1)*by**dble(t1)*bz**dble(n-2*r1-s1-t1) &
                                    &*2.d0**dble(n-2*r1)/fact(r1)/fact(s1)/fact(t1)/fact(n-2*r1-s1-t1)
!
!                              Sum over d, f, g, h
                               Do d=0,a
                                  Do f=0,(d+s+s1)/2
                                     Do g=0,d+s+s1-2*f
                                        Do h=0,d+s+s1-2*f
                                           ctmp3=czero
                                           If ((g+h) .le. (d+s+s1-2*f)) Then
!
!                                             The 2nd & 3rd terms of (S28)
                                              ctmp3=cone*fact(a)/fact(d)/fact(a-d)*fact(d+s+s1)/fact(f)/fact(g)/fact(h)/fact(d+s+s1-2*f-g-h) &
                                                   &*cdsqrt(pi/(cone*expo+cone*cxx))*(-xp)**dble(a-d) &
                                                   &*2.d0**dble(d+s+s1-2*f)/(4.d0*(cone*expo+cone*cxx))**dble(d+s+s1-f) &
                                                   &*(cone*2.d0*expo*xp-eye*dsin(alpha)*dsin(beta)*k0)**dble(g) &
                                                   &*(-cxy)**dble(h)*(-czx)**dble(d+s+s1-2*f-g-h)
!                                              
!                                             Sum over d1(d'), f1(f'), g1(g')  
                                              Do d1=0,b
                                                 Do f1=0,(t+t1+d1+h)/2
                                                    Do g1=0,t+t1+d1+h-2*f1
!
!                                                      The 4th and 5th terms of (S28) 
                                                       ctmp4=cdsqrt(pi/AA)*fact(b)/fact(d1)/fact(b-d1) &
                                                            &*fact(t+t1+d1+h)/fact(f1)/fact(g1)/fact(t+t1+d1+h-2*f1-g1) &
                                                            &*(-yp)**dble(b-d1)*2.d0**dble(t+t1+d1+h-2*f1)/(4.d0*AA)**dble(t+t1+d1+h-f1) &
                                                            &*BB**dble(g1)*CC**dble(t+t1+d1+h-2*f1-g1)
!
!                                                      Sum over d2(d") and f2(f") 
                                                       Do d2=0,c
                                                          Do f2=0,(m-2*r+n-2*r1+d2+d-2*f-g+d1-2*f1-g1)/2
!
!                                                            The last two terms of (S28)                                         
                                                             ctmp5=cdsqrt(pi/DD)*fact(c)/fact(d2)/fact(c-d2) &
                                                                  &*(-zp)**dble(c-d2) &
                                                                  &*fact(m-2*r+n-2*r1+d2+d-2*f-g+d1-2*f1-g1)/fact(f2)&
                                                                  &/fact(m-2*r+n-2*r1+d2+d-2*f-g+d1-2*f1-g1-2*f2) &
                                                                  &*(2.d0*FF)**dble(m-2*r+n-2*r1+d2+d-2*f-g+d1-2*f1-g1-2*f2) &
                                                                  &/(4.d0*DD)**dble(m-2*r+n-2*r1+d2+d-2*f-g+d1-2*f1-g1-f2)
!
                                                             gau_her=gau_her+ctmp1*ctmp2*ctmp3*ctmp4*ctmp5*ctmp7
!
                                                          End Do ! f2 loop                                            
                                                       End Do ! d2 loop                                              
                                                    End Do ! g1 loop                                                
                                                 End Do ! f1 loop                                                  
                                              End Do ! d1 loop                                                    
                                           End If
                                        End Do ! h loop                                                         
                                     End Do ! g loop                                                           
                                  End Do ! f loop                                                             
                               End Do ! d loop                                                               
                            End If
                         End Do ! t1 loop                                                                  
                      End Do ! s1 loop                                                                    
                   End Do ! r1 loop                                                                      
                End If
             End Do ! t loop                                                                           
          End Do ! s loop                                                                             
       End Do ! r loop                                                                               
       !                                                                                            
    End If 
!                                                                                            
    Return
!*************************************************************************
!*************************************************************************    
  End SubRoutine OAM_xyz
!*************************************************************************
! Function to calculate the factorial
!------------------------------------  
  Function fact(n)
!
    Implicit None
    Integer :: n
    Real*8 :: fact
    Integer :: i
!
    If (n .eq. 0) Then
       fact=1.d0
    Else If (n .lt. 0) Then
       Write(*,*) 'n should .ge. 0 ', n
       Stop
    Else
       fact=1.d0
       Do i=1,n
          fact=fact* dble(i) 
       End Do
    End If
!    
    Return
  End function fact
!
!=========================================================================
! Subroutine to calculate the xy-integral -- M-R Version -- 
!=========================================================================
!     SubRoutine OAM_xy(Vxyz,Sxyz,nZeta,la,lb,nOrdOp,Alpha,Beta,Final,nComp)

!       Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp)
!       Complex*16 Vxyz(nZeta,3,0:la,0:lb,2)
!       Complex*16 Sxyz(nZeta,3,0:la+nOrdOp,0:lb+nOrdOp)
!       Real*8 Alpha(nZeta), Beta(nZeta)
!       Integer iZeta, i_x, i_y, i_z, j_x, j_y, j_z
!       Integer ix, iz, ixyz, ipa, ipb
!       Complex*16 Value1111, Value2111, Value0111, Value1211, Value1011, &
!                  Value1121, Value1101, Value1112, Value1110
!       Complex*16 Value1111_xA, Value1111_xB,                            &
!                  Value1111_yA, Value1111_yB
!       Real*8 Fact, rTemp
!       Complex*16 Temp1, Temp2, Temp
! !     Complex*16, Pointer :: zQxyz(:),zVxyz(:)
!       Integer nZeta, la, lb, nComp, nOrdOp, Ind

! !     Statement function for Cartesian index
! !
!       Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
! !
! !
! !     The integration over the z-subspace is done using the normal code
! !     for the exponetial operator.
! !
! #ifdef _DEBUG_
!       Write (6,*) 'OAM_xy: nOrdOp,nComp=',nOrdOp,nComp
!       Call RecPrt('Alpha',' ',Alpha,nZeta,1)
!       Call RecPrt('Beta',' ',Beta,nZeta,1)
!       Call RecPrt('Vxyz',' ',Vxyz,nZeta*2,3*(la+1)*(lb+1)*2)
!       Call RecPrt('Sxyz',' ',Sxyz,nZeta*2,3*(la+1+nOrdOp)*(lb+1+nOrdOp)*2)
!       Call RecPrt('Final',' ',Final,nZeta,                               &
!                   ((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)*nComp)
! #endif
!       Do iBeta = 1, nBeta
!          Do iAlpha = 1, nAlpha
!             iZeta = nAlpha*(iBeta-1) + iAlpha
! !
!             Do i_x =  0, la
!                Do i_y = 0, la-i_x
!                   i_z = la - i_x - i_y
!                   ipa=Ind(la,i_x,i_z)
! !
!                   Do j_x =  0, lb
!                      Do j_y = 0, lb-j_x
!                         j_z = lb - j_x - j_y
!                         ipb=Ind(lb,j_x,j_z)
! !
!                         rTemp=KVector(1)**2 + kVector(2)**2 + kVector(3)**2
!                         rTemp=rTemp/(Four*Zeta(iZeta))
! !
! !
! !      Compute the x-y part of the integral
! !
!                         Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                     Alpha(iZeta), Beta(iZeta),               &
!                                     i_x, i_y,                                &
!                                     j_x, j_y, lAng, Value1111)

!                         If (nOrdOp.eq.1) Then
!                            Fact = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two) * Exp(-rTemp)

!                            Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                        Alpha(iZeta), Beta(iZeta),               &
!                                        i_x+1, i_y,                              &
!                                        j_x, j_y, lAng, Value2111)

!                            Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                        Alpha(iZeta), Beta(iZeta),               &
!                                        i_x-1, i_y,                              &
!                                        j_x, j_y, lAng, Value0111)


!                            Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                        Alpha(iZeta), Beta(iZeta),               &
!                                        i_x, i_y,                                &
!                                        j_x+1,j_y, lAng, Value1211)

!                            Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                        Alpha(iZeta), Beta(iZeta),               &
!                                        i_x, i_y,                                &
!                                        j_x-1, j_y, lAng, Value1011)

!                            Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                        Alpha(iZeta), Beta(iZeta),               &
!                                        i_x, i_y+1,                              &
!                                        j_x, j_y, lAng, Value1121)
!                            Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                        Alpha(iZeta), Beta(iZeta),               &
!                                        i_x ,i_y-1,                              &
!                                        j_x, j_y, lAng, Value1101)

!                            Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                        Alpha(iZeta), Beta(iZeta),               &
!                                        i_x, i_y,                                &
!                                        j_x, j_y+1, lAng, Value1112)

!                            Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
!                                        Alpha(iZeta), Beta(iZeta),               &
!                                        i_x, i_y,                                &
!                                        j_x, j_y-1, lAng, Value1110)

! !                          Let us form the Cartesian intermediate intergral

!                            If (i_x.ne.0) Then
!                            Value1111_xA=DBLE(i_x) * Value0111                   &
!                                        - Two*Alpha(iZeta) * Value2111
!                            Else
!                            Value1111_xA= - Two*Alpha(iZeta) * Value2111
!                            End If
!                            If (j_x.ne.0) Then
!                            Value1111_xB=DBLE(j_x) * Value1011                   &
!                                        - Two* Beta(iZeta) * Value1211
!                            Else
!                            Value1111_xB= - Two* Beta(iZeta) * Value1211

!                            End If
!                            If (i_y.ne.0) Then
!                            Value1111_yA=DBLE(i_y) * Value0111                   &
!                                        - Two*Alpha(iZeta) * Value2111
!                            Else
!                            Value1111_yA= - Two*Alpha(iZeta) * Value2111

!                            End If
!                            If (j_y.ne.0) Then
!                            Value1111_yB=DBLE(j_y) * Value1011                   &
!                                        - Two* Beta(iZeta) * Value1211
!                            Else
!                            Value1111_yB= - Two* Beta(iZeta) * Value1211

!                            End If
! !
!                            Temp1 = Fact *                                       &
!                                    Value1111_xA * Sxyz(iZeta,3,i_z,j_z)
!                            Temp2 = Fact *                                       &
!                                    Value1111_xB * Sxyz(iZeta,3,i_z,j_z)
!                            Final(iZeta,ipa,ipb,1) = DBLE((Temp1+Temp2)*Half)
!                            Final(iZeta,ipa,ipb,4) = DBLE((Temp1-Temp2)*Half)
!                            Final(iZeta,ipa,ipb,7) = DIMAG((Temp1+Temp2)*Half)
!                            Final(iZeta,ipa,ipb,10)= DIMAG((Temp1-Temp2)*Half)
!                            Temp1 = Fact *                                       &
!                                    Value1111_yA * Sxyz(iZeta,3,i_z,j_z)
!                            Temp2 = Fact *                                       &
!                                    Value1111_yB * Sxyz(iZeta,3,i_z,j_z)
!                            Final(iZeta,ipa,ipb,2) = DBLE((Temp1+Temp2)*Half)
!                            Final(iZeta,ipa,ipb,5) = DBLE((Temp1-Temp2)*Half)
!                            Final(iZeta,ipa,ipb,8) = DIMAG((Temp1+Temp2)*Half)
!                            Final(iZeta,ipa,ipb,11)= DIMAG((Temp1-Temp2)*Half)
!                            Temp1 = Fact *                                       &
!                                    Value1111    * Vxyz(iZeta,3,i_z,j_z,1)
!                            Temp2 = Fact *                                       &
!                                    Value1111    * Vxyz(iZeta,3,i_z,j_z,2)
!                            Final(iZeta,ipa,ipb,3) = DBLE((Temp1+Temp2)*Half)
!                            Final(iZeta,ipa,ipb,6 )= DBLE((Temp1-Temp2)*Half)
!                            Final(iZeta,ipa,ipb,9 )= DIMAG((Temp1+Temp2)*Half)
!                            Final(iZeta,ipa,ipb,12)= DIMAG((Temp1-Temp2)*Half)
!                         Else
!                            Fact = rKappa(iZeta) * (1.0D0/Sqrt(Zeta(iZeta)**3))  &
!                                 * Exp(-rTemp)
!                            Temp1 = Fact *                                       &
!                                    Value1111    * Sxyz(iZeta,3,i_z,j_z)
!                            Final(iZeta,ipa,ipb,1) = DBLE(Temp)
!                            Final(iZeta,ipa,ipb,2) = DIMAG(Temp)
!                         End If
!                      End Do  ! j_y
!                   End Do     ! j_x
!                End Do  ! i_y
!             End Do     ! i_x
!          End Do  ! iAlpha
!       End Do     ! iBeta
! #ifdef _DEBUG_
!       Call RecPrt('Final',' ',Final,nZeta,                               &
!                   ((la+1)*(la+2)/2)*((lb+1)*(lb+2)/2)*nComp)
! #endif
! !
! !**********************************************************************
! !**********************************************************************
!     End SubRoutine OAM_xy
!======================================================================
End SubRoutine TWLInt
