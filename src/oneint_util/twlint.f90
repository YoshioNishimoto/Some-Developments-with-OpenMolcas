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
             Zeta(nZeta), Alpha(nAlpha), Beta(nBeta),                  &
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
      Integer iAlpha, iBeta, ixyz
      Integer ipA, ipAOff, ipAxyz, ipB, ipBOff, ipBxyz, ipQxyz, ipRes, &
              ipVxyz, nip, icomp, iOper, lDCRT, llOper, LmbdT, nDCRT,  &
              nIrrep, nOp, nStabO, ipP, lAng, ipScr, iOff
      Real*8 Zero, Half, One, Two, Three, Four, Rxy, Fi1, Rxyz, Fi2, Fi3
      Real*8, Allocatable:: TransM(:,:)
      Real*8 kVector_local(3)
      Real*8 A_Local(3), RB_Local(3)

      Integer :: i, l
      Real*8 :: dmm

!      Real(kind=8), dimension(dim,dim) :: a_mtrx, inv_mtrx, RSph_inv
      Integer :: i, im, jm
!
      Zero =0.0D0
      Half =0.5D0
      One  =1.0D0
      Two  =2.0D0
      Three=3.0D0
      Four =4.0D0
      lAng= 1             ! Temporary set value
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
      If (lAng.eq.0) Then
         kVector_Local(1)=kVector(1)
         kVector_Local(2)=kVector(2)
         kVector_Local(3)=kVector(3)
         Go To 114
      End If
!
!     Use Euler angles (z-x-z). We skip the last rotation.
!     Rotate around the z-axis so that the x-component becomes zero.
      Rxy=Sqrt(kvector(1)**2 + kvector(2)**2)
      Fi1= ATAN2(kvector(2),kvector(1))
      kVector_Local(1)=Rxy
      kVector_Local(2)=Zero
      kVector_Local(3)=kVector(3)
!
!     Rotate around the x-axis so that the y-component becomes zero.
      Rxyz=Sqrt(kVector_Local(1)**2 + kVector_Local(3)**2)
      Fi2 = ATAN2(kVector_Local(3),kVector_Local(1))
      kVector_Local(1)=Zero
      kVector_Local(2)=Zero
      kVector_Local(3)=Rxyz
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
      Call DGEMM_('N','T',3,3,nZeta,                                    &
                   1.0D0,P,nZeta,                                       &
                         TransM,3,                                      &
                   0.0D0,Array(ipP),3)
      Call mma_deallocate(TransM)
!
 114  Continue
      If (lAng.eq.0) Then
         Call TWLInt_Internal(Array,A,RB,P)
      Else
         Call TWLInt_Internal(Array,A_Local,RB_Local,Array(ipP))
      End If
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!====================================================================
!
      If (lAng.eq.0) Go To 263
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
      Call mma_Allocate(TransM,(la+1)*(la+2)/2,(la+1)*(la+2)/2,Label="TransM")
      TransM(:,:)=0.0D0
      iOff=1
      Do i = la,0,-2
!==============================================================================
!        Generate the blocks of the transformation matrix and put them into
!        TransM
         Call dmm_transform(i,-Fi3,-Fi2,-Fi1,TransM(iOff,iOff),            &
                            (la+1)*(la+2)/2))
         iOff = iOff + 2*i + 1
!
!==============================================================================
      End Do
      Call DGEMM_('T','T',                                                 &
                 (la+1)*(la+2)/2,((lb+1)*(lb+2)/2)*nZeta,(la+1)*(la+2)/2,  &
                 1.0D0,TransM,(la+1)*(la+2)/2,                             &
                       Array(ipRes),((lb+1)*(lb+2)/2)*nZeta,               &
                 0.0D0,Array(ipScr),(la+1)*(la+2)/2)
      Call mma_deallocate(TransM)
!
!        B(new),ij,A(orig) -> ij,A(orig),B(orig)
      Call mma_Allocate(TransM,(lb+1)*(lb+2)/2,(lb+1)*(lb+2)/2,Label="TransM")
      Do i = lb,0,-2
!==============================================================================
!        Generate the blocks of the transformation matrix and put them into
!        TransM
         Call dmm_transform(...)
!==============================================================================
      End Do
      Call DGEMM_('T','T',                                                 &
                 (lb+1)*(lb+2)/2,((la+1)*(la+2)/2)*nZeta,(lb+1)*(lb+2)/2,  &
                  1.0D0,TransM,(lb+1)*(lb+2)/2,                            &
                        Array(ipScr),nZeta*((la+1)*(la+2)/2),              &
                 0.0D0,Array(ipRes),(lb+1)*(lb+2)/2)

      Call mma_deallocate(TransM)
!
!     3) transform the spherical harmonics to the Cartesians.
!
!        Find inverse for RSph(ipSph(ld))
!====================================================================
!       ij,A,B -> b,ij,A
        Call DGEMM_('T','T',                                                &
                   (ld+1)*(ld+2)/2,nZeta*((la+1)*(la+2)/2),(lb+1)*(lb+2)/2, &
                   1.0D0,RSph(ipSph(lb)),((lb+1)*(lb+2)/2),                 &
                         Array(ipRes),nZeta*((la+1)*(la+2)/2),              &
                   0.0D0,Array(ipScr),(lb+1)*(lb+2)/2 )

!       Find inverse for RSph(ipSph(lc))
!       b,ij,A -> a,b,ij
        Call DGEMM_('T','T',                                                &
                   (lc+1)*(lc+2)/2,nZeta*((ld+1)*(ld+2)/2),(lc+1)*(lc+2)/2, &
                   1.0D0,RSph(ipSph(lc)),((lc+1)*(lc+2)/2),                 &
                         Array(ipScr),((ld+1)*(ld+2)/2)*nZeta,              &
                   0.0D0,Array(ipRes),(lc+1)*(lc+2)/2 )
!
!     4) transpose to the correct order
!        c,d,ij -> ij,c,d  make sure that it is in Array(ipRes)
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
      Real*8 Alpha_, Beta_, rKappa(nZeta)
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
         Call WarningMessage(2,'EMFInt: nip-1.gt.nArr*nZeta')
         Write (6,*) ' nArr is Wrong! ', nip-1,' > ',nArr*nZeta
         Write (6,*) ' Abend in EMFInt'
         Call Abend()
      End If
!
#ifdef _DEBUG_
      Call RecPrt(' In EMFInt: A',' ',A,1,3)
      Call RecPrt(' In EMFInt: RB',' ',RB,1,3)
      Call RecPrt(' In EMFInt: KVector',' ',kvector_Local,1,3)
      Call RecPrt(' In EMFInt: P',' ',P,nZeta,3)
      Write (6,*) ' In EMFInt: la,lb=',la,lb
#endif
!
!**********************************************************************
!    Computation of intermediate integrals
!
!    Compute all Cartesian components with the OAM free code. Then
!    compute the x and y components in the OAM formalism replacing the
!    OAM-free x and y components.
!
      Call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
!
!
!     Compute the Cartesian values of the basis functions angular part
!     Note that these arrays are complex.
!
      Call C_F_Pointer(C_Loc(Array(ipAxyz)),zAxyz,[nZeta*3*nHer*(la+nOrdOp+1)])
      Call CCrtCmp(Zeta,P,nZeta,A,zAxyz,la+nOrdOp,HerR(iHerR(nHer)),nHer,ABeq,kvector_Local)
      Call C_F_Pointer(C_Loc(Array(ipBxyz)),zBxyz,[nZeta*3*nHer*(lb+nOrdOp+1)])
      Call CCrtCmp(Zeta,P,nZeta,RB,zBxyz,lb+nOrdOp,HerR(iHerR(nHer)),nHer,ABeq,kvector_Local)
      Nullify(zAxyz,zBxyz)
!
!     Compute the Cartesian components for the multipole moment
!     integrals. The integrals are factorized into components.
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
!
!     Compute the cartesian components for the velocity integrals.
!     The velocity components are linear combinations of overlap
!     components.
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
      If (lAng.ne.0) Then
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
!
!     At this point the integrals are store in Array starting at the position
!     ipRes. The size of the block is nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2
!
    End SubRoutine TWLInt_Internal
!=========================================================================
!
    Integer Function nElem(ixyz)
      Integer ixyz
      nElem = (ixyz+1)*(ixyz+2)/2
    End Function nElem
!
!=========================================================================
    SubRoutine OAM_xy(Vxyz,Sxyz,nZeta,la,lb,nOrdOp,Alpha,Beta,Final,nComp)

      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nComp)
      Complex*16 Vxyz(nZeta,3,0:la,0:lb,2)
      Complex*16 Sxyz(nZeta,3,0:la+nOrdOp,0:lb+nOrdOp)
      Real*8 Alpha(nZeta), Beta(nZeta)
      Integer iZeta, i_x, i_y, i_z, j_x, j_y, j_z
      Integer ix, iz, ixyz, ipa, ipb
      Complex*16 Value1111, Value2111, Value0111, Value1211, Value1011, &
                 Value1121, Value1101, Value1112, Value1110
      Complex*16 Value1111_xA, Value1111_xB,                            &
                 Value1111_yA, Value1111_yB
      Real*8 Fact, rTemp
      Complex*16 Temp1, Temp2, Temp
      Complex*16, Pointer :: zQxyz(:),zVxyz(:)
      Integer nZeta, la, lb, nComp, nOrdOp, Ind

!     Statement function for Cartesian index
!
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!
!
!     The integration over the z-subspace is done using the normal code
!     for the exponetial operator.
!
      Do iBeta = 1, nBeta
         Do iAlpha = 1, nAlpha
            iZeta = nAlpha*(iBeta-1) + iAlpha
!
            Do i_x =  0, la
               Do i_y = 0, la-i_x
                  i_z = la - i_x - i_y
                  ipa=Ind(la,i_x,i_z)
!
                  Do j_x =  0, lb
                     Do j_y = 0, lb-j_x
                        j_z = lb - j_x - j_y
                        ipb=Ind(lb,j_x,j_z)
!
                        rTemp=KVector(1)**2 + kVector(2)**2 + kVector(3)**2
                        rTemp=rTemp/(Four*Zeta(iZeta))
!
!
!      Compute the x-y part of the integral
!
                        Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                    Alpha(iZeta), Beta(iZeta),               &
                                    i_x, i_y,                                &
                                    j_x, j_y, lAng, Value1111)

                        If (nOrdOp.eq.1) Then
                           Fact = rKappa(iZeta) * Zeta(iZeta)**(-Three/Two) * Exp(-rTemp)

                           Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                       Alpha(iZeta), Beta(iZeta),               &
                                       i_x+1, i_y,                              &
                                       j_x, j_y, lAng, Value2111)

                           Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                       Alpha(iZeta), Beta(iZeta),               &
                                       i_x-1, i_y,                              &
                                       j_x, j_y, lAng, Value0111)


                           Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                       Alpha(iZeta), Beta(iZeta),               &
                                       i_x, i_y,                                &
                                       j_x+1,j_y, lAng, Value1211)

                           Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                       Alpha(iZeta), Beta(iZeta),               &
                                       i_x, i_y,                                &
                                       j_x-1, j_y, lAng, Value1011)

                           Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                       Alpha(iZeta), Beta(iZeta),               &
                                       i_x, i_y+1,                              &
                                       j_x, j_y, lAng, Value1121)
                           Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                       Alpha(iZeta), Beta(iZeta),               &
                                       i_x ,i_y-1,                              &
                                       j_x, j_y, lAng, Value1101)

                           Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                       Alpha(iZeta), Beta(iZeta),               &
                                       i_x, i_y,                                &
                                       j_x, j_y+1, lAng, Value1112)

                           Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                                       Alpha(iZeta), Beta(iZeta),               &
                                       i_x, i_y,                                &
                                       j_x, j_y-1, lAng, Value1110)

!                          Let us form the Cartesian intermediate intergral

                           If (i_x.ne.0) Then
                           Value1111_xA=DBLE(i_x) * Value0111                   &
                                       - Two*Alpha(iZeta) * Value2111
                           Else
                           Value1111_xA= - Two*Alpha(iZeta) * Value2111
                           End If
                           If (j_x.ne.0) Then
                           Value1111_xB=DBLE(j_x) * Value1011                   &
                                       - Two* Beta(iZeta) * Value1211
                           Else
                           Value1111_xB= - Two* Beta(iZeta) * Value1211

                           End If
                           If (i_y.ne.0) Then
                           Value1111_yA=DBLE(i_y) * Value0111                   &
                                       - Two*Alpha(iZeta) * Value2111
                           Else
                           Value1111_yA= - Two*Alpha(iZeta) * Value2111

                           End If
                           If (j_y.ne.0) Then
                           Value1111_yB=DBLE(j_y) * Value1011                   &
                                       - Two* Beta(iZeta) * Value1211
                           Else
                           Value1111_yB= - Two* Beta(iZeta) * Value1211

                           End If
!
                           Temp1 = Fact *                                       &
                                   Value1111_xA * Sxyz(iZeta,3,i_z,j_z)
                           Temp2 = Fact *                                       &
                                   Value1111_xB * Sxyz(iZeta,3,i_z,j_z)
                           Final(iZeta,ipa,ipb,1) = DBLE((Temp1+Temp2)*Half)
                           Final(iZeta,ipa,ipb,4) = DBLE((Temp1-Temp2)*Half)
                           Final(iZeta,ipa,ipb,7) = DIMAG((Temp1+Temp2)*Half)
                           Final(iZeta,ipa,ipb,10)= DIMAG((Temp1-Temp2)*Half)
                           Temp1 = Fact *                                       &
                                   Value1111_yA * Sxyz(iZeta,3,i_z,j_z)
                           Temp2 = Fact *                                       &
                                   Value1111_yB * Sxyz(iZeta,3,i_z,j_z)
                           Final(iZeta,ipa,ipb,2) = DBLE((Temp1+Temp2)*Half)
                           Final(iZeta,ipa,ipb,5) = DBLE((Temp1-Temp2)*Half)
                           Final(iZeta,ipa,ipb,8) = DIMAG((Temp1+Temp2)*Half)
                           Final(iZeta,ipa,ipb,11)= DIMAG((Temp1-Temp2)*Half)
                           Temp1 = Fact *                                       &
                                   Value1111    * Vxyz(iZeta,3,i_z,j_z,1)
                           Temp2 = Fact *                                       &
                                   Value1111    * Vxyz(iZeta,3,i_z,j_z,2)
                           Final(iZeta,ipa,ipb,3) = DBLE((Temp1+Temp2)*Half)
                           Final(iZeta,ipa,ipb,6 )= DBLE((Temp1-Temp2)*Half)
                           Final(iZeta,ipa,ipb,9 )= DIMAG((Temp1+Temp2)*Half)
                           Final(iZeta,ipa,ipb,12)= DIMAG((Temp1-Temp2)*Half)
                        Else
                           Fact = rKappa(iZeta) * (1.0D0/Sqrt(Zeta(iZeta)**3))  &
                                * Exp(-rTemp)
                           Temp1 = Fact *                                       &
                                   Value1111    * Sxyz(iZeta,3,i_z,j_z)
                           Final(iZeta,ipa,ipb,1) = DBLE(Temp)
                           Final(iZeta,ipa,ipb,2) = DIMAG(Temp)
                        End If
                     End Do  ! j_y
                  End Do     ! j_x
               End Do  ! i_y
            End Do     ! i_x
         End Do  ! iAlpha
      End Do     ! iBeta
!
!**********************************************************************
!**********************************************************************
    End SubRoutine OAM_xy
!======================================================================
End SubRoutine TWLInt
