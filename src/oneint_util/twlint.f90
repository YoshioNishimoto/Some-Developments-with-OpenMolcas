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
! Copyright (C) 2019 Marjan Khamesian and Roland Lindh                 *
!***********************************************************************
      SubRoutine TWLInt(Alpha,nAlpha,Beta, nBeta,Zeta,ZInv,rKappa,P,   &
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
      Implicit None
!
!     External Arrays and integers
!
      Integer nZeta, la, lb, nIC, nAlpha, nBeta, nArr, nComp, nOrdOp,  &
              nStabM
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,nIC),         &
             Zeta(nZeta), ZInv(nZeta), Alpha(nAlpha), Beta(nBeta),     &
             rKappa(nZeta), P(nZeta,3), A(3), RB(3),                   &
             Array(nZeta*nArr), kvector(3)
      Integer iStabM(0:nStabM-1), iDCRT(0:7),                          &
              iStabO(0:7), lOper(nComp), iChO(nComp)
      Integer nHer, nGrid, iAddPot
      Real*8  PtChrg
      Logical ABeq(3)
!
!     Local Arrays and integers
!
      Integer iZeta, iAlpha, iBeta, i_x, i_y, i_z, j_x, j_y, j_z,      &
              lAng, ixyz
      Integer ipA, ipAOff, ipAxyz, ipB, ipBOff, ipBxyz, ipQxyz, ipRes, &
              ipVxyz, nip
      Real*8 Value, Zero, Rxy, Fi1, Rxyz, Fi2
      Real*8 TransM(3,3)
      Real*8 kVector_local(3)
!
      Zero=0.0D0
      lAng= 1             ! Temporary set value
!
!     We need here a transformation of the Cartesian coordinates in which
!     the k-vector coinsides with the z-vector direction.
!     In particular, we will transform P to the new coordinate system.
!
!     Use Euler angles (z-x-z). We skip the last rotation.
!     Rotate around the z-axis so that the x-component becomes zero.
      Rxy=Sqrt(kvector(1)**2 + kvector(2)**2)
      Fi1= - ATAN2(kvector(2),kvector(1))
      kVector_Local(1)=Rxy
      kVector_Local(2)=Zero
      kVector_Local(3)=kVector(3)
!
!     Rotate around the x-axis so that the y-componen becomes zero.
      Rxyz=Sqrt(kVector_Local(1)**2 + kVector_Local(3)**2)
      Fi2 = - ATAN2(kVector_Local(3),kVector_Local(1))
      kVector_Local(1)=Zero
      kVector_Local(2)=Zero
      kVector_Local(3)=Rxyz
!
      TransM(1,1)= Cos(Fi1)
      TransM(2,1)= Sin(Fi1)
      TransM(3,1)= Zero
      TransM(1,2)=-Sin(Fi1)*Cos(Fi2)
      TransM(2,2)= Cos(Fi1)*Cos(Fi2)
      TransM(2,2)=          Sin(Fi2)
      TransM(1,3)= Sin(Fi1)*Sin(Fi2)
      TransM(2,3)=-Cos(Fi1)*Sin(Fi2)
      TransM(3,3)=          Cos(Fi2)

!     ... more to come ...
!
      Call TWLInt_Internal(Array)
!
!**********************************************************************

!      Now when all Cartesian components have been computed we
!      transform back to the coordinate system of the molecule.

!      ... more to come ...
!
      Return

! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(PtChrg)
         Call Unused_integer(nGrid)
         Call Unused_integer(iAddPot)
      End If
!**********************************************************************
!**********************************************************************
!
!     This is to allow type punning without an explicit interface
      Contains
      SubRoutine TWLInt_Internal(Array)
      Use Iso_C_Binding
      Real*8, Target :: Array(*)
      Complex*16, Pointer :: zAxyz(:),zBxyz(:),zQxyz(:),zVxyz(:)
!
      ABeq(1) = A(1).eq.RB(1)
      ABeq(2) = A(2).eq.RB(2)
      ABeq(3) = A(3).eq.RB(3)
!
      nip = 1
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
      Call RecPrt(' In EMFInt: KVector',' ',kvector,1,3)
      Call RecPrt(' In EMFInt: P',' ',P,nZeta,3)
      Write (6,*) ' In EMFInt: la,lb=',la,lb
#endif
!
!**********************************************************************
!    Computation of intermediate integrals
!
!    Compute all Cartesian components with the OAM free code. Then
!    compute the z-component only in the OAM formalism replacing the
!    OAM-free z-components.
!
      call dcopy_(nZeta*nElem(la)*nElem(lb)*nIC,[Zero],0,Final,1)
!
!     Compute the cartesian values of the basis functions angular part
!     Note that these arrays are complex.
!
      Call C_F_Pointer(C_Loc(Array(ipAxyz)),zAxyz,[nZeta*3*nHer*(la+nOrdOp+1)])
      Call CCrtCmp(Zeta,P,nZeta,A,zAxyz,la+nOrdOp,HerR(iHerR(nHer)),nHer,ABeq,kvector)
      Call C_F_Pointer(C_Loc(Array(ipBxyz)),zBxyz,[nZeta*3*nHer*(lb+nOrdOp+1)])
      Call CCrtCmp(Zeta,P,nZeta,RB,zBxyz,lb+nOrdOp,HerR(iHerR(nHer)),nHer,ABeq,kvector)
      Nullify(zAxyz,zBxyz)
!
!     Compute the cartesian components for the multipole moment
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
!     Now compute the OAM z-component.
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

            Do j_x =  0, lb
               Do j_y = 0, lb-j_x
                  j_z = lb - j_x - j_y

!                 Compute the x-y part of the integral

                  Call twlprm(Zeta(iZeta),P(iZeta,1),P(iZeta,2),       &
                              Alpha(iAlpha), Beta(iBeta),              &
                              i_x,i_y,                                 &
                              j_x,j_y, lAng, Value)

!                 Compute the z part of the integral

!                 ... more to come ...

!                 Assemble to the complete integral

!                 ... more to come ...


               End Do ! j_y
            End Do    ! j_x

               End Do    ! i_y
            End Do    ! i_x

         End Do ! iAlpha
      End Do    ! iBeta
!**********************************************************************
!
!     Do not forget to merge the z-component into the right place.
!
!**********************************************************************
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
!
!     Combine the cartesian components to the full one electron
!     integral.
!
         Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,[nZeta*3*(la+1)*(lb+1)])
         Call C_F_Pointer(C_Loc(Array(ipVxyz)),zVxyz,[nZeta*3*(la+1)*(lb+1)*2])
         Call CCmbnVe(zQxyz,nZeta,la,lb,Zeta,rKappa,Array(ipRes),nComp,zVxyz,kvector)
         Nullify(zQxyz,zVxyz)
      Else
         Call C_F_Pointer(C_Loc(Array(ipQxyz)),zQxyz,[nZeta*3*(la+1)*(lb+1)*(nOrdOp+1)])
         Call CCmbnMP(zQxyz,nZeta,la,lb,nOrdOp,Zeta,rKappa,Array(ipRes),nComp,kvector)
         Nullify(zQxyz)
      End If
!
      End SubRoutine TWLInt_Internal
      Integer Function nElem(ixyz)
      Integer ixyz
!
      nElem = (ixyz+1)*(ixyz+2)/2
!
      End Function nElem
!**********************************************************************
!**********************************************************************
!
      End Subroutine twlint
