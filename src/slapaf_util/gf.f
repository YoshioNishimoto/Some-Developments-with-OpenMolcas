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
      Subroutine GF(nX,nDoF,nInter,EVec,EVal,RedM,iNeg,dDipM,mTR,nAtom,
     &              DipM)
      use Slapaf_Info, only: Smmtrc
      use Slapaf_parameters, only: nDimBC
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 dDipM(3,nInter+mTR), DipM(3),
     &       EVec(2*nDoF,nDoF), EVal(2*nDoF), RedM(nDoF)
      Real*8, Allocatable:: G(:), GInv(:), F(:), Tmp1(:), Tmp2(:)
*                                                                      *
************************************************************************
*                                                                      *
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Call RecPrt('GF: dDipM',' ',dDipM,3,nInter)
      Call RecPrt('GF: DipM',' ',DipM,3,1)
#endif
      Call mma_allocate(Tmp1,nX**2,Label='Tmp1')
      Call mma_allocate(Tmp2,nX**2,Label='Tmp2')
*                                                                      *
************************************************************************
*                                                                      *
*     Note that all calculations will be done in the Cartesian basis!  *
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute harmonic frequencies
*                                                                      *
************************************************************************
*                                                                      *
*     Generate the G matrix (mass tensor in cartesians)
*
      Call mma_allocate(G,nDimBC**2,Label='G')
      Call mma_allocate(GInv,nDimBC**2,Label='GInv')
      Call Mk_G(G,GInv,nDimBC)
*
*
*     Get the force constant matrix in cartesians
*
      Call mma_allocate(F,nX**2,Label='F')
      Call Get_H(F,nX)
*                                                                      *
************************************************************************
*                                                                      *
*     Form the GF-matrix (actually G^(1/2)FG^(1/2))
*
      Call GF_Mult(G,F,Tmp2,nDoF)  ! Result in Tmp2
      Call mma_deallocate(F)
*
*     Compute the frequencies and harmonic eigenfunctions in
*     Cartesians.
*
      Call GF_Harmonic_Frequencies(G,GInv,Tmp1,Tmp2,EVec,EVal,RedM,
     &                             iNeg,nX,nDoF)
*
      Call mma_deallocate(G)
      Call mma_deallocate(GInv)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Compute the dipole moment derivative in Cartesians.
*
      Call Get_dDipM(dDipM,DipM,nDoF,nInter)
*                                                                      *
************************************************************************
*                                                                      *
*    Transform from cartesian to normal coordinates
*
      Do iNC = 1, nDoF
         call dcopy_(nDoF,EVec(1,iNC),2,Tmp2,1)
         ix = (iNC-1)*3 + 1
         iy = (iNC-1)*3 + 2
         iz = (iNC-1)*3 + 3
         Tmp1(ix) = Zero
         Tmp1(iy) = Zero
         Tmp1(iz) = Zero
         i = 0
         Do iAtom = 1, nAtom
            Do ixyz = 1, 3
               If (Smmtrc(ixyz,iAtom)) Then
                  i = i + 1
                  Tmp1(ix) = Tmp1(ix) + dDipM(1,i)*Tmp2(i)
                  Tmp1(iy) = Tmp1(iy) + dDipM(2,i)*Tmp2(i)
                  Tmp1(iz) = Tmp1(iz) + dDipM(3,i)*Tmp2(i)
               End If
            End Do
         End Do
      End Do
      call dcopy_(3*nDoF,Tmp1,1,dDipM,1)
#ifdef _DEBUGPRINT_
      Call RecPrt('dDipM(normal coord.)',' ',dDipM,3,nDoF)
#endif
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Tmp1)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
