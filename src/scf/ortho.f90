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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************
      SubRoutine Ortho(AMat,nAMat,Ovlp,nOvlp)
!***********************************************************************
!                                                                      *
!     purpose: Transform to orthonormal basis (one symmetry block      *
!              at a time)                                              *
!     input:                                                           *
!       AMat    : transformation matrix to non-orthonormal basis of    *
!                 length nAMat                                         *
!       Ovlp    : overlap in AO basis of length nOvlp                  *
!                                                                      *
!     output:                                                          *
!       AMat    : transformation matrix to orthonormal basis           *
!                                                                      *
!***********************************************************************
      use InfSCF, only: MaxBas, MaxBxO, MaxOrb, nSym, nOrb, nBas
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nAMat, nOvlp
      Real*8 AMat(nAMat),Ovlp(nOvlp)

      Real*8, Dimension(:), Allocatable:: OvlT, OvlH, OvlS
      Integer iiBO, iiBT, ij, im, iSym
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!---- Allocate memory for transformed overlap matrix
      Call mma_allocate(OvlT,MaxOrb**2,Label='OvlT')
!
!---- Allocate memory for half-transformed overlap matrix
      Call mma_allocate(OvlH,MaxBxO,Label='OvlH')
!
!---- Allocate memory for squared overlap matrix
      Call mma_allocate(OvlS,MaxBas**2,Label='OvlS')
!
      ij = 1
      im = 1
      Do iSym = 1, nSym
!
         iiBO = nBas(iSym)*nOrb(iSym)
         iiBT = nBas(iSym)*(nBas(iSym) + 1)/2
!
         If (nOrb(iSym).gt.0) Then
!
!---------- Square overlap and transform to the basis given by AMat
!           Call Square(Ovlp(ij),nBas(iSym),OvlS,nBas(iSym))
            Call Square(Ovlp(ij),OvlS,1,nBas(iSym),nBas(iSym))
            Call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym),   &
                        One,OvlS,nBas(iSym),                        &
                              AMat(im),nBas(iSym),                  &
                        Zero,OvlH,nBas(iSym))
            Call DGEMM_('T','N',nOrb(iSym),nOrb(iSym),nBas(iSym),   &
                        One,AMat(im),nBas(iSym),                    &
                              OvlH,nBas(iSym),                      &
                        Zero,OvlT,nOrb(iSym))
!
!---------- Orthogonalize (Gram-Schmidt)
            Call Orthox(OvlT,AMat(im),nOrb(iSym),nBas(iSym))
!
         End If
!
!------- Update pointers
         im = im + iiBO
         ij = ij + iiBT
!
      End Do
!
!---- Deallocate memory
      Call mma_deallocate(OvlT)
      Call mma_deallocate(OvlH)
      Call mma_deallocate(OvlS)
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
      End
