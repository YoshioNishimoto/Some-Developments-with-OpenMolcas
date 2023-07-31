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
      Subroutine Get_DEcorr(nh1,Grad,nGrad,DFTFOCK)
      use SCF_Arrays, only: CMO
      use SpinAV
      use InfSCF
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8  Grad(nGrad), Ec_AB(2)
      Character*4 DFTFOCK
#include "addcorr.fh"
      Real*8, Allocatable :: F_DFT(:,:), D_DS(:,:)
      nD=2
*
      Call mma_allocate(F_DFT,nBT,nD,Label='F_DFT')
      Call mma_allocate(D_DS ,nBT,nD,Label='D_DS')
*
      Do iAB=1,2
       iOff=1
       jOff=1
       lOff=0
       Do iSym=1,nSym
          If (iAB.eq.1) Then
             nXoX=nOcc(iSym,1)
             iXoX=0
          Else
             nXoX=nConstr(iSym)
             iXoX=nOcc(iSym,1)-nConstr(iSym)
          EndIf
          mAdCMOO=iOff+nBas(iSym)*iXoX
          Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nXoX,
     &                     1.0d0,CMO(mAdCMOO,1),nBas(iSym),
     &                           CMO(mAdCMOO,1),nBas(iSym),
     &                     0.0d0,D_DS(jOff,1),nBas(iSym))
          If (iAB.eq.1) Then
             nXoX=nOcc(iSym,2)
             iXoX=0
          Else
             nXoX=nConstr(iSym)
             iXoX=nOcc(iSym,2)-nConstr(iSym)
          EndIf
          mAdCMOO=iOff+nBas(iSym)*iXoX
          Call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nXoX,
     &                     1.0d0,CMO(mAdCMOO,2),nBas(iSym),
     &                           CMO(mAdCMOO,2),nBas(iSym),
     &                     0.0d0,D_DS(jOff,2),nBas(iSym))
*
          If (Do_SpinAV) Then
             Do j=1,nBas(iSym)
                Do i=1,j
                   iDSc=nBas(iSym)*(j-1)+i
                   ji=j*(j-1)/2+i
                   iDij=jOff-1+ji
                   D_DS(iDij,1)=D_DS(iDij,1)-DSc(iDSc)
                   D_DS(iDij,2)=D_DS(iDij,2)+DSc(iDSc)
                End Do
             End Do
             lOff=lOff+nBas(iSym)**2
          EndIf
*
          Do j=1,nBas(iSym)
             Do i=1,j-1
                ji=j*(j-1)/2+i
                iDij=jOff-1+ji
                D_DS(iDij,1)=2.0d0*D_DS(iDij,1)
                D_DS(iDij,2)=2.0d0*D_DS(iDij,2)
             End Do
          End Do
          iOff=iOff+nBas(iSym)*nOrb(iSym)
          jOff=jOff+nBas(iSym)*(nBas(iSym)+1)/2
       End Do
*
       Call Get_Ecorr_dft(nh1,Grad,nGrad,DFTFOCK,F_DFT,D_DS,nBT,nD,
     &                        ADDC_KSDFT,Ec_AB(iAB))
      End Do
*----------------------------------------------------------------------*
      DE_KSDFT_c=Ec_AB(1)-Ec_AB(2)
*----------------------------------------------------------------------*
*
      Call mma_deallocate(D_DS)
      Call mma_deallocate(F_DFT)
      Return
      End
*                                                                      *
************************************************************************
*                                                                      *
      Subroutine Get_Ecorr_dft(nh1,Grad,nGrad,DFTFOCK,F_DFT,D_DS,
     &                         nBT,nD,KSDFT,Ec_AB)
      use OFembed, only: dFMD, Do_Core
      use nq_Info, only: Dens_I, Grad_I, Tau_I
      Implicit Real*8 (a-h,o-z)

#include "real.fh"
#include "debug.fh"
      Real*8  Grad(nGrad)
      Logical Do_MO,Do_TwoEl,Do_Grad
      Character(LEN=4) DFTFOCK
      Character(LEN=80)  KSDFT
      Real*8 :: F_DFT(nBT,nD), D_DS(nBT,nD)

      Debug=.False.
*                                                                      *
************************************************************************
*                                                                      *
*     DFT functionals, compute integrals over the potential
*
      Func            =Zero
      Dens_I          =Zero
      Grad_I          =Zero
      Tau_I           =Zero
      Do_MO           =.False.
      Do_TwoEl        =.False.
      Do_Grad=.false.
*
      nFckDim=2
      dFMD_=dFMD
      dFMD=1.0d0
*                                                                      *
************************************************************************
*                                                                      *
      Do_Core=.True.
      Call Driver(KSDFT,Do_Grad,Func,Grad,nGrad,
     &            Do_MO,Do_TwoEl,D_DS,F_DFT,nh1,nFckDim,DFTFOCK)
      Do_Core=.False.
*                                                                      *
************************************************************************
*                                                                      *
      Ec_AB=Func
*
#ifdef _DEBUGPRINT_
      write(6,*) ' Correlation energy: ',Ec_AB
      write(6,*)
      write(6,*) ' Correlation potentials: (itri,F_alpha,F_beta)'
      write(6,*)
      Do i=1,nh1
        Write(6,'(i4,3f22.16)') i,F_DFT(i,1),F_DFT(i,2)
      End Do
#endif
*
      dFMD=dFMD_
      Return
      End
