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
      SubRoutine InpOne()
      use Arrays, only: CMO, Int1, KAIN1
      use OneDat, only: sOpSiz
      use rctfld_module
      use PCM_grad, only: iStpPCM,potnuc_pcm,DSCFAO,
     * PCM_grad_dens,PCM_grad_dens2,dscfmo,iCharge_PCM,
     * RFPERT
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Files_mclr.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "real.fh"
      Logical Do_ESPF,First,Dff,Do_DFT,NonEq
      Character*8 Label
      Integer iComp, idum(1)
      Real*8  rdum(1)
      Real*8, Allocatable:: D1ao(:), Nuc(:)
      Real*8, Allocatable:: Temp1(:), Temp2(:), Temp3(:)
      Real*8, Allocatable:: HTmp(:), GTmp(:)
      logical, external :: RF_On
      logical :: Found
*
      iRc=-1
      iOpt=ibset(0,sOpSiz)
      ndens2=0
      iisym=2**0
      Do iS=1,nSym
        nDens2=nDens2+Nbas(is)**2
      End Do
      Label='ONEHAM'
      iComp=1
      Call iRdOne(iRc,iOpt,Label,iComp,idum,iisym)
      leng=idum(1)
      If (iRC.ne.0)  Then
         Write (6,*) 'InpOne: Error reading ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
      iisym=2**0
      iRc=-1
      iOpt=0
      Call mma_allocate(Int1,ndens2,Label='Int1')
      kain1=>Int1

      Call mma_allocate(Temp1,leng+10,Label='Temp1')
      Call mma_allocate(Temp2,ndens2,Label='Temp2')
      Call mma_allocate(Temp3,ndens2,Label='Temp3')

      Call RdOne(iRc,iOpt,Label,iComp,Temp1,iisym)
      If (iRC.ne.0)  Then
         Write (6,*) 'InpOne: Error reading ONEINT'
         Write (6,'(A,A)') 'Label=',Label
         Call Abend()
      End If
C     write (6,*) "ONEHAM"
C     do i = 1, leng
C       write (6,'(i3,f20.10)') i,temp1(i)
C     end do
cnf
*
*     Modify the one electron Hamiltonian for reaction
*     field and ESPF calculations
*
      Tot_Nuc_Charge=Zero
      Call mma_allocate(Nuc,nAtoms,Label='Nuc')
      Call Get_dArray('Effective nuclear Charge',Nuc,nAtoms)
      Do iNuc = 1, nAtoms
        Tot_Nuc_Charge = Tot_Nuc_Charge + Nuc(iNuc)
      End Do
      Call mma_deallocate(Nuc)
      Tot_El_Charge = Zero
      Do iSym = 1, nSym
         Tot_El_Charge = Tot_El_Charge
     &                 - Two*DBLE(nFro(iSym)+nIsh(iSym))
      End Do
      Tot_El_Charge = Tot_El_Charge - DBLE(nActEl)
      Tot_Charge = Tot_Nuc_Charge + Tot_El_Charge
      iCharge = Int(Tot_Charge)
      Call DecideOnESPF(Do_ESPF)
      If (RF_On()) iCharge_PCM = iCharge
      If ( Do_ESPF .or. (lRF.and.iStpPCM==1)) then
*       If (lRF) Then
*          Write(6,*) 'Sorry, MCLR+RF NYI'
*          Call Quit_OnUserError()
*       End If
*
*------ Scratch for one- and two-electron type contributions
*------ + variational density-matrix
*
        If (RFPERT) then
*
*       Read the reaction field from RunFile or RunOld
*
          Call f_Inquire('RUNOLD',Found)
          If (Found) Call NameRun('RUNOLD')
          Call mma_allocate(Htmp,leng,Label='RCTFLD')
          Call Get_dScalar('RF Self Energy',ERFX)
          potnuc = potnuc + ERFX
          Call Get_dArray('Reaction field',Htmp,leng)
          Temp1 = Temp1 + Htmp
          Call mma_deallocate(Htmp)
          If (Found) Call NameRun('#Pop')
        Else
          Call mma_allocate(Htmp,leng,Label='Htmp')
          Call mma_allocate(Gtmp,leng,Label='Gtmp')
          Htmp(:)=Zero
          Gtmp(:)=Zero
          Call mma_allocate(D1ao,leng,Label='D1ao')
          Call Get_dArray_chk('D1ao',D1ao,leng)
          if (lRF) then
            !! D1ao is the state-specific (RlxRoot) density,
            !! but we need the density used for polarizing ASCs in SCF.
C           write (6,*) "d1ao before"
C           do i = 1, 10
C           write (6,'(i3,f20.10)') i,d1ao(i)
C           end do
            call PCM_grad_dens(1) ! SCF
            call PCM_grad_dens2(1,DSCFMO,DSCFAO)
            call fold(nSym,nBas,DSCFAO,D1ao)
            !! save the density for gradient
            Call Put_dArray('D1ao_PCM',D1ao,leng)
C           write (6,*) "d1ao after"
C           do i = 1, 10
C           write (6,'(i3,f20.10)') i,d1ao(i)
C           end do
          end if
*
          NonEq=.False.
          First=.True.
          Dff=.False.
          Do_DFT=.True.
          Call Get_dScalar('PotNuc',PotNuc)
          potnucsav = potnuc
          !! Htmp: nuclear-electron contributions
          !! Gtmp: electron-electron contributions
          Call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF,
*
*------ Don't care about the last arguments: no (CAS-)DFT here I guess)
*
     &               'SCF',Zero,iCharge,iSpin,rdum,rdum,0,'1234',Do_DFT)
          potnuc_pcm = potnuc - potnucsav
          !! not quite sure why e-e contributions are omitted
          !! in the original implementation
C         Call Daxpy_(leng,One,Htmp,1,Temp1,1)
          Temp1 = Temp1 + Htmp + Gtmp
*
*------ Hum, where the hell is FI (Fock Inactive) ???
*
*         Call Daxpy_(leng,One,Gtmp,1,FI,1)
          Call mma_deallocate(Gtmp)
          Call mma_deallocate(Htmp)
          Call mma_deallocate(D1ao)
        End If
      End If
C     write (6,*) "ONEHAM after DrvXV"
C     do i = 1, leng
C       write (6,'(i3,f20.10)') i,temp1(i)
C     end do
cnf
      ip=1
      ip2=1
      Do iS=1,nSym
        If (nBas(is).ne.0 .AND. nOrb(iS).ne.0) Then
           Call Square(Temp1(ip),
     &                   Temp2,
     &                   1,nBas(is),nBas(is))
           ip=ip+nBas(is)*(nBas(iS)+1)/2
           Call DGEMM_('T','N',
     &                 nOrb(iS),nBas(iS),nBas(iS),
     &                 1.0d0,CMO(ip2),nBas(iS),
     &                 Temp2,nBas(iS),
     &                 0.0d0,Temp3,nOrb(iS))
           Call DGEMM_('N','N',
     &                 nOrb(is),nOrb(iS),nBas(iS),
     &                 1.0d0,Temp3,nOrb(iS),
     &                 CMO(ip2),nBas(iS),
     &                 0.0d0,Int1(ip2),nOrb(iS))
C     write (6,*) "int1 in inpone.f"
C     call sqprt(int1(ip2),nbas(1))
           ip2=ip2+nBas(is)**2
        End If
      End Do
      Call mma_deallocate(Temp1)
      Call mma_deallocate(Temp2)
      Call mma_deallocate(Temp3)

      Return
      End
