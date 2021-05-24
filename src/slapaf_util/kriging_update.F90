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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************
Subroutine Kriging_Update(nQQ,iter,qInt,E_Disp)
Use Slapaf_Info, only: Energy, dqInt, Energy0, BMx_Save, Gx, Gx0, Degen, NAC
Use Kriging_Mod, only: nSet, iter_actual
Use Slapaf_Parameters, only: Curvilinear
Implicit None
Integer nQQ, iter
Real*8  qInt(nQQ), E_Disp

#include "real.fh"
#include "stdalloc.fh"
Integer iSet, nAtoms, iAtom, ixyz
Real*8  Temp(3), Demp(3)
Real*8, Allocatable:: Aux(:,:)


Call mma_allocate(Aux,nQQ,nSet,Label='Aux')
!                                                                      *
!***********************************************************************
!                                                                      *
!   This routine will interface between the relaxation code and the
!   kriging code. Give a new set of internal coordinates, qInt(:), it
!   will update lists containing values, and corresponding gradients
!   as if an iteration on the parent model had been performed. Note that
!   depending on the case BOTH gradients in internal and Cartesian
!   coordinates are done. All this to mimic that a true macro iteration
!   had been performed.
!                                                                      *
!***********************************************************************
!
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
Write (6,*) 'Kriging_Update, iter, iter_actual=', iter, iter_actual
Call RecPrt('Kriging_Update: qInt',' ',qInt,nQQ,1)
#endif

Call Energy_Kriging_layer(qInt(:),Temp,nQQ)

Call Dispersion_Kriging_Layer(qInt(:),Demp,nQQ)

Call Gradient_Kriging_layer(qInt(:),Aux,nQQ)

#ifdef _DEBUGPRINT_
Call RecPrt('Kriging_Update: Temp',' ',Temp,1,nSet)
Call RecPrt('Kriging_Update: Demp',' ',Demp,1,nSet)
Call RecPrt('Kriging_Update: Aux',' ',Aux,nQQ,nSet)
#endif

Aux(:,:) = - Aux(:,:) ! Change the sign of the gradient
!                                                                      *
!***********************************************************************
!
!  The default case
!
iSet = 1
Energy(iter) = Temp(iSet)

E_Disp = Demp(iSet)

dqInt(:,iter) = Aux(:,iSet)

If (nSet==1) Then
   Call mma_deallocate(Aux)
   Return
End If
!                                                                      *
!***********************************************************************
!                                                                      *
!  For the energy difference
!
iSet = 2
Energy0(iter+iter_actual-1) = Temp(iSet)

! Right now we do not use the dispersion for the constraints

! The computation of the gradients of the constraints are always done in
! Cartesian coordinates. The GEK, however, predict them in internal coordinates.
! Thus, we need to transfrom the GEK predicted gradients to Cartesian and put
! them at the place where the code to compute the value and the Cartesian gradient
! find them, that is, in Gx, Gx0, and NAC.

! Note, in the case of nSet==3 we are directly working with GEK predictions of the
! energy difference, the corresponding gradient is  stored in Gx0,
! while for nSet==2 we will have GEK predictions for the individual
! energies and gradients. Hence, for nSet==2 we will put the estimates of the individual
! states into both Gx and Gx0.

If (.NOT.Allocated(BMx_Save)) Call Abend()

!
! dE/dx = dq/dx dE/dq
!

nAtoms = Size(Gx,2)


If (nSet==2)                                          &
call DGEMM_('N','N',                                  &
            3*nAtoms,1,nQQ,                           &
            One,BMx_Save,3*nAtoms,                    &
                Aux(:,iSet-1),nQQ,                    &
           Zero,Gx (:,:,iter+iter_actual-1),3*nAtoms)
call DGEMM_('N','N',                                  &
            3*nAtoms,1,nQQ,                           &
            One,BMx_Save,3*nAtoms,                    &
                Aux(:,iSet),nQQ,                      &
           Zero,Gx0(:,:,iter+iter_actual-1),3*nAtoms)


! Modify with degeneracy factors.

if (Curvilinear) then
   do iAtom=1,nAtoms
      do ixyz=1,3
         If (nSet==2)    &
         Gx (ixyz,iAtom,iter+iter_actual-1) = Gx (ixyz,iAtom,iter+iter_actual-1)/Degen(ixyz,iAtom)
         Gx0(ixyz,iAtom,iter+iter_actual-1) = Gx0(ixyz,iAtom,iter+iter_actual-1)/Degen(ixyz,iAtom)
      end do
   end do
end if
#ifdef _DEBUGPRINT_
Call RecPrt('Kriging_Update: Gx ',' ',Gx (:,:,iter+iter_actual-1),3,nAtoms)
Call RecPrt('Kriging_Update: Gx0',' ',Gx0(:,:,iter+iter_actual-1),3,nAtoms)
#endif



If (nSet==2) Then
   Call mma_deallocate(Aux)
   Return
End If
!                                                                      *
!***********************************************************************
!                                                                      *
!  For the non-adiabatic vector
!
iSet = 3
! The code will assume that the value is Zero by convention.
! ?? = Temp(iSet)

call DGEMM_('N','N',                                  &
            3*nAtoms,1,nQQ,                           &
            One,BMx_Save,3*nAtoms,                    &
                Aux(:,iSet),nQQ,                      &
           Zero,NAC(:,:,iter+iter_actual-1),3*nAtoms)


! Modify with degeneracy factors.

if (Curvilinear) then
   do iAtom=1,nAtoms
      do ixyz=1,3
         NAC(ixyz,iAtom,iter+iter_actual-1) = NAC(ixyz,iAtom,iter+iter_actual-1)/Degen(ixyz,iAtom)
      end do
   end do
end if
#ifdef _DEBUGPRINT_
Call RecPrt('Kriging_Update: NAC',' ',NAC(:,:,iter+iter_actual-1),3,nAtoms)
#endif

! Right now we do not use the dispersion for the constraints

!                                                                      *
!***********************************************************************
!                                                                      *
Call mma_deallocate(Aux)
!
End Subroutine Kriging_Update

