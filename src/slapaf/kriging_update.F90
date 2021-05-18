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
Use Slapaf_Info, only: Energy, dqInt, Energy0, BMx_Save
Use Kriging_Mod, only: nSet, iter_actual
Implicit None
Integer nQQ, iter
Real*8  qInt(nQQ), E_Disp

#include "stdalloc.fh"
Integer iSet
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

Call Energy_Kriging_layer(qInt(:),Temp,nQQ)

Call Dispersion_Kriging_Layer(qInt(:),Demp,nQQ)

Call Gradient_Kriging_layer(qInt(:),Aux,nQQ)

!                                                                      *
!***********************************************************************
!
!  The default case
!
iSet = 1
Energy(iter) = Temp(iSet)

E_Disp = Demp(iSet)

dqInt(:,iter) = - Aux(:,iSet)

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
Energy0(iter_actual+1) = Temp(iSet)

! Right now we do not use the dispersion for the constraints

! We should not have the gradient of the energy difference with respect
! to the internal coordinates. We now need to transform those to
! Cartesians and store them at the correct place, that is, in Gx0.

If (.NOT.Allocated(BMx_Save)) Call Abend()

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

! Right now we do not use the dispersion for the constraints

!                                                                      *
!***********************************************************************
!                                                                      *
Call mma_deallocate(Aux)
!
End Subroutine Kriging_Update

