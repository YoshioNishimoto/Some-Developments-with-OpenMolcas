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
Use Slapaf_Info, only: Energy, dqInt
Implicit None
Integer nQQ, iter
Real*8  qInt(nQQ), E_Disp
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
!  The default case
!
Call Energy_Kriging_layer(qInt(:),Energy(iter),nQQ)

Call Dispersion_Kriging_Layer(qInt(:),E_Disp,nQQ)

Call Gradient_Kriging_layer(qInt(:),dqInt(:,iter),nQQ)

dqInt(:,iter) = - dqInt(:,iter)
!                                                                      *
!***********************************************************************
!                                                                      *
!  For the energy difference
!                                                                      *
!***********************************************************************
!                                                                      *
!  For the non-adiabatic vector
!                                                                      *
!***********************************************************************
!                                                                      *
End Subroutine Kriging_Update

