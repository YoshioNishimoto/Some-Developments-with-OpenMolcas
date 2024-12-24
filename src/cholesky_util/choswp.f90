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
Module ChoSwp
Implicit none
Private
Public:: iQuAB, iQuAB_L, iQuAB_Hidden, iQuAB_L_Hidden, pTemp, iQuAB_here, &
         nnBstRSh, nnBstRSh_Hidden, nnBstRSh_G, nnBstRSh_L_Hidden, pTemp3, &
         iiBstRSh, iiBstRSh_Hidden, iiBstRSh_G, iiBstRSh_L_Hidden, &
         IndRSh, IndRSh_Hidden, IndRSh_G, IndRSh_G_Hidden, pTemp1, &
         InfRed, InfRed_Hidden, InfRed_G, InfRed_G_Hidden, &
         InfVec, InfVec_Hidden, InfVec_G, InfVec_G_Hidden, InfVec_Bak, &
         IndRed, IndRed_Hidden, IndRed_G, IndRed_G_Hidden, &
         Diag, Diag_Hidden, Diag_G, Diag_G_Hidden


Integer, Allocatable, Target:: iQuAB_Hidden(:,:), iQuAB_L_Hidden(:,:), iQuAB_here(:,:)
Integer, Pointer:: iQuAB(:,:)=>Null() , iQuAB_L(:,:)=>Null(), pTemp(:,:)=>Null()

Integer, Allocatable, Target:: nnBstRSh_Hidden(:,:,:), nnBstRSh_L_Hidden(:,:,:)
Integer, Allocatable, Target:: iiBstRSh_Hidden(:,:,:), iiBstRSh_L_Hidden(:,:,:)
Integer, Pointer:: nnBstRSh(:,:,:)=>Null(), nnBstRSh_G(:,:,:)=>Null(), pTemp3(:,:,:)=>Null()
Integer, Pointer:: iiBstRSh(:,:,:)=>Null(), iiBstRSh_G(:,:,:)=>Null()
Integer, Pointer:: IndRSh(:)=>Null(), IndRSh_G(:)=>Null(), pTemp1(:)=>Null()
Integer, Allocatable, Target:: IndRSh_Hidden(:), IndRSh_G_Hidden(:)
Integer, Pointer:: InfRed(:)=>Null(), InfRed_G(:)=>Null()
Integer, Allocatable, Target:: InfRed_Hidden(:), InfRed_G_Hidden(:)
Integer, Pointer:: InfVec(:,:,:)=>Null(), InfVec_G(:,:,:)=>Null()
Integer, Allocatable, Target:: InfVec_Hidden(:,:,:), InfVec_G_Hidden(:,:,:)
Integer, Allocatable, Target:: InfVec_Bak(:,:,:)
Integer, Pointer:: IndRed(:,:)=>Null(), IndRed_G(:,:)=>Null()
Integer, Allocatable, Target:: IndRed_Hidden(:,:), IndRed_G_Hidden(:,:)
Real*8, Pointer:: Diag(:)=>Null(), Diag_G(:)=>Null()
Real*8, Allocatable, Target:: Diag_Hidden(:), Diag_G_Hidden(:)
End Module ChoSwp
