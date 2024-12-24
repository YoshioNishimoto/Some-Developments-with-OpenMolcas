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
! Copyright (C) 2024, Yoshio Nishimoto                                 *
!***********************************************************************
!-----------------------------------------------------------------------
!
  Subroutine CGS_pre(nDensC,nConf1,nRoots,Kappa,dKappa,ipST,ipCId,delta)

  use definitions, only: iwp,wp
  use ipPage, only: W
  use ISRotation, only: ISR,InvSCF
  use cgs_mod, only: CGSvec

  implicit none

  integer(kind=iwp), intent(in) :: nDensC,nConf1,nRoots,ipST,ipCId
  real(kind=wp), intent(inout) :: Kappa(*),dKappa(*)
  real(kind=wp), intent(out) :: delta

  real(kind=wp) :: deltaC,deltaK
  real(kind=wp), external :: ddot_

  ! R0
  call dcopy_(nDensC,Kappa,1,CGSvec%R0,1)
  call dcopy_(nConf1*nRoots,W(ipST)%Vec,1,W(CGSvec%ipR0)%Vec,1)
  if (.not.InvSCF) ISR%R0 = +ISR%Rvec

  ! u_k
  call dcopy_(nDensC,Kappa,1,CGSvec%Uvec,1)
  call dcopy_(nConf1*nRoots,W(ipST)%Vec,1,W(CGSvec%ipUvec)%Vec,1)
  if (.not.InvSCF) ISR%Uvec = +ISR%Rvec

  ! p_k
  call dcopy_(nDensC,CGSvec%Uvec,1,dKappa,1)
  Call DCopy_(nDensC,dKappa,1,CGSvec%Pvec,1)
  call dcopy_(nConf1*nRoots,W(CGSvec%ipUvec)%Vec,1,W(ipCId)%Vec,1)
  Call DCopy_(nconf1*nRoots,W(ipCId)%Vec,1,W(CGSvec%ipPvec)%Vec,1)
  if (.not.InvSCF) then
    ISR%Pvec = ISR%Uvec
    ISR%Rvec = ISR%R0
    ISR%Xvec = 0.0d+00
  end if

  deltaC = ddot_(nConf1*nroots,W(CGSvec%ipR0)%Vec,1,W(CGSvec%ipR0)%Vec,1)
  if (.not.InvSCF) deltaC = deltaC + ddot_(nRoots**2,ISR%R0,1,ISR%R0,1)
  deltaK = ddot_(nDensC,CGSvec%R0,1,CGSvec%R0,1)
  Kappa(1:nDensC)=0.0d+00
  delta  = deltaC + deltaK

  End Subroutine CGS_pre
