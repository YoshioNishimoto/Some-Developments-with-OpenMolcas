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
!
! See cgs_mod.f90 for details
!
  Subroutine CGS_x(nDensC,nDens2,nConf1,nRoots,iSym,jspin,iter, &
                   ipCI,ipS1,ipS2,ipST,ipCIT,ipCId,ipDia,ipPre,ipPre2, &
                   reco,Fancy,Kappa,dKappa,Sigma,Temp4,Sc1,delta,resk,resci,deltak,deltac)

  use definitions, only: iwp,wp
  use ipPage, only: W
  use ISRotation, only: ISR,InvSCF,DMInvISR
  use cgs_mod, only: CGSvec

  implicit none

  integer(kind=iwp), intent(in) :: nDensC,nDens2,nConf1,nRoots,iSym,jspin
  integer(kind=iwp), intent(inout) :: iter
  integer(kind=iwp), intent(in) :: ipCI,ipS1,ipS2,ipST,ipCIT,ipCId,ipDia,ipPre,ipPre2
  real(kind=wp), intent(inout) :: reco,Fancy(*),Kappa(*),dKappa(*),Sigma(*),Temp4(*),Sc1(*), &
                                  delta,resk,resci,deltak,deltac

  integer(kind=iwp) :: irc
  integer(kind=iwp),external :: ipIn,opOut
  real(kind=wp) :: rdum(1),ralpha,rAlphaK,rAlphaC,rbeta
  real(kind=wp), external :: ddot_

  !! precondition p (K-1*p)
  irc=ipIn(ipS2)
  Call DMinvCI_SA(CGSvec%ipPvec,W(ipCId)%Vec,rdum(1),isym,Fancy)
  irc=opOut(ipci)
  irc=opOut(ipdia)

  irc=ipIn(ipPre2)
  Call DMInvKap(W(ipPre2)%Vec,CGSvec%Pvec,nDens2+6,dKappa,nDens2+6,Sc1,nDens2+6,iSym,iter)
  irc=opOut(ippre2)

  if (.not.InvSCF) Call DMInvISR(ISR%Pvec,ISR%p)

  !! A*K-1*p
  Call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

  !! compute alpha
  rAlphaK=ddot_(nDensC,CGSvec%R0,1,Temp4,1)
  rAlphaC=ddot_(nConf1*nroots,W(CGSvec%ipR0)%Vec,1,W(ipS1)%Vec,1)
  if (.not.InvSCF) rAlphaC = rAlphaC + ddot_(nRoots**2,ISR%R0,1,ISR%Ap,1)
  rAlpha=delta/(rAlphaK+rAlphaC)

  ! q_k = u_k - alpha*(A*K-1*p)
  call dcopy_(nDensC,CGSvec%Uvec,1,CGSvec%Qvec,1)
  call daxpy_(nDensC,-ralpha,Temp4,1,CGSvec%Qvec,1)
  call dcopy_(nConf1*nRoots,W(CGSvec%ipUvec)%Vec,1,W(CGSvec%ipQvec)%Vec,1)
  call daxpy_(nConf1*nRoots,-ralpha,W(ipS1)%Vec,1,W(CGSvec%ipQvec)%Vec,1)
  if (.not.InvSCF) ISR%Qvec = ISR%Uvec - ralpha*ISR%Ap

  !! precondition u + q
  call daxpy_(nDensC,+1.0d+00,CGSvec%Qvec,1,CGSvec%Uvec,1)
  call daxpy_(nConf1*nRoots,+1.0d+00,W(CGSvec%ipQvec)%Vec,1,W(CGSvec%ipUvec)%Vec,1)
  if (.not.InvSCF) ISR%Uvec = ISR%Uvec + ISR%Qvec

  irc=ipIn(ipS2)
  Call DMinvCI_SA(CGSvec%ipUvec,W(ipCId)%Vec,rdum(1),isym,Fancy)
  irc=opOut(ipci)
  irc=opOut(ipdia)

  irc=ipIn(ipPre2)
  Call DMInvKap(W(ipPre2)%Vec,CGSvec%Uvec,nDens2+6,dKappa,nDens2+6,Sc1,nDens2+6,iSym,iter)
  irc=opOut(ippre2)

  if (.not.InvSCF) Call DMInvISR(ISR%Uvec,ISR%p)

  !! x_k+1 = x_k + alpha*K-1*(u+q)
  call daxpy_(nDensC,+ralpha,dKappa,1,Kappa,1)
  Call DaXpY_(nConf1*nroots,ralpha,W(ipCId)%Vec,1,W(ipCIT)%Vec,1)
  if (.not.InvSCF) ISR%Xvec = ISR%Xvec + ralpha*ISR%p

  !! A*K-1*(u+q)
  Call TimesE2(dKappa,ipCId,1,reco,jspin,ipS2,Temp4,ipS1)

  !! r_k+1 = r_k - alpha*(A*K-1*(u+q))
  Call DaxPy_(nDensC,-ralpha,Temp4,1,Sigma,1)
  resk=sqrt(ddot_(nDensC,Sigma,1,Sigma,1))
  Call DaXpY_(nConf1*nroots,-ralpha,W(ipS1)%Vec,1,W(ipST)%Vec,1)
  resci=sqrt(ddot_(nconf1*nroots,W(ipST)%Vec,1,W(ipST)%Vec,1))
  if (.not.InvSCF) then
    ISR%Rvec = ISR%Rvec - ralpha*ISR%Ap
    resci = resci + sqrt(ddot_(nRoots**2,ISR%Rvec,1,ISR%Rvec,1))
  end if

  !! compute beta
  deltaK=ddot_(nDensC,CGSvec%R0,1,Sigma,1)
  deltaC=ddot_(nConf1*nroots,W(CGSvec%ipR0)%Vec,1,W(ipST)%Vec,1)
  if (.not.InvSCF) deltaC = deltaC + ddot_(nRoots**2,ISR%R0,1,ISR%Rvec,1)
  rbeta=(deltaC+deltaK)/delta
  delta=deltac+deltaK

  !! u_k = r_k + beta*q
  call dcopy_(nDensC,Sigma,1,CGSvec%Uvec,1)
  call daxpy_(nDensC,+rbeta,CGSvec%Qvec,1,CGSvec%Uvec,1)
  call dcopy_(nConf1*nRoots,W(ipST)%Vec,1,W(CGSvec%ipUvec)%Vec,1)
  Call DaXpY_(nConf1*nroots,+rbeta,W(CGSvec%ipQvec)%Vec,1,W(CGSvec%ipUvec)%Vec,1)
  if (.not.InvSCF) ISR%Uvec = ISR%Rvec + rbeta*ISR%Qvec

  !! p = u + beta*q + beta*beta*p
  call dscal_(nDensC,rbeta**2,CGSvec%Pvec,1)
  call daxpy_(nDensC,+1.0d+00,CGSvec%Uvec,1,CGSvec%Pvec,1)
  call daxpy_(nDensC,+rbeta,CGSvec%Qvec,1,CGSvec%Pvec,1)
  call dscal_(nConf1*nRoots,rbeta**2,W(CGSvec%ipPvec)%Vec,1)
  Call DaXpY_(nConf1*nroots,+1.0d+00,W(CGSvec%ipUvec)%Vec,1,W(CGSvec%ipPvec)%Vec,1)
  Call DaXpY_(nConf1*nroots,+rbeta,W(CGSvec%ipQvec)%Vec,1,W(CGSvec%ipPvec)%Vec,1)
  if (.not.InvSCF) ISR%Pvec = ISR%Uvec + rbeta*ISR%Qvec + rbeta*rbeta*ISR%Pvec

  End Subroutine CGS_x
