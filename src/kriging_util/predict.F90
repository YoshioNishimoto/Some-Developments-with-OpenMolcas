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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

subroutine predict(gh)

use kriging_mod, only: cv, full_R, gpred, hpred, Kv, m_t, nInter, ordinary, pred, Rones, sb, sigma, var, variance
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: gh
integer(kind=iwp) :: INFO, i, k
integer(kind=iwp), allocatable :: IPIV(:) ! ipiv the pivot indices that define the permutation matrix
real(kind=wp) :: tsum
real(kind=wp), allocatable :: B(:), A(:,:)

call mma_allocate(B,m_t,label='B')

if (gh == 0) then ! calculate the energy and dispersion

  !A contains the factors L and U from the factorization A = P*L*U as computed by DGETRF
  call mma_allocate(A,m_t,m_t,label='A')
  call mma_allocate(IPIV,m_t,label='IPIV')
  ! calculations of Energy and dispersion
  A(:,:) = full_R(:,:)
  B(:)   = cv(:,1,1)    ! the covariance vector
  pred   = sb(1)+dot_product(B,Kv(:,1))   ! compute the energy

! solve R x = c,  x = R^{-1} c
  call DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO)

  var    = One-dot_product(B,cv(:,1,1))

  if (ordinary) then
    tsum = sum(rones(1:m_t))
    var = max(var+(One-dot_product(cv(:,1,1),rones))**2/tsum,Zero)
  end if

  sigma = sqrt(var*variance(1))
  call mma_deallocate(A)
  call mma_deallocate(IPIV)

else if (gh == 1) then ! calculate the gradient

  do k=1,nInter
    B(:) = cv(:,k,1)
    gpred(k) = dot_product(B,Kv(:,1))
  end do

else if (gh == 2) then ! calculate the Hessian

  do k=1,nInter
    do i=k,nInter
      B(:) = cv(:,i,k)
      hpred(k,i) = dot_product(B,Kv(:,1))
      if (i /= k) hpred(i,k) = hpred(k,i)
    end do
  end do

end if

call mma_deallocate(B)

end subroutine predict
