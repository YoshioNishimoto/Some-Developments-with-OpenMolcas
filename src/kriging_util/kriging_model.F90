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

subroutine kriging_model()

use kriging_mod, only: blAI, blaAI, blavAI, blvAI, detR, dy, full_R, Index_PGEK, Kv, lh, m_t, mblAI, nPoints, nD, nInter_Eff, &
                       ordinary, Rones, sb, sbmev, sbO, variance, y
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, r8, u6

!#define _DPOSV_
!#ifdef _DPOSV_
!use Constants, only: Two
!#endif

implicit none
integer(kind=iwp) :: i_eff, is, ie, ise, iee, i, INFO ! ipiv the pivot indices that define the permutation matrix
integer(kind=iwp), allocatable :: IPIV(:)
real(kind=wp), allocatable :: B(:), A(:,:)
real(kind=r8), external :: dDot_


call mma_Allocate(B,m_t,label='B')            ! the f vector
call mma_Allocate(A,m_t,m_t,label='A')        ! the correlation matrix
call mma_Allocate(IPIV,m_t,label='IPIV')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

B(1:nPoints) = One       ! the f vector, the value part
B(nPoints+1:) = Zero     ! the f vector, the gradient part

A(:,:) = full_r(:,:)     ! the correlation matrix

#ifdef _DEBUGPRINT_
call RecPrt('f',' ',B,1,m_t)
call RecPrt('PSI',' ',A,m_t,m_t)
#endif

!
!  Now form (A^{-1} f) (to be used for the computation of the dispersion)

#ifdef _DPOSV_
call DPOSV_('U',m_t,1,A,m_t,B,m_t,INFO)
#else
call DGESV_(m_t,1,A,m_t,IPIV,B,m_t,INFO)
#endif

if (INFO /= 0) then
  write(u6,*) 'kriging_model: INFO.ne.0'
  write(u6,*) 'kriging_model: INFO=',INFO
  call Abend()
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'Info=',Info
write(u6,*) 'nPoints=',nPoints
call RecPrt('PSI^{-1}',' ',A,m_t,m_t)
call RecPrt('X=PSI^{-1}f',' ',B,1,m_t)
#endif

rones(:) = B(:)     ! Move result over to storage for later use, (R^{-1} f)

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Compute contribution to the expression for the liklihood function.
!
! Now A contains the factors L and U from the factorization A = P*L*U as computed by DGESV
! Where L in the lower triangular matrix with 1 in the diagonal and U is the upper
! triangular matrix of A thus the determinant of A is giving by multipling its diagonal

detR = Zero
do i=1,m_t
#ifdef _DPOSV_
  detR = detR+Two*log(A(i,i))
#else
  detR = detR+log(abs(A(i,i)))
#endif
end do

!
!  Now work on the trend function and the vector with the values and gradients, the generalized y vector
!
!  sb   : mu, the trend function
!  rones: (R^{-1} f)
!  kv   :  R^{-1} (y - mu f)

!**********************************************************************
!
! 1) Establish the trend function (baseline), mu.
!
!**********************************************************************

if (blaAI) then    ! This is the default.

  ! Make sure the base line is above any data point, this to make sure the surrogate model is bound.

  sb = -huge(sb)
  do i=1,nPoints
    sb = max(sb,y(i,1)+blavAI)
  end do

else if (mblAI) then

  sb = sbmev

else if (blAI) then

  sb = blvAI

else

  ! Derived according to ordinary kriging

  ordinary = .true.

  ! Note dy only containes gradients for the set of points which we are considering
  ! B(:) = [y(:),dy(:)]  original code before PGEK implementation

  B(1:nPoints) = y(1:nPoints,1)
  do i_eff=1,nInter_eff
    i = Index_PGEK(i_eff)
    is = (i-1)*(nPoints-nD)+1
    ise = nPoints+(i_eff-1)*(nPoints-nD)+1
    ie = is+(nPoints-nD)-1
    iee = ise+(nPoints-nD)-1
    B(ise:iee) = dy(is:ie,1)
  end do

#ifdef _DEBUGPRINT_
  write(u6,*) DDot_(m_t,rones,1,B,1),DDot_(nPoints,rones,1,[One],0)
#endif
  ! mu =  (f R^{-1} y) /(f R^{-1} f)
  sbO = DDot_(m_t,rones,1,B,1)/DDot_(nPoints,rones,1,[One],0)

  sb = sbO
end if

!**********************************************************************
!
! 2) form the value vector, (y - mu f)
!
!**********************************************************************

!     B(1:m_t) = [y(1:nPoints)-sb,dy(1:nInter*(nPoints-nD))]

! the values
B(1:nPoints) = y(1:nPoints,1)-sb
! the gradients
do i_eff=1,nInter_eff
  i = Index_PGEK(i_eff)
  is = +(i-1)*(nPoints-nD)+1
  ise = nPoints+(i_eff-1)*(nPoints-nD)+1
  ie = is+(nPoints-nD)-1
  iee = ise+(nPoints-nD)-1
  B(ise:iee) = dy(is:ie,1)
end do

#ifdef _DEBUGPRINT_
write(u6,*) 'sb,ln(det|PSI|)=',sb,detR
call RecPrt('[y-sb,dy]','(12(2x,E9.3))',B,1,m_t)
#endif

Kv(:) = B(:)             ! The value vector

#ifdef _DEBUGPRINT_
call RecPrt('Kv',' ',Kv,1,m_t)
#endif

A(:,:) = full_R(:,:)

!
! Solve R x = (y - mu f), i.e. x = R^{-1} (y - mu f)
!

#ifdef _DPOSV_
call DPOSV_('U',m_t,1,A,m_t,Kv,m_t,INFO)
#else
call DGESV_(m_t,1,A,m_t,IPIV,Kv,m_t,INFO)
#endif

!   Compute the dispersion
!
!   s^2 = (y - mu f) R^{-1} (y - mu f) / n

variance = dot_product(B,Kv)/real(m_t,kind=wp)

! compute the value of the likelihood function

lh = variance*exp(detR/real(m_t,kind=wp))

#ifdef _DEBUGPRINT_
write(u6,*) 'Variance=',Variance
write(u6,*) 'Info=',Info
write(u6,*) 'lh=',lh
call RecPrt('X=A^{-1}Kv','(5(E15.7,2X))',Kv,1,m_t)
#endif

call mma_Deallocate(B)
call mma_Deallocate(A)
call mma_Deallocate(IPIV)

end subroutine kriging_model
