!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

Module DWSol

  use Definitions, only: iwp, u6, wp
  use stdalloc, only: mma_allocate, mma_deallocate

  integer(kind=iwp) :: DWType    = 0_iwp
  integer(kind=iwp) :: IPCMROOT  = 0_iwp
  integer(kind=iwp) :: nRoots    = 0_iwp
  real(kind=wp)     :: DWZeta    = 0.0_wp

  real(kind=wp), allocatable :: W_SOLV(:)

contains
!
!-----------------------------------------------------------------------
!
! dummy call because of an initialization thing
Subroutine DWSol_dummy()
return
End Subroutine DWSol_dummy
!
!-----------------------------------------------------------------------
!
! should be called somewhere in RASSCF or MCLR after nRoots and
! IPCMROOT is determined
Subroutine DWSol_init(IPCMROOT_,nRoots_)

  use rctfld_module, only: iSLPar,rSlPar

  implicit none

  integer(kind=iwp), intent(in) :: IPCMROOT_,nRoots_

  integer(kind=iwp) :: i,istate,jstate
  real(kind=wp)     :: wgt

  IPCMROOT  = IPCMROOT_
  nRoots    = nRoots_
  DWZeta    = RSlPar(54)
  DWType    = ISlPar(17)

  call mma_allocate(W_SOLV,nRoots,Label='W_SOLV')
  W_SOLV(:) = 0.0D+00

  if (DWZeta == 0.0d+00) then
    if (IPCMROOT > 0) then
      W_SOLV(IPCMROOT) = 1.0d+00
    else
      W_SOLV(:) = 1.0d+00/dble(nRoots)
    end if
  else if (DWZeta < 0.0d+00) then
    !! fixed averaging
    Call DWSol_fixed(istate,jstate)
    do i = 1, nRoots
      if (i == istate .or. i == jstate) then
        wgt = 0.5d+00
      else
        wgt = 0.0d+00
      end if
      W_SOLV(i) = wgt
    end do
  end if

End Subroutine DWSol_init
!
!-----------------------------------------------------------------------
!
Subroutine DWSol_final()

  implicit none

  if (allocated(W_SOLV)) call mma_deallocate(W_SOLV)

End Subroutine DWSol_final
!
!-----------------------------------------------------------------------
!
Subroutine DWSol_wgt(ENER)

  implicit none

  real(kind=wp),intent(in) :: ENER(*)

  real(kind=wp)     :: Ealpha,Ebeta,Egamma,wgt,Wtot,xi_ab,xi_ag
  integer(kind=iwp) :: i,istate,jstate

  if (DWZeta <= 0.0d+00) return

  if (DWZeta > 0.0d+00) then
    !! dynamical weighting
    Wtot = 0.0d+00
    Ealpha = ENER(IPCMROOT)
    do i = 1, nRoots
      Egamma = ENER(i)
      xi_ag = (Ealpha-Egamma)**2
      wgt = DWSol_func(xi_ag)
      Wtot = Wtot + wgt
    end do

    do i = 1, nRoots
      Ebeta = ENER(i)
      xi_ab = (Ealpha-Ebeta)**2
      wgt = DWSol_func(xi_ab)
      W_SOLV(i) = wgt/Wtot
    end do
  end if

End Subroutine DWSol_wgt
!
!-----------------------------------------------------------------------
!
Function DWSol_func(xx)

  implicit none

  real(kind=wp) :: DWSol_func
  real(kind=wp), intent(in)  :: xx

  if (DWType == 1) then
    DWSol_func = exp(-DWZeta*xx)
  else
    write (u6,*) 'Unrecognized DWType in DWSol...'
  end if

  return

End Function DWSol_func
!
!-----------------------------------------------------------------------
!
Subroutine DWSol_fixed(istate,jstate)

  implicit none

  integer(kind=iwp), intent(out) :: istate,jstate

  if (DWZeta == -2.0d+00) then
    istate = 1
    jstate = 2
  else if (DWZeta == -6.0d+00) then
    istate = 2
    jstate = 3
  else if (DWZeta == -12.0d+00) then
    istate = 3
    jstate = 4
  else if (DWZeta == -20.0d+00) then
    istate = 4
    jstate = 5
  else if (DWZeta == -30.0d+00) then
    istate = 5
    jstate = 6
  else if (DWZeta == -42.0d+00) then
    istate = 6
    jstate = 7
  else if (DWZeta == -56.0d+00) then
    istate = 7
    jstate = 8
  else if (DWZeta == -72.0d+00) then
    istate = 8
    jstate = 9
  else if (DWZeta == -90.0d+00) then
    istate = 9
    jstate = 10
  else
    write (u6,*) 'Unrecognized negative DWZeta ...'
    istate = 0
    jstate = 0
  end if

  return

End Subroutine DWSol_fixed
!
!-----------------------------------------------------------------------
!
Subroutine DWSol_der(DEROMG,DERHII,ENER)

  implicit none

  real(kind=wp), intent(in) :: DEROMG(:),ENER(:)
  !! partial derivative (pseudo-density, more precisely) of H_{II}
  real(kind=wp), intent(inout) :: DERHII(:)

  integer(kind=iwp) :: i,j,k
  real(kind=wp)     :: Ealpha,Ebeta,Egamma,Scal

  if (DWZeta <= 0.0d+00) Return

  i = IPCMROOT
  Ealpha = ENER(i)

  do j = 1, nRoots
    if (DWType == 1) then
      Ebeta = ENER(j)
      Scal = -2.0d+00*DWZeta*W_SOLV(j)*(Ealpha-Ebeta)*DEROMG(j)
      DERHII(i) = DERHII(i) + Scal
      DERHII(j) = DERHII(j) - Scal

      do k = 1, nRoots
        Egamma = ENER(k)
        Scal = 2.0d+00*DWZeta*W_SOLV(j)*W_SOLV(k)*(Ealpha-Egamma)*DEROMG(j)
        DERHII(i) = DERHII(i) + Scal
        DERHII(k) = DERHII(k) - Scal
      end do
    else
    end if
  end do

  return

End Subroutine DWSol_Der

End Module DWSol
