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

  Type type_DW
    logical(kind=iwp) :: do_DW   = .false.
    integer(kind=iwp) :: DWType  = -99
    integer(kind=iwp) :: DWRoot  = 0
    real(kind=wp)     :: DWZeta  = 0.0d+00
  End Type type_DW

  Type (type_DW), target :: DWSolv, DWSCF

  !! For dynamically weighted solvation
! logical(kind=iwp) :: DWSolv    = .false.
! integer(kind=iwp) :: DWType    = 0_iwp
! real(kind=wp)     :: DWZeta    = 0.0_wp
  integer(kind=iwp) :: IPCMROOT  = 0_iwp

  !! For dynamically weighted MCSCF
! logical(kind=iwp) :: DWSCF     = .false.
! integer(kind=iwp) :: DWTypeSCF = 0_iwp   ! DWTY in &RASSCF
! real(kind=wp)     :: DWZetaSCF = 0.0_wp  ! DWZE in &RASSCF
! integer(kind=iwp) :: IDWROOT   = 0_iwp   ! DWRO in &RASSCF

  integer(kind=iwp) :: nRoots    = 0_iwp


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
Subroutine DWSol_init(IPCMROOT_,nRoots_,NonEq_)

  use rctfld_module, only: iSLPar,rSlPar

  implicit none

  integer(kind=iwp), intent(in) :: IPCMROOT_,nRoots_
  logical(kind=iwp), intent(in) :: NonEq_

  integer(kind=iwp) :: i,istate,jstate
  real(kind=wp)     :: wgt

  IPCMROOT      = IPCMROOT_
  nRoots        = nRoots_

  DWSolv%DWType = ISlPar(17)
  DWSolv%DWRoot = IPCMROOT
  DWSolv%DWZeta = RSlPar(54)

  if (NonEq_) DWSolv%DWZeta = 0.0d+00

  call mma_allocate(W_SOLV,nRoots,Label='W_SOLV')
  W_SOLV(:) = 0.0D+00

  if (DWSolv%DWZeta == 0.0d+00) then
    if (IPCMROOT > 0) then
      W_SOLV(IPCMROOT) = 1.0d+00
    else
      W_SOLV(:) = 1.0d+00/dble(nRoots)
    end if
  else if (DWSolv%DWZeta < 0.0d+00) then
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
  else
    DWSolv%do_DW = .true.
  end if

End Subroutine DWSol_init
!
!-----------------------------------------------------------------------
!
! should be called somewhere in RASSCF or MCLR after nRoots and
! IPCMROOT is determined
Subroutine DWSCF_init(mode,nRoots_)

  implicit none

  integer(kind=iwp), intent(in) :: mode,nRoots_
  integer(kind=iwp) :: idum
  real(kind=wp) :: rdum

  nRoots    = nRoots_

  if (mode==1) then
    !! call from RASSCF, nothing is done
  else if (mode==2) then
    !! call from MCLR, read from RunFile(?)
    Call Get_iScalar('DWTypeSCF',idum)
    Call Get_dScalar('DWZetaSCF',rdum)
    DWSCF%DWRoot = ibits(idum,16,16)
    DWSCF%DWType = ibits(idum, 0,16)
    DWSCF%DWZeta = rdum
!   write (u6,*) "parameters have been read in DWSCF_init"
!   write (u6,*) DWSCF%DWRoot,DWSCF%DWtype,DWSCF%DWZeta
  end if

  DWSCF%do_DW = .false.
  if (DWSCF%DWROOT /= 0) then
    write (u6,*) 'DW-MCSCF is enabled'
    DWSCF%do_DW = .true.
    if (DWSCF%DWType == -99) then
      write (u6,*) 'DWTYpe has not been specified.'
      write (u6,*) 'The default value DWTYpe = 0 will be used.'
      DWSCF%DWType = 0
    end if
    if (DWSCF%DWZeta == 0.0d+00) then
      if (DWSCF%DWType ==  -1) DWSCF%DWZeta = 1.0d+00
      if (DWSCF%DWType ==   0) DWSCF%DWZeta = 1.0d+00
      if (DWSCF%DWType ==   1) DWSCF%DWZeta = 1.0d+00
    end if
  end if

End Subroutine DWSCF_init
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
Subroutine DWSCF_final()

  implicit none

  integer(kind=iwp) :: idum
  real(kind=wp) :: rdum

  if (DWSCF%do_DW) then
!   idum = ishft(DWSCF%DWRoot,8) + DWSCF%DWType
    idum = DWSCF%DWType
    call mvbits(DWSCF%DWRoot,0,16,idum,16)
    rdum = DWSCF%DWZeta
!   write (*,*) dwscf%dwroot, dwscf%dwtype
!   write (*,*) "idum = ", idum
!   write (*,*) "idum = ", ibits(DWSCF%DWRoot,-8,32) + ibits(DWSCF%DWType,0,8)
!   write (*,*) ibits(idum,0,8)
!   write (*,*) ibits(idum,8,8)
  else
    idum = 0
    rdum = 0.0d+00
  end if

  Call Put_iScalar('DWTypeSCF',idum)
  Call Put_dScalar('DWZetaSCF',rdum)

End Subroutine DWSCF_final
!
!-----------------------------------------------------------------------
!
Subroutine DWSol_wgt(mode,ENER,weight)

  implicit none

  integer(kind=iwp),intent(in) :: mode
  real(kind=wp),intent(in) :: ENER(*)
  real(kind=wp),intent(out),optional :: weight(*)

  real(kind=wp), allocatable :: W_local(:)
  real(kind=wp)     :: Ealpha,Ebeta,Egamma,wgt,Wtot,xi_ab,xi_ag
  integer(kind=iwp) :: i,istate,jstate

  Type (type_DW), pointer :: DWlocal

  !! DW-MCSCF
  if (mode==1) DWlocal => DWSCF
  !! DW solvation
  if (mode==2) DWlocal => DWSolv
  !! DW CASPT2?
! if (mode==3)

  if (DWlocal%DWZeta <= 0.0d+00) return

! write (*,*) dwlocal%do_dw
! write (*,*) dwlocal%dwtype
! write (*,*) dwlocal%dwroot
! write (*,*) dwlocal%dwzeta
! write (*,*) ener(1:nroots)

  !! dynamical weighting
  call mma_allocate(W_local,nRoots,Label='W_local')

  !! Normalization
  Wtot = 0.0d+00
  Ealpha = ENER(DWlocal%DWRoot)
  do i = 1, nRoots
    Egamma = ENER(i)
    xi_ag = Ealpha-Egamma
    wgt = DWSol_func(xi_ag,DWlocal)
    Wtot = Wtot + wgt
  end do

  W_local = 0.0d+00
  do i = 1, nRoots
    Ebeta = ENER(i)
    xi_ab = Ealpha-Ebeta
    wgt = DWSol_func(xi_ab,DWlocal)
    W_local(i) = wgt/Wtot
  end do

  if (mode==1) weight(1:nRoots) = W_local(1:nRoots)
  if (mode==2) W_SOLV(1:nRoots) = W_local(1:nRoots)
! write (*,*) "nroots = ", nroots
! write (*,*) "dw weight", weight(1:nroots)

  call mma_deallocate(W_local)

End Subroutine DWSol_wgt
!
!-----------------------------------------------------------------------
!
Function DWSol_func(xx,DWlocal)

  implicit none

  real(kind=wp) :: DWSol_func
  real(kind=wp), intent(in) :: xx
  type(type_DW), intent(in) :: DWlocal

  if (DWlocal%DWType == -1) then
    !! the oldest (?) DW-CASSCF (J. Chem. Phys. 120, 7281 (2004)
    DWSol_func = cosh(xx/DWlocal%DWZeta)**2
  else if (DWlocal%DWType == 0) then
    !! new (?) DW-CASSCF: J. Chem. Phys. 141, 171102 (2014)
    if (xx.le.0.0d+00) then
      DWSol_func = 1.0d+00/nRoots
    else if (xx.le.DWlocal%DWZeta) then
      DWSol_func = 1.0d+00 &
                 - 3.0d+00*(xx/DWlocal%DWZeta)**2 &
                 + 2.0d+00*(xx/DWlocal%DWZeta)**3
    else
      DWSol_func = 0.0d+00
    end if
  else if (DWlocal%DWType == 1) then
    DWSol_func = exp(-DWlocal%DWZeta*xx*xx)
!   write (*,*) "dwsol_func = ", dwsol_func,dwlocal%dwzeta
  else
    write (u6,*) 'Unrecognized DWType in DWSol...'
    call abend()
  end if

  return

End Function DWSol_func
!
!-----------------------------------------------------------------------
!
Subroutine DWSol_fixed(istate,jstate)

  implicit none

  integer(kind=iwp), intent(out) :: istate,jstate

  if (DWSolv%DWZeta == -2.0d+00) then
    istate = 1
    jstate = 2
  else if (DWSolv%DWZeta == -6.0d+00) then
    istate = 2
    jstate = 3
  else if (DWSolv%DWZeta == -12.0d+00) then
    istate = 3
    jstate = 4
  else if (DWSolv%DWZeta == -20.0d+00) then
    istate = 4
    jstate = 5
  else if (DWSolv%DWZeta == -30.0d+00) then
    istate = 5
    jstate = 6
  else if (DWSolv%DWZeta == -42.0d+00) then
    istate = 6
    jstate = 7
  else if (DWSolv%DWZeta == -56.0d+00) then
    istate = 7
    jstate = 8
  else if (DWSolv%DWZeta == -72.0d+00) then
    istate = 8
    jstate = 9
  else if (DWSolv%DWZeta == -90.0d+00) then
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
Subroutine DWSol_der(mode,DEROMG,DERHII,ENER,weight)

  implicit none

  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: DEROMG(:),ENER(:)
  real(kind=wp), intent(in), optional :: weight(1:nRoots)
  !! partial derivative (pseudo-density, more precisely) of H_{II}
  real(kind=wp), intent(inout) :: DERHII(:)

  integer(kind=iwp) :: i,j,k
  real(kind=wp)     :: Ealpha,Ebeta,Egamma,Scal
  real(kind=wp), allocatable :: W_local(:)
  Type (type_DW), pointer :: DWlocal

  !! DW-MCSCF
  if (mode==1) DWlocal => DWSCF
  !! DW solvation
  if (mode==2) DWlocal => DWSolv

  DERHII(1:nRoots) = 0.0d+00
  if (DWlocal%DWZeta <= 0.0d+00) return

  call mma_allocate(W_local,nRoots,Label='W_local')

  !! DW-MCSCF
  if (mode==1) W_local(1:nRoots) = weight(1:nRoots)
  !! DW solvation
  if (mode==2) W_local(1:nRoots) = W_SOLV(1:nRoots)

! i = IPCMROOT
  i = DWlocal%DWRoot
  Ealpha = ENER(i)

  do j = 1, nRoots
    if (DWlocal%DWType == 1) then
      Ebeta = ENER(j)
      Scal = -2.0d+00*DWlocal%DWZeta*W_local(j)*(Ealpha-Ebeta)*DEROMG(j)
      DERHII(i) = DERHII(i) - Scal
      DERHII(j) = DERHII(j) + Scal

      do k = 1, nRoots
        Egamma = ENER(k)
        Scal = 2.0d+00*DWlocal%DWZeta*W_local(j)*W_local(k)*(Ealpha-Egamma)*DEROMG(j)
        DERHII(i) = DERHII(i) - Scal
        DERHII(k) = DERHII(k) + Scal
      end do
    else
    end if
  end do

! write (*,*) "derhii = ", derhii(1:nroots)
! derhii(1:nroots) = 0.0d+00

! do j = 1, nRoots
!   if (DWlocal%DWType == 1) then
!     Ebeta = ENER(j)
!     Scal = -2.0d+00*DWlocal%DWZeta*W_local(j)*(Ealpha-Ebeta)*DEROMG(j)
!     DERHII(i) = DERHII(i) - Scal
!     DERHII(j) = DERHII(j) + Scal

!     do k = 1, nRoots
!       Egamma = ENER(k)
!       Scal = 2.0d+00*DWlocal%DWZeta*W_local(j)*W_local(k)*(Ealpha-Egamma)*DEROMG(j)
!       DERHII(i) = DERHII(i) - Scal
!       DERHII(k) = DERHII(k) + Scal
!     end do
!   else
!   end if
! end do

! write (*,*) "derhii = ", derhii(1:nroots)

  call mma_deallocate(W_local)

  return

End Subroutine DWSol_Der

End Module DWSol
