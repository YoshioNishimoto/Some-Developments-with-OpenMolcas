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
!
! OpenMolcas wrapper for the soMPSoo program
! sompsoo_init initializes some settings for soMPSoo and OpenMolcas
!
!*******************

subroutine sompsoo_init(refwfnfile)

! set global variables directly from the soMPSoo program
use global_control  

use Definitions, only: iwp, u6

implicit none

!> Cholesky stuff
#include "chlcas.fh"
#include "chodensity.fh"
#include "chotime.fh"
#include "cholk.fh"
#include "choscreen.fh"
#include "chopar.fh"


character(len=*), intent(inout) :: refwfnfile
integer(kind=iwp)               :: i, j, m, n


! Initial values
ifcore           = .false.
do_soMPSoo_srdft = .false.
sopri            = u6;

! If (default) further SA-DMRG-SCF is needed for
!     Coupled-Perturbation/Lagrange calculations
DMRG_SCF= .true.
! include active-active rotations
aaRotations = .false.

!> number of point-group irreps
orb%nsub = -1

!> number of states to optimize
dmrg_nstates = 1

!> always allow a coupling between orbitals and MPS in WMK solver
CP_integrals=.true.
ith_INTE = 2

!> default NL solver is AH
do i=1,MAX_macro_iterations
  method(i)="AH        "
end do

! initialize optimization parameters
thrs%r         = 1.0d-8   ! rotation
thrs%s         = 1.0d-8   ! symmetry
thrs%e         = 1.0d-8   ! energy
thrs%d         = 1.0d-12  ! davidson
thrs%c         = 1.0d-4   ! coupling
Nmps_update    = 10
Nint_update    = 0
!> max micro-iterations
Max_iter             = 10
!> max macro-iterations
num_macro_iterations = -1


!> OpenMolcas initial values
refwfnfile = 'JOBIPH'

!> set default values for orbitals
!>
allocate(orb%frozen(max_nsym_openm)); orb%frozen=0 !define
allocate(orb%closed(max_nsym_openm)); orb%closed=0
allocate(orb%   occ(max_nsym_openm)); orb%   occ=0
allocate(orb%   act(max_nsym_openm)); orb%   act=0
allocate(orb% extnl(max_nsym_openm)); orb% extnl=0
allocate(orb% total(max_nsym_openm)); orb% total=0

allocate(orbbk%frozen(max_nsym_openm)); orbbk%frozen=0 !define
allocate(orbbk%closed(max_nsym_openm)); orbbk%closed=0
allocate(orbbk%   occ(max_nsym_openm)); orbbk%   occ=0
allocate(orbbk%   act(max_nsym_openm)); orbbk%   act=0
allocate(orbbk% extnl(max_nsym_openm)); orbbk% extnl=0
allocate(orbbk% total(max_nsym_openm)); orbbk% total=0

allocate(orb%grouptable(max_nsym_openm,max_nsym_openm));orb%grouptable=0
allocate(orb%invelm(max_nsym_openm)); orb%invelm=0
allocate(orb%adjsym(max_nsym_openm)); orb%adjsym=0

allocate(off%iorb(max_nsym_openm)); off%iorb=0
allocate(off%icmo(max_nsym_openm)); off%icmo=0
allocate(off%ibas(max_nsym_openm)); off%ibas=0

allocate(bas%nbas(max_nsym_openm)); bas%nbas=0

!> for consistency within OpenMolcas we use the same initialization as in rasscf
!> SET UP SYMMETRY MULTIPLICATION TABLE:
      orb%grouptable(1,1)=1
      M=1
      DO  N=1,3
        DO  I=1,M
          DO  J=1,M
            orb%grouptable(I+M,J)   = M + orb%grouptable(I,J)
            orb%grouptable(I,J+M)   =     orb%grouptable(I+M,J)
            orb%grouptable(I+M,J+M) =     orb%grouptable(I,J)
          END DO
         END DO
        M=2*M
      END DO

!> Cholesky-related settings:
      Call DecideOnCholesky(DoCholesky)
      ALGO  = 1
      DensityCheck=.false.
      Deco=.true.
      timings=.false.
      DoLock=.true.
      Nscreen=10
      dmpk=1.0d-1
      Update=.true.
      Estimate=.false.

#if defined (_MOLCAS_MPP_)
      ChFracMem=0.3d0
#else
      ChFracMem=0.0d0
#endif

end subroutine sompsoo_init
