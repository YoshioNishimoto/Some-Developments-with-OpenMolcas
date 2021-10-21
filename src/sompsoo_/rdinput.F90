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
! rdinput reads input from the OpenMolcas input file
!
!*******************

subroutine rdinput(refwfnfile)

! set global variables directly from the soMPSoo program
use global_control  

use Definitions, only: iwp, u6

implicit none
character(len=*), intent(inout) :: refwfnfile
character(len=180) :: Line, key
integer(kind=iwp) :: LuSpool, iError
integer(kind=iwp), external :: isFreeUnit
logical(kind=iwp), external :: next_non_comment
character(len=180), external :: Get_Ln
integer(kind=iwp) :: total_weight, nr_of_wf, i, j
integer(kind=iwp) :: istep, estep, k, l, ierr, nstep
character(len=10) :: local_method


istep = 1; estep = 0; i = 0; j=0; k = 0; l = 0;

LuSpool = isFreeUnit(18)
call SpoolInp(LuSpool)

rewind(LuSpool)
call RdNLst(LuSpool,'soMPSoo')

!> scan for # of WF types
nr_of_wf = 0
do
  key = Get_Ln(LuSpool)
  call LeftAd(key)
  Line = key
  if (Line(1:1) == '*') cycle
  if (Line == ' ') cycle
  call UpCase(Line)
  if(Line(1:4) == 'END ') exit
  if(Line(1:2).eq."WF") nr_of_wf = nr_of_wf + 1;
end do

allocate(wf(nr_of_wf),stat=iError); if(iError /= 0) call error(2)

!> position is now again at &soMPSoo
rewind(LuSpool)
call RdNLst(LuSpool,'soMPSoo')

j = 0

do
  key = Get_Ln(LuSpool)
  write(sopri,*) 'next key is ',trim(key); call flush(sopri);
  call LeftAd(key)
  Line = key
  if (Line(1:1) == '*') cycle
  if (Line == ' ') cycle
  call UpCase(Line)
  select case ((Line(1:4)))

    case ('FILE')
      !========= FILE =============
      !> sets the reference wave function for DMRG-SCF with soMPSoo
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      call LeftAd(line)
      call fileorb(Line,refwfnfile)

    case('PUNC')
      !========= PUNCHMOs =============
      !> save intermediate-optimized MOs to file
      drop_intermediate_orbs = .true.;

    case('CONV')
      !========= CONV =============
      !> set the convergence (energy) threshold for the MPS optimization
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      read(Line,*,iostat=iError) mpsopt_param%conv_thresh
      write(sopri,*) 'MPS convergence threshold: ',mpsopt_param%conv_thresh; call flush(sopri);
      if (iError /= 0) call error(0)

    case('MAX_')
      !========= MAX_BOND_DIMENSION =============
      !> set the upper bond dimension M in the MPS optimization
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      read(Line,*,iostat=iError) mpsopt_param%M
      write(sopri,*) 'MPS max bond dimension: ',mpsopt_param%M; call flush(sopri);
      if (iError /= 0) call error(0)

    case('NSWE')
      !========= NSWEEPS =============
      !> set the number of sweeps for the MPS optimization
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      read(Line,*,iostat=iError) mpsopt_param%num_sweeps
      write(sopri,*) 'MPS number of sweeps: ',mpsopt_param%num_sweeps; call flush(sopri);
      if (iError /= 0) call error(0)

    case('ENTR')
      !========= ENTROPY =============
      !> enables the calculation of the entropy measures
      mpsopt_param%meas_entropy = .true.

    case('FIED')
      !========= FIEDLER =============
      !> enables the Fiedler vector ordering for the MPS optimization
      mpsopt_param%fiedler = .true.
    case('WF  ')
      !========= WF card =============
      !> WF,elec,sym,spin
      !> where

      ! @elec is the number of electrons
      ! @sym is the number of the irreducible representation
      ! @spin defines the spin symmetry, with spin=2*S (singlet=0, doublet=1, triplet=2, etc.)
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      j = j + 1;
      read(Line,*,iostat=iError) wf(j)%nelec,wf(j)%sym,wf(j)%spin
      write(sopri,*) 'WF(',j,'): parameters',wf(j)%nelec,wf(j)%sym,wf(j)%spin; call flush(sopri);
      if (iError /= 0) call error(0)  

    case('STAT')
      !========= STATE =============
      !> sets the number of states to convergence for each WF type
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      read(Line,*,iostat=iError) wf(j)%nstates
      if (iError /= 0) call error(0)
      if(.not.allocated(wf(j)%weights))then
          allocate(wf(j)%weights(wf(j)%nstates),stat=iError); 
          if( iError /= 0 ) call error(2)
          wf(j)%weights = 1;
      end if

    case('WEIG')
      !========= WEIGHTS =============
      !> determine the weights for each state
      if(.not.allocated(wf(j)%weights))then
          allocate(wf(j)%weights(wf(j)%nstates),stat=iError);
          if( iError /= 0 ) call error(2)
          wf(j)%weights = 1;
      end if
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      read(Line,*,iostat=iError) (wf(i)%weights(j),j=1,wf(i)%nstates)
      if (iError /= 0) call error(0)
   
    case('NO_2')
      !========= NO_2ndorder =============
      !> disables the coupling of orbital coefficients and MPS parameters in the optimization
      CP_integrals = .false.; ith_INTE     = -1;
    case('N_MI')
      !========= N_MICROiteration =============
      !> sets the number of micro iterations
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      read(Line,*,iostat=iError) MAX_iter;
      if (iError /= 0) call error(0);
    case('ALGO')
      !========= ALGORITHM =============
      !> sets NL solver algorithm
      if (.not. next_non_comment(LuSpool,Line)) call error(1)
      read(Line,*,iostat=iError) local_method,nstep
      if (iError /= 0) call error(0)

      estep = estep + nstep
      do j=istep,estep
        method(j)=local_method
      end do
      istep = istep + nstep
    case('A-A ')
      !========= a-a rotations =============
      !> enables the optimization of active-active rotations
      aaRotations = .true.
    case ('END ')

      exit

    case default
      write(u6,*) 'Unidentified key word  : '
      call FindErrorLine()
      call Quit_OnUserError()

  end select

end do
! END of Input

!> to do: set the total number of states (dmrg_nstates) + correct weights...
total_weight = 0; dmrg_nstates = 0;
write(sopri,*) ' size of wf',size(wf);call flush(sopri);
do j = 1, size(wf)
  total_weight = total_weight + sum(wf(j)%weights);
  dmrg_nstates = dmrg_nstates + wf(j)%nstates
end do
write(sopri,*) ' nstates=',dmrg_nstates;call flush(sopri);

allocate(dmrg_weight(dmrg_nstates),stat=ierr); if(iError/=0) call error(2)
!> set the global weights for each state
i = 1; 
do j = 1, size(wf)
  do k = 1, wf(j)%nstates
    dmrg_weight(i) = dble(wf(j)%weights(k))/dble(total_weight);
    i = i + 1;
  end do
end do

write(sopri,*) ' weights=',dmrg_weight(1:dmrg_nstates);call flush(sopri);


!> consistency check
num_macro_iterations = max(num_macro_iterations,estep)
num_macro_iterations = min(num_macro_iterations,MAX_macro_iterations)

!> allocate space for the state-specific energies of each WF type
do j = 1, size(wf)
  allocate(wf(j)%ss_energy(wf(j)%nstates)); wf(j)%ss_energy = 0;
end do

write(sopri,*) 'NL solver: number of iterations: ',num_macro_iterations
do j=1,num_macro_iterations
  write(sopri,*) ' iteration(',j,') => ',trim(method(j))
end do
call flush(sopri);

contains

subroutine error(code)

  integer(kind=iwp), intent(in) :: code

  select case(code)
  case(1)
    call WarningMessage(2,'Premature end of input file.')
  case(2) 
    call WarningMessage(2,'Problem with allocation.')
  case default
    call WarningMessage(2,'Read error during input preprocessing.')
  end select

  call Quit_OnUserError()
  call ABEND()

end subroutine

end subroutine rdinput
