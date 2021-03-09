************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2014, Giovanni Li Manni                                *
*               2019-2021, Oskar Weser                                 *
*               2021, Werner Dobrautz                                  *
************************************************************************
#include "macros.fh"
      module fciqmc
#ifdef _MOLCAS_MPP_
      use mpi
      use definitions, only: MPIInt
      use Para_Info, only: Is_Real_Par
#endif
      use definitions, only: wp
      use Para_Info, only: MyRank
#ifdef _NECI_
      use filesystem, only: chdir_
      use stdalloc, only: mxMem
      use fortran_strings, only: str
#endif
      use filesystem, only: getcwd_, get_errno_, strerror_,
     &    real_path, basename
      use linalg_mod, only: abort_
      use stdalloc, only: mma_allocate, mma_deallocate

      use rasscf_data, only: nAcPar, nAcPr2
      use general_data, only: nSym, nConf

      use CI_solver_util, only: wait_and_read, RDM_to_runfile

      use generic_CI, only: CI_solver_t

      implicit none
      save
      private
      public :: DoNECI, DoEmbdNECI, fciqmc_solver_t, tGUGA_in

      logical :: DoEmbdNECI = .false., DoNECI = .false.,
     &  tGUGA_in  = .false.

#ifdef _NECI_
      interface
        subroutine NECImain(fcidmp, input_name, MemSize, NECIen)
          use, intrinsic :: iso_fortran_env, only: int64, real64
          character(len=*), intent(in) :: fcidmp, input_name
          integer(int64), intent(in) :: MemSize
          real(real64), intent (out) :: NECIen
        end subroutine
      end interface
#endif


      type, extends(CI_solver_t) :: fciqmc_solver_t
        private
        logical :: tGUGA
      contains
        private
        procedure, public :: run => fciqmc_ctl
        procedure, public :: cleanup
      end type

      interface fciqmc_solver_t
        module procedure construct_FciqmcSolver_t
      end interface

      contains

      function construct_FciqmcSolver_t(tGUGA) result(res)
        logical, intent(in) :: tGUGA
        type(fciqmc_solver_t) :: res
        res%tGUGA = tGUGA
        write(6,*) ' NECI activated. List of Confs might get lengthy.'
        write(6,*) ' Number of Configurations computed by GUGA: ', nConf
        write(6,*) ' nConf variable is set to zero to avoid JOBIPH i/o'
        nConf= 0
      end function

!>  @brief
!>    Start and control FCIQMC.
!>
!>  @author Giovanni Li Manni, Oskar Weser
!>
!>  @details
!>  For meaning of global variables NTOT1, NTOT2, NACPAR
!>  and NACPR2, see src/Include/general.inc and src/Include/rasscf.inc.
!>  This routine will replace CICTL in FCIQMC regime.
!>  Density matrices are generated via double-run procedure in NECI.
!>  They are then dumped on arrays DMAT, DSPN, PSMAT, PAMAT to replace
!>  what normally would be done in DavCtl if NECI is not used.
!>  F_In is still generated by SGFCIN... in input contains
!>  only two-electron terms as computed in TRA_CTL2.
!>  In output it contains also the one-electron contribution
!>
!>  @paramin[in] actual_iter The actual iteration number starting at 0.
!>      This means 0 is 1A, 1 is 1B, 2 is 2 and so on.
!>  @paramin[in] CMO MO coefficients
!>  @paramin[in] DIAF Diagonal of Fock matrix useful for NECI
!>  @paramin[in] D1I_MO Inactive 1-dens matrix
!>  @paramin[in] TUVX Active 2-el integrals
!>  @paramin[inout] F_In Fock matrix from inactive density
!>  @paramin[inout] D1S_MO Average spin 1-dens matrix
!>  @paramin[out] DMAT Average 1 body density matrix
!>  @paramin[out] PSMAT Average symm. 2-dens matrix
!>  @paramin[out] PAMAT Average antisymm. 2-dens matrix
!>  @paramin[in] fake_run  If true the NECI run is not performed, but
!>    the RDMs are read from previous runs.
      subroutine fciqmc_ctl(
     &      this, actual_iter, CMO, DIAF, D1I_AO, D1A_AO,
     &      TUVX, F_IN, D1S_MO, DMAT, PSMAT, PAMAT)
      use general_data, only : iSpin, ntot, ntot1, ntot2, nAsh
      use rasscf_data, only : iter, lRoots, S, KSDFT, EMY,
     &    rotmax, Ener, Nac, nAcPar, nAcpr2

      use gas_data, only : ngssh, iDoGas, nGAS, iGSOCCX

      use fcidump_reorder, only : get_P_GAS, get_P_inp,ReOrFlag,ReOrInp
      use fcidump, only : make_fcidumps, transform
#include "rctfld.fh"
      class(fciqmc_solver_t), intent(in) :: this
      integer, intent(in) :: actual_iter
      real(wp), intent(in) ::
     &    CMO(nTot2), DIAF(nTot),
     &    D1I_AO(nTot2), D1A_AO(nTot2), TUVX(nAcpr2)
      real(wp), intent(inout) :: F_In(nTot1), D1S_MO(nAcPar)
      real(wp), intent(out) :: DMAT(nAcpar),
     &    PSMAT(nAcpr2), PAMAT(nAcpr2)

      real(wp) :: NECIen
      integer, allocatable :: permutation(:),
     &  GAS_spaces(:, :), GAS_particles(:, :)
      real(wp) :: orbital_E(nTot), folded_Fock(nAcPar)
#ifdef _MOLCAS_MPP_
      integer(MPIInt) :: error
#endif
      character(len=*), parameter ::
     &  ascii_fcidmp = 'FCIDUMP', h5_fcidmp = 'H5FCIDUMP'


! SOME DIRTY SETUPS
      S = 0.5_wp * dble(iSpin - 1)

      call check_options(lRoots, lRf, KSDFT)

! Produce a working FCIDUMP file
      if (ReOrFlag /= 0) then
        allocate(permutation(sum(nAsh(:nSym))))
        if (ReOrFlag >= 2) permutation(:) = get_P_inp(ReOrInp)
        if (ReOrFlag == -1) permutation(:) = get_P_GAS(nGSSH)
      end if

! This call is not side effect free, sets EMY and modifies F_IN
      call transform(actual_iter, CMO, DIAF, D1I_AO, D1A_AO, D1S_MO,
     &               F_IN, orbital_E, folded_Fock)

! Fortran Standard 2008 12.5.2.12:
! Allocatable actual arguments that are passed to
! non-allocatable, optional dummy arguments are **not** present.
      call make_fcidumps(ascii_fcidmp, h5_fcidmp,
     &                   orbital_E, folded_Fock, TUVX, EMY, permutation)

! Run NECI
#ifdef _MOLCAS_MPP_
      if (is_real_par()) call MPI_Barrier(MPI_COMM_WORLD, error)
#endif

      if (iDoGAS) then
        call mma_allocate(GAS_spaces, nGAS, nSym)
        GAS_spaces(:, :) = nGSSH(: nGAS, : nSym)
        call mma_allocate(GAS_particles, nGAS, nGAS)
        GAS_particles(:, :) = iGSOCCX(: nGAS, : nGAS)
      end if

      call run_neci(DoEmbdNECI, actual_iter == 1,
     &  ascii_fcidmp, h5_fcidmp,
     &  GAS_spaces=GAS_spaces, GAS_particles=GAS_particles,
     &  reuse_pops=actual_iter >= 5 .and. abs(rotmax) < 1d-2,
     &  NECIen=NECIen,
     &  D1S_MO=D1S_MO, DMAT=DMAT, PSMAT=PSMAT, PAMAT=PAMAT,
     &  tGUGA=this%tGUGA)
! NECIen so far is only the energy for the GS.
! Next step it will be an array containing energies for all the optimized states.
      ENER(1 : lRoots, iter) = NECIen

      if (nAsh(1) /= nac) call dblock(dmat)


      if (allocated(GAS_spaces)) then
          call mma_deallocate(GAS_spaces)
          call mma_deallocate(GAS_particles)
      end if
      end subroutine fciqmc_ctl


      subroutine run_neci(DoEmbdNECI, fake_run,
     &      ascii_fcidmp, h5_fcidmp,
     &      reuse_pops,
     &      NECIen, D1S_MO, DMAT, PSMAT, PAMAT, tGUGA,
     &      GAS_spaces, GAS_particles)
        use fciqmc_make_inp, only: make_inp
        logical, intent(in) :: DoEmbdNECI, fake_run, reuse_pops
        character(len=*), intent(in) :: ascii_fcidmp, h5_fcidmp
        real(wp), intent(out) :: NECIen, D1S_MO(nAcPar), DMAT(nAcpar),
     &      PSMAT(nAcpr2), PAMAT(nAcpr2)
        logical, intent(in) :: tGUGA
        integer, intent(in), optional ::
     &      GAS_spaces(:, :), GAS_particles(:, :)
        real(wp), save :: previous_NECIen = 0.0_wp

        character(len=*), parameter :: input_name = 'FCINP',
     &    energy_file = 'NEWCYCLE'

#ifdef _WARNING_WORKAROUND_
        ! there are wrong maybe-unitialized warnings otherwise.
        ! (unitialized codepaths lead to abortion).
        NECIen = huge(NECIen)
#endif

        if (fake_run) then
          NECIen = previous_NECIen
        else if (DoEmbdNECI) then
            call make_inp(input_name, readpops=reuse_pops, tGUGA=tGUGA,
     &          GAS_spaces=GAS_spaces, GAS_particles=GAS_particles)
#ifdef _NECI_
            write(6,*) 'NECI called automatically within Molcas!'
            if (myrank /= 0) call chdir_('..')
            call necimain(real_path(ascii_fcidmp),
     &                    real_path(input_name),
     &                    MxMem, NECIen)
            if (myrank /= 0) call chdir_('tmp_'//str(myrank))
#else
            call WarningMessage(2, 'EmbdNECI is given in input, '//
     &'so the embedded NECI should be used. Unfortunately MOLCAS was '//
     &'not compiled with embedded NECI. Please use -DNECI=ON '//
     &'for compiling or use an external NECI.')
#endif
        else
            call make_inp(input_name, readpops=.false., tGUGA=tGUGA,
     &              FCIDUMP_name=basename(real_path(ascii_fcidmp)),
     &              GAS_spaces=GAS_spaces, GAS_particles=GAS_particles)
            if (myrank == 0) then
              call write_ExNECI_message(input_name, ascii_fcidmp,
     &                                  h5_fcidmp, energy_file, tGUGA)
            end if
            call wait_and_read(energy_file, NECIen)
        end if
        previous_NECIen = NECIen
        call get_neci_RDM(D1S_MO, DMAT, PSMAT, PAMAT, tGUGA)
        call RDM_to_runfile(DMAT, D1S_MO, PSMAT, PAMAT)
      end subroutine run_neci

      subroutine cleanup(this)
        use fciqmc_make_inp, only : make_inp_cleanup => cleanup
        use fciqmc_read_RDM, only : read_RDM_cleanup => cleanup
        use fcidump, only : fcidump_cleanup => cleanup
        class(fciqmc_solver_t), intent(inout) :: this
        unused_var(this)
        call make_inp_cleanup()
        call read_RDM_cleanup()
        call fcidump_cleanup()
      end subroutine cleanup

      subroutine check_options(lroots, lRf, KSDFT)
        integer, intent(in) :: lroots
        logical, intent(in) :: lRf
        character(len=*), intent(in) :: KSDFT
        logical :: Do_ESPF
        if (lroots > 1) then
          call abort_('FCIQMC does not support State Average yet!')
        end if
        call DecideOnESPF(Do_ESPF)
        if ( lRf .or. KSDFT /= 'SCF' .or. Do_ESPF) then
          call abort_('FCIQMC does not support Reaction Field yet!')
        end if
      end subroutine check_options

      subroutine write_ExNECI_message(
     &      input_name, ascii_fcidmp, h5_fcidmp, energy_file, tGUGA)
        character(len=*), intent(in) :: input_name, ascii_fcidmp,
     &          h5_fcidmp, energy_file
        logical, intent(in) :: tGUGA
        character(len=1024) :: WorkDir
        integer :: err

        call getcwd_(WorkDir, err)
        if (err /= 0) write(6, *) strerror_(get_errno_())

        if (tGUGA) then
            write(6,'(A)')'Run spin-free GUGA NECI externally.'
        else
            write(6,'(A)')'Run NECI externally.'
        end if
        write(6,'(A)')'Get the (example) NECI input:'
        write(6,'(4x, A, 1x, A, 1x, A)')
     &    'cp', real_path(input_name), '$NECI_RUN_DIR'
        write(6,'(A)')'Get the ASCII formatted FCIDUMP:'
        write(6,'(4x, A, 1x, A, 1x, A)')
     &    'cp', real_path(ascii_fcidmp), '$NECI_RUN_DIR'
        write(6,'(A)')'Or the HDF5 FCIDUMP:'
        write(6,'(4x, A, 1x, A, 1x, A)')
     &    'cp', real_path(h5_fcidmp), '$NECI_RUN_DIR'
        write(6, *)
        write(6,'(A)') "When finished do:"
        if (tGUGA) then
          write(6,'(4x, A)') 'cp PSMAT PAMAT DMAT '//trim(WorkDir)
        else
          write(6,'(4x, A)')
     &      'cp TwoRDM_aaaa.1 TwoRDM_abab.1 TwoRDM_abba.1 '//
     &      'TwoRDM_bbbb.1 TwoRDM_baba.1 TwoRDM_baab.1 '//trim(WorkDir)
        end if
        write(6,'(4x, A)')
     &    'echo $your_RDM_Energy > '//real_path(energy_file)
        call xflush(6)
      end subroutine write_ExNECI_message

!> Generate density matrices for Molcas
!>   Neci density matrices are stored in Files TwoRDM_**** (in spacial orbital basis).
!>   I will be reading them from those formatted files for the time being.
!>   Next it will be nice if NECI prints them out already in Molcas format.
      subroutine get_neci_RDM(D1S_MO, DMAT, PSMAT, PAMAT, tGUGA)
        use fciqmc_read_RDM, only : read_neci_RDM, read_neci_GUGA_RDM
        real*8, intent(out) ::
     &      D1S_MO(nAcPar), DMAT(nAcpar),
     &      PSMAT(nAcpr2), PAMAT(nAcpr2)
        logical, intent(in) :: tGUGA
        if (tGUGA) then
            ! for GUGA we only need the spin-free 1 and 2-RDM.
            ! The spin density D1S_MO is set to zero.
            call read_neci_GUGA_RDM(DMAT, D1S_MO, PSMAT, PAMAT)
        else
            call read_neci_RDM(DMAT, D1S_MO, PSMAT, PAMAT)
        end if
      end subroutine get_neci_RDM
      end module fciqmc