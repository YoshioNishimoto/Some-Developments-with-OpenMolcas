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
*               2019, Oskar Weser                                      *
************************************************************************
      module fciqmc
#ifdef _MOLCAS_MPP_
      use mpi
#endif
      use definitions, only: wp, MPIInt, int64, real64
      use filesystem, only: chdir_, getcwd_, get_errno_, strerror_,
     &    real_path
      use fortran_strings, only: str
      use stdalloc, only : mma_allocate, mma_deallocate, mxMem

      use rasscf_data, only: iter, lRoots, nRoots, iRoot, EMY,
     &    S, KSDFT, rotmax, Ener, iAdr15, Weight, nAc, nAcPar, nAcPr2
      use general_data, only: iSpin, nSym, nConf, JobIPH,
     &    ntot, ntot1, ntot2, nAsh, nBas
      use gugx_data, only: IfCAS
      use gas_data, only: ngssh, iDoGas, nGAS, iGSOCCX

      use CI_solver_util, only: wait_and_read, abort_, RDM_to_runfile
      use fciqmc_read_RDM, only: read_neci_RDM

      use generic_CI, only: CI_solver_t, unused

      implicit none
      save
      private
      public :: DoNECI, DoEmbdNECI, fciqmc_solver_t
      logical :: DoEmbdNECI = .false., DoNECI = .false.
#include "para_info.fh"

      interface
          integer function isfreeunit(iseed)
              integer, intent(in) :: iseed
          end function

          subroutine NECImain(fcidmp, input_name, MemSize, NECIen)
            import :: int64, real64
            character(*), intent(in) :: fcidmp, input_name
            integer(int64), intent(in) :: MemSize
            real(real64), intent (out) :: NECIen
          end subroutine
      end interface


      type, extends(CI_solver_t) :: fciqmc_solver_t
      contains
        procedure, nopass :: init
        procedure, nopass :: run => fciqmc_ctl
        procedure, nopass :: cleanup
      end type


      contains

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
      subroutine fciqmc_ctl(actual_iter, CMO, DIAF, D1I_AO, D1A_AO,
     &                      TUVX, F_IN, D1S_MO, DMAT, PSMAT, PAMAT)
      use fcidump_reorder, only : get_P_GAS, get_P_inp,ReOrFlag,ReOrInp
      use fcidump, only : make_fcidumps, transform
#include "output_ras.fh"
#include "rctfld.fh"
#include "timers.fh"
      integer, intent(in) :: actual_iter
      real(wp), intent(in) ::
     &    CMO(nTot2), DIAF(nTot),
     &    D1I_AO(nTot2), D1A_AO(nTot2), TUVX(nAcpr2)
      real(wp), intent(inout) :: F_In(nTot1), D1S_MO(nAcPar)
      real(wp), intent(out) :: DMAT(nAcpar),
     &    PSMAT(nAcpr2), PAMAT(nAcpr2)
      real(wp) :: NECIen
      integer :: jRoot
      integer, allocatable :: permutation(:)
      real(wp) :: orbital_E(nTot), folded_Fock(nAcPar)
#ifdef _MOLCAS_MPP_
      integer(MPIInt) :: error
#endif

      parameter(ROUTINE = 'FCIQMC_clt')
      character(*), parameter ::
     &  ascii_fcidmp = 'FCIDUMP', h5_fcidmp = 'H5FCIDUMP'

      call qEnter(routine)

! SOME DIRTY SETUPS
      S = 0.5_wp * dble(iSpin - 1)

      call check_options(lRoots, lRf, KSDFT, iDoGAS, iGSOCCX, nGAS)

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

      if (iDoGAS) call write_GASORB(nGSSH, permutation)

! Run NECI
      call Timing(Rado_1, Swatch, Swatch, Swatch)
#ifdef _MOLCAS_MPP_
      if (is_real_par()) call MPI_Barrier(MPI_COMM_WORLD, error)
#endif

      call run_neci(DoEmbdNECI, actual_iter == 1,
     &  ascii_fcidmp, h5_fcidmp, doGAS=iDoGAS,
     &  reuse_pops=actual_iter >= 5 .and. abs(rotmax) < 1d-2,
     &  NECIen=NECIen,
     &  D1S_MO=D1S_MO, DMAT=DMAT, PSMAT=PSMAT, PAMAT=PAMAT)
! NECIen so far is only the energy for the GS.
! Next step it will be an array containing energies for all the optimized states.
      do jRoot = 1, lRoots
        ENER(jRoot, ITER) = NECIen
      end do

      if (nAsh(1) /= nac) call dblock(dmat)


      call Timing(Rado_2, Swatch, Swatch, Swatch)
      Rado_2 = Rado_2 - Rado_1
      Rado_3 = Rado_3 + Rado_2

      call qExit(routine)
      end subroutine fciqmc_ctl


      subroutine run_neci(DoEmbdNECI, fake_run,
     &      ascii_fcidmp, h5_fcidmp,
     &      reuse_pops,
     &      NECIen, D1S_MO, DMAT, PSMAT, PAMAT, doGAS)
        use fciqmc_make_inp, only: make_inp
        logical, intent(in) :: DoEmbdNECI, fake_run, reuse_pops
        character(*), intent(in) :: ascii_fcidmp, h5_fcidmp
        real(wp), intent(out) :: NECIen, D1S_MO(nAcPar), DMAT(nAcpar),
     &      PSMAT(nAcpr2), PAMAT(nAcpr2)
        logical, intent(in), optional :: doGAS
        logical :: doGAS_
        real(wp), save :: previous_NECIen = 0.0_wp

        character(*), parameter :: input_name = 'FCINP',
     &    energy_file = 'NEWCYCLE'

        if (present(doGAS)) then
          doGAS_ = doGAS
        else
          doGAS_ = .false.
        end if

        if (fake_run) then
          NECIen = previous_NECIen
        else
          if (DoEmbdNECI) then
            call make_inp(input_name, doGAS=doGAS_, readpops=reuse_pops)
#ifdef _NECI_
            write(6,*) 'NECI called automatically within Molcas!'
            if (myrank /= 0) call chdir_('..')
            call necimain(
     &        real_path(ascii_fcidmp), real_path(input_name),
     &                  MxMem, NECIen)
            if (myrank /= 0) call chdir_('tmp_'//str(myrank))
#else
            call WarningMessage(2, 'EmbdNECI is given in input, '//
     &'so the embedded NECI should be used. Unfortunately MOLCAS was '//
     &'not compiled with embedded NECI. Please use -DNECI=ON '//
     &'for compiling or use an external NECI.')
#endif
          else
            call make_inp(input_name, doGAS=doGAS_)
            if (myrank == 0) then
              call write_ExNECI_message(input_name, ascii_fcidmp,
     &                                  h5_fcidmp, energy_file)
            end if
            call wait_and_read(energy_file, NECIen)
          end if
          previous_NECIen = NECIen
        end if
        call read_neci_RDM(DMAT, D1S_MO, PSMAT, PAMAT)
        call RDM_to_runfile(DMAT, D1S_MO, PSMAT, PAMAT)
      end subroutine run_neci

      subroutine cleanup()
        use fciqmc_make_inp, only : make_inp_cleanup => cleanup
        use fciqmc_read_RDM, only : read_RDM_cleanup => cleanup
        use fcidump, only : fcidump_cleanup => cleanup
        call make_inp_cleanup()
        call read_RDM_cleanup()
        call fcidump_cleanup()
      end subroutine cleanup

      subroutine init()
! For NECI only orbital related arrays are allowed to be stored.
! Arrays of nConf size need to be avoided.
! For this reason set nConf to zero.
        write(6,*) ' NECI activated. List of Confs might get lengthy.'
        write(6,*) ' Number of Configurations computed by GUGA: ', nConf
        write(6,*) ' nConf variable is set to zero to avoid JOBIPH i/o'
        nConf= 0
      end subroutine


      subroutine check_options(lroots, lRf, KSDFT,
     &      DoGAS, iGSOCCX, nGAS)
        integer, intent(in) :: lroots, iGSOCCX(:, :),nGAS
        logical, intent(in) :: lRf, DoGAS
        character(*), intent(in) :: KSDFT
        logical :: Do_ESPF
        if (lroots > 1) then
          call abort_('FCIQMC does not support State Average yet!')
        end if
        call DecideOnESPF(Do_ESPF)
        if ( lRf .or. KSDFT /= 'SCF' .or. Do_ESPF) then
          call abort_('FCIQMC does not support Reaction Field yet!')
        end if
        if (DoGAS) then
          if (.not. all(iGSOCCX(:nGAS, 1) == iGSOCCX(:nGAS, 2))) then
            call abort_('Only disconnected GAS spaces are '//
     &        'currently supported in FCIQMC.')
          end if
        end if
      end subroutine check_options

      subroutine write_ExNECI_message(
     &      input_name, ascii_fcidmp, h5_fcidmp, energy_file)
        character(*), intent(in) :: input_name, ascii_fcidmp, h5_fcidmp,
     &          energy_file
        character(1024) :: WorkDir
        integer :: err

        call getcwd_(WorkDir, err)
        if (err /= 0) write(6, *) strerror_(get_errno_())

        write(6,'(A)')'Run NECI externally.'
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
        write(6,'(4x, A)')
     &    'cp TwoRDM_aaaa.1 TwoRDM_abab.1 TwoRDM_abba.1 '//
     &    'TwoRDM_bbbb.1 TwoRDM_baba.1 TwoRDM_baab.1 '//trim(WorkDir)
        write(6,'(4x, A)')
     &    'echo $your_RDM_Energy > '//real_path(energy_file)
        call xflush(6)
      end subroutine write_ExNECI_message

      subroutine write_GASORB(GAS_spaces, permutation)
        integer, intent(in) :: GAS_spaces(:, :)
        integer, intent(in), optional :: permutation(:)
        integer, parameter :: arbitrary_magic_number = 42
        integer :: i, GAS_ORB(sum(GAS_spaces)), iGAS, iSym, file_id

        GAS_ORB(:) = [(((iGAS, i = 1, GAS_spaces(iGAS, iSym)),
     &                 iGAS = 1, size(GAS_spaces, 1)), iSym = 1, nSym)]

        if (present(permutation)) GAS_ORB = GAS_ORB(permutation)

        file_id = arbitrary_magic_number
        file_id = isfreeunit(file_id)
        call molcas_open(file_id, 'GASOrbs')
          do i = 1, size(GAS_ORB)
            write(file_id,'(I0, A)', advance='no') GAS_ORB(i), ','
          end do
        close(file_id)
      end subroutine

      end module fciqmc
