module soMPSoo_interface

implicit none

public soMPSoo_init
public soMPSoo_driver

contains

subroutine soMPSoo_init(              &
                        nsym_ext,     &
                        nactel_ext,   &
                        nstates_ext,  &
                        spin_ext,     &
                        lsym_ext,     &
                        max_macro_ext,&
                        max_micro_ext,&
                        ogra_thr_ext, &
                        nfro_ext,     &
                        nish_ext,     &
                        nash_ext,     &
                        nssh_ext,     &
                        nnorbt_ext,   &
                        nnashx_ext,   &
                        nbas_ext,     &
                        lwrk_ext,     &
                        weights_ext,  &
                        potnuc_ext,   &
                        method_ext,   &
                        lupri_ext,    &
                        srdft_ext     &
                       )
  use global_control
  use orbopt_header

  integer, intent(in)          :: nsym_ext, nactel_ext, nstates_ext, spin_ext, lsym_ext
  integer, intent(in)          :: max_macro_ext, max_micro_ext, lupri_ext
  integer, intent(in)          :: nnorbt_ext, nnashx_ext, lwrk_ext
  integer, intent(in)          :: nfro_ext(nsym_ext)
  integer, intent(in)          :: nish_ext(nsym_ext)
  integer, intent(in)          :: nash_ext(nsym_ext)
  integer, intent(in)          :: nssh_ext(nsym_ext)
  integer, intent(in)          :: nbas_ext(nsym_ext)
  real*8,  intent(in)          :: weights_ext(nstates_ext), potnuc_ext
  real*8,  intent(in)          :: ogra_thr_ext
  logical, intent(in)          :: srdft_ext

  character(len=7), intent(in) :: method_ext

  integer             :: i
  real*8              :: d1

  !> print unit
  sopri = lupri_ext

  !> logical flags
  if(sum(nfro_ext) > 0) ifcore  = .true.

  do_soMPSoo_srdft = srdft_ext
  DMRG_SCF         = .true.
  LMCSCF           = .false.
  !> LR stuff
  DMRG_LAG         = .false.
  LRirreps         = .true.
  CKP_LR           = .true.
  Lcorrected_phase = .false.

  !> integer
  ntdms         = 0
  !> target state for SA-gradient is state # 1
  dmrg_rlxstate = 1

  FCIDUMP = ''; FCIDUMP(1:7)  = 'FCIDUMP'

  !> number of point group irreps
  orb%nsub   = nsym_ext

  !> allocate space for arrays
  allocate(orb%frozen(orb%nsub));     orb%frozen=0
  allocate(orb%closed(orb%nsub));     orb%closed=0
  allocate(orb%   occ(orb%nsub));     orb%   occ=0
  allocate(orb%   act(orb%nsub));     orb%   act=0
  allocate(orb% extnl(orb%nsub));     orb% extnl=0
  allocate(orb% total(orb%nsub));     orb% total=0
                                      orb%   ele=0
  allocate(orbbk%frozen(orb%nsub)); orbbk%frozen=0
  allocate(orbbk%closed(orb%nsub)); orbbk%closed=0
  allocate(orbbk%   occ(orb%nsub)); orbbk%   occ=0
  allocate(orbbk%   act(orb%nsub)); orbbk%   act=0
  allocate(orbbk% extnl(orb%nsub)); orbbk% extnl=0
  allocate(orbbk% total(orb%nsub)); orbbk% total=0
                                    orbbk%   ele=0

  allocate(orb%grouptable(orb%nsub,orb%nsub)); orb%grouptable=0

  allocate(off%iish(orb%nsub)); off%iish=0
  allocate(off%iash(orb%nsub)); off%iash=0
  allocate(off%issh(orb%nsub)); off%issh=0
  allocate(off%iocc(orb%nsub)); off%iocc=0
  allocate(off%iorb(orb%nsub)); off%iorb=0
  allocate(off%icmo(orb%nsub)); off%icmo=0
  allocate(off%ibas(orb%nsub)); off%ibas=0

  allocate(bas%nbas(orb%nsub)); bas%nbas=0

  !> orbital data
  orb%frozen = nfro_ext
  orb%closed = nish_ext
  orb%act    = nash_ext
  orb%occ    = orb%frozen + orb%closed + orb%act
  orb%extnl  = nssh_ext
  orb%total  = orb%frozen + orb%closed + orb%act + orb%extnl

  orb%ele    = nactel_ext

  nnorbt     = nnorbt_ext
  nnashx     = nnashx_ext

  !> basis functions
  bas%nbas   = nbas_ext

  !> offset arrays
  do i=2,orb%nsub
   off%iish(I)  =off%iish(i-1)  +orb%closed(i-1)
   off%iash(I)  =off%iash(i-1)  +orb%act   (i-1)
   off%issh(I)  =off%issh(i-1)  +orb%extnl (i-1)
   off%iocc(I)  =off%iocc(i-1)  +orb%occ   (i-1)
   off%iorb(I)  =off%iorb(i-1)  +orb%total (i-1)
   off%ibas(I)  =off%ibas(i-1)  +bas%nbas  (i-1)
   off%icmo(I)  =off%icmo(i-1)  +orb%total (i-1) * bas%nbas(i-1)
  end do

  !> various data
  lwrk       = lwrk_ext
  potnuc     = potnuc_ext

  !> wave function data
  tspin        = spin_ext
  tirep        = lsym_ext
  dmrg_nstates = nstates_ext
  dmrg_weight  = weights_ext
  d1=0.0d0
  do i=1,dmrg_nstates
    d1 = d1 + dmrg_weight(i)
  end do
  dmrg_weight = dmrg_weight / d1

  !> DMRG driver data
  dmrg_binary_name  = 'maquis'
  dmrg_reorder_name = 'MASORB.orb'
  dmrg_output_name  = 'dmrg.out'
  ifmpirun          = .true.
  nproc_in_mpirun   = ''; nproc_in_mpirun = '1'
  if(dmrg_nstates > 10) stop 'interface only for < 11 states'
  do i=1,dmrg_nstates
    write(dmrg_inputs(i),'(a11,i1)') 'dmrg-input-',i-1
  end do

  !> Solver defaults
  thrs%r         = 1.0d-5   ! rotation (use norm)
  if(ogra_thr_ext /= 0.0d0) thrs%r = ogra_thr_ext
  thrs%s         = 1.0d-8   ! symmetry
  thrs%e         = 1.0d-8   ! energy
  thrs%d         = 1.0d-12  ! davidson
  thrs%c         = 1.0d-4   ! coupling
  thrs%g         = 1.0d-2   ! gradients
  Nmps_update    = 1000
  Nint_update    = 0
  step_damping   = 0.0d0
  restart_maquis = .false.
  CP_integrals   = .true.
  ith_INTE       = 2

  !> micro / macro iterations
  Max_iter         = max_micro_ext
  ncycle           = max_macro_ext

  do i = 1, ncycle
    method(i) = method_ext
  end do

  !> set active / occupied orbitals
  nact=0
  nocc=0
  do i=1,orb%nsub
    nact    = nact    + orb%act(i)
    nocc    = nocc    + orb%occ(i)
  end do

  !> setup multiplication table
  call group_table()

  !> say hello
  call print_orbopt_header(sopri)

  if(orb%nsub > 1) stop "soMPSoo interface only tested for C1 symmetry"

end subroutine soMPSoo_init

subroutine soMPSoo_driver(             &
                          task,        &
                          dv,          &
                          pv,          &
                          intT,        &
                          intU,        &
                          CMO,         &
                          work,        &
                          norb_ext,    &
                          core_energy, &
                          fname        &
                         )
  use global_control
  use integral_interface
  use dmrgscf, only: dmrgscf_main

  real*8, optional,          intent(in   ) :: core_energy
  real*8, optional,          intent(in   ) :: dv(*)
  real*8, optional,          intent(in   ) :: pv(*)
  real*8, optional,          intent(in   ) :: intT(*)
  real*8, optional,          intent(in   ) :: intU(*)
  real*8, optional,          intent(inout) :: CMO(*)
  real*8, optional,          intent(inout) :: work(*)
  integer,optional,          intent(in   ) :: norb_ext
  character(len=8),          intent(in   ) :: task
  character(len=18),optional,intent(in   ) :: fname

  select case(task)
    case('optimize')
      call dmrgscf_main( relOPT = .false.,&
                         CMO    = CMO,    &
                         work   = work    &
                       )
    case('fcidump ')
        call fcidump_write_full(                             &
                                ndim        = norb_ext,      &
                                T           = intT,          &
                                U           = intU,          &
                                core_energy = core_energy,   &
                                fcidump1    = fname,         &
                                nsub        = orb%nsub,      &
                                tot         = orb%total,     &
                                mtab        = orb%grouptable &
                               )
    case('get rdms')
    case default
      stop 'no task assigned for soMPSoo'
  end select

end subroutine soMPSoo_driver

end module soMPSoo_interface
