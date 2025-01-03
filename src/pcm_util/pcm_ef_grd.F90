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

subroutine PCM_EF_grd(Grad,nGrad)

use Basis_Info, only: dbsc, nCnttp
use Center_Info, only: dc
use PCM_arrays, only: dCntr, dPnt, PCM_SQ, PCMiSph, PCMTess
use Gateway_global, only: PrPrt
use Symmetry_Info, only: nIrrep
use rctfld_module, only: nS, nTS
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6
use PCM_alaska, only: lSA,DSA_AO,PCM_SQ_ind,PCM_SQ_ori,def_solv,DSS_ori
use NAC, only: isNAC

implicit none
integer(kind=iwp), intent(in) :: nGrad
real(kind=wp), intent(inout) :: Grad(nGrad)
integer(kind=iwp) :: i, iTile, jCnt, jCnttp, MaxAto, mCnt, nc, nComp, ndc, nDens, nOrdOp
real(kind=wp) :: EF_Temp(3), Z
logical(kind=iwp) :: Save_tmp, Found
integer(kind=iwp), allocatable :: lOper(:)
real(kind=wp), allocatable :: Chrg(:), Cord(:,:), D1ao(:), EF(:,:,:), FactOp(:), tmpchg(:)

!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
Save_tmp = PrPrt
PrPrt = .true.
nOrdOp = 1
nComp = (nOrdOp+1)*(nOrdOp+2)/2
call mma_allocate(EF,nComp,2,nTs,label='EF')
!                                                                      *
!***********************************************************************
!                                                                      *
call Get_nAtoms_All(MaxAto)

call mma_allocate(Cord,3,MaxAto)
call mma_allocate(Chrg,MaxAto)

ndc = 0
nc = 1
do jCnttp=1,nCnttp
  if (dbsc(jCnttp)%Aux) cycle
  Z = dbsc(jCnttp)%Charge
  mCnt = dbsc(jCnttp)%nCntr
  do jCnt=1,mCnt
    ndc = ndc+1
    do i=0,nIrrep/dc(ndc)%nStab-1
      call OA(dc(ndc)%iCoSet(i,0),dbsc(jCnttp)%Coor(1:3,jCnt),Cord(1:3,nc))
      Chrg(nc) = Z
      nc = nc+1
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the the electric field on the tiles

! 1) The nuclear contribution

do iTile=1,nTs
  call EFNuc(PCMTess(1,iTile),Chrg,Cord,MaxAto,EF_temp,nOrdOp)
  EF(:,1,iTile) = EF_Temp(:)
  EF(:,2,iTile) = Zero
end do

call mma_deallocate(Cord)
call mma_deallocate(Chrg)

! 2) The electronic contribution

! Get the total 1st order AO density matrix

call Qpg_dArray('D1ao',Found,nDens)
if (Found .and. (nDens /= 0)) then
  call mma_allocate(D1ao,nDens,Label='D1ao')
else
  write(u6,*) 'pcm_ef_grd: D1ao not found.'
  call Abend()
end if
call Get_dArray_chk('D1ao',D1ao,nDens)

if (lSA) then
  !! For SA-CASSCF, one thing we need here (for def_solv=3) is
  !! V^N*q^N + (V^N+V^SA)*q^SS + V^SS*(q^N+q^SA)
  !! with V^N = EF(:,1,:), V^S[A,S] = EF(:,2,:), q^N = PCM_SQ(1,:) = PCM_SQ_ind(1,:),
  !! q^SA = PCM_SQ(2,:), q^SS = PCM_SQ_ind(2,:)
  call mma_allocate(tmpchg,nTs,Label='tmpchg')
  !! V^N*q^N contribution
  if (.not.isNAC) then
    tmpchg(:) = PCM_SQ(2,:)
    PCM_SQ(2,:) = 0.0d+00 !! No electron contributions
    call Cmbn_EF_DPnt(EF,nTs,dPnt,MaxAto,dCntr,nS,PCMiSph,PCM_SQ,Grad,nGrad)
    PCM_SQ(2,:) = tmpchg(:) !! Restore electron
  end if

  !! The first Cmbn_EF_DPnt will compute (V^N+V^SA)*q^SS, so computes the ESP with D^SA,
  !! and q^SS and q^SA should be swapped
  call dcopy_(nDens,DSA_AO,1,D1ao,1)
  call dswap_(2*nTS,PCM_SQ,1,PCM_SQ_ind,1)
  tmpchg(:) = PCM_SQ(1,:)
  PCM_SQ(1,:) = 0.0d+00 ! we do not need q^N here
end if

call mma_allocate(FactOp,nTs)
call mma_allocate(lOper,nTs)
FactOp(:) = One
lOper(:) = 255

!! this drv1_pcm is not parallelized
call Drv1_PCM(FactOp,nTs,D1ao,nDens,PCMTess,lOper,EF,nOrdOp)

call mma_deallocate(lOper)
call mma_deallocate(FactOp)
call mma_deallocate(D1ao)
!                                                                      *
!***********************************************************************
!                                                                      *
! Now form the correct combinations

call Cmbn_EF_DPnt(EF,nTs,dPnt,MaxAto,dCntr,nS,PCMiSph,PCM_SQ,Grad,nGrad)

if (lSA) then
  !! restore q^SS and q^SA correctly
  PCM_SQ(1,:) = tmpchg(:)
  call dswap_(2*nTS,PCM_SQ,1,PCM_SQ_ind,1)

  call mma_allocate(FactOp,nTs)
  call mma_allocate(lOper,nTs)
  FactOp(:) = One
  lOper(:) = 255
  call mma_allocate(D1ao,nDens,Label='D1ao')

  !! Here computes V^SS*(q^N+q^SA)
  EF = 0.0d+00
  call Get_D1ao_Var(D1ao,nDens)
  call Drv1_PCM(FactOp,nTs,D1ao,nDens,PCMTess,lOper,EF,nOrdOp)
  call Cmbn_EF_DPnt(EF,nTs,dPnt,MaxAto,dCntr,nS,PCMiSph,PCM_SQ,Grad,nGrad)

  if (.not.isNAC) then
    EF = 0.0d+00
    if (def_solv==1 .or. def_solv==3) then
      ! -V^SA*q^{e,SA}/2 term (= -V^SA*q^SA)
      call dcopy_(nDens,DSA_AO,1,D1ao,1)
      call Drv1_PCM(FactOp,nTs,D1ao,nDens,PCMTess,lOper,EF,nOrdOp)
      EF = -EF

      tmpchg(:) = PCM_SQ(1,:)
      PCM_SQ(1,:) = 0.0d+00
      call Cmbn_EF_DPnt(EF,nTs,dPnt,MaxAto,dCntr,nS,PCMiSph,PCM_SQ,Grad,nGrad)
      PCM_SQ(1,:) = tmpchg(:)
    else if (def_solv==4) then
      ! -V^SS*q^{e,SA}/2 term
      call dcopy_(nDens,DSA_AO,1,D1ao,1)
      call Drv1_PCM(FactOp,nTs,D1ao,nDens,PCMTess,lOper,EF,nOrdOp) ! V^SA
      EF = -EF

      tmpchg(:) = PCM_SQ_ori(1,:)
      PCM_SQ_ori(1,:) = 0.0d+00
      PCM_SQ_ori(2,:) = PCM_SQ_ori(2,:)*0.5d+00
      ! V^SA*q^SS
      call Cmbn_EF_DPnt(EF,nTs,dPnt,MaxAto,dCntr,nS,PCMiSph,PCM_SQ_ori,Grad,nGrad)
      PCM_SQ_ori(1,:) = tmpchg(:)
      PCM_SQ_ori(2,:) = PCM_SQ_ori(2,:)*2.0d+00

      EF = 0.0D+00
      call dcopy_(nDens,DSS_ori,1,D1ao,1)
      call Drv1_PCM(FactOp,nTs,D1ao,nDens,PCMTess,lOper,EF,nOrdOp) ! V^SS
      EF = -EF

      tmpchg(:) = PCM_SQ(1,:)
      PCM_SQ(1,:) = 0.0d+00
      PCM_SQ(2,:) = PCM_SQ(2,:)*0.5d+00
      ! V^SS*q^SA
      call Cmbn_EF_DPnt(EF,nTs,dPnt,MaxAto,dCntr,nS,PCMiSph,PCM_SQ,Grad,nGrad)
      PCM_SQ(1,:) = tmpchg(:)
      PCM_SQ(2,:) = PCM_SQ(2,:)*2.0d+00
    else if (def_solv==5) then
    end if
  end if

  call mma_deallocate(tmpchg)
  call mma_deallocate(lOper)
  call mma_deallocate(FactOp)
  call mma_deallocate(D1ao)
end if

!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(EF)
PrPrt = Save_tmp
call Free_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine PCM_EF_grd
