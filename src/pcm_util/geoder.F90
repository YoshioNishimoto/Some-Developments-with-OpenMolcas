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

subroutine GeoDer(nAt,Cond,nTs,nS,Eps,Sphere,ISphe,NOrd,Tessera,Q1,Q2,DerDM,Grd,DerTes,DerPunt,DerRad,DerCentr)

use Constants, only: Zero, Half, Two, one
use Definitions, only: wp, iwp
use PCM_alaska, only: lSA,PCM_SQ_ori,def_solv
use NAC, only: isNAC

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs, nS, ISphe(nTs), NOrd(nS)
logical(kind=iwp), intent(in) :: Cond
real(kind=wp), intent(in) :: Eps, Sphere(4,nS), Tessera(4,nTs), Q1(2,nTs), Q2(2,nTs), DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3), &
                             DerRad(nS,nAt,3), DerCentr(nS,nAt,3,3)
real(kind=wp), intent(out) :: DerDM(nTs,nTs), Grd(3,nAt)
integer(kind=iwp) :: IAtom, iGo, iTs, IXYZ, jTs
real(kind=wp) :: GeoGrd, Qi, Qj

! Compute the PCM geometric contribution to gradients

Grd(:,:) = Zero
DerDM(:,:) = Zero
do IAtom=1,nAt
  do IXYZ=1,3
    ! Dielectric model
    if (.not. Cond) then
      call Over(IAtom,IXYZ,GeoGrd,nAt,nTs,nS,Eps,Sphere,ISphe,NOrd,Tessera,Q1,DerRad,DerCentr)
    ! Conductor model
    else if (Cond) then
      GeoGrd = Zero
      call DerD(IAtom,IXYZ,Tessera,ISphe,DerDM,DerTes,DerPunt,DerCentr,nTs,nAt,nS)
      if (lSA) then
        !! For SA-CASSCF, q^N*dC/dR*q^N/2 + q^SS*dC/dR*q^SA + q^SS*dC/dR*q^N - q^SA*dC/dR*q^SA/2
        !! Q1(1,:) = Q2(1,:) = q^N
        !! Q1(2,:) = q^SA
        !! Q2(2,:) = q^SS (plus relaxed ASCs)
        if (isNAC) then
          do iTs=1,nTs
            do jTs=1,nTs
              GeoGrd = GeoGrd+Two*(Q1(1,iTs)+Q1(2,iTs))*Q2(2,jTs)*DerDM(iTs,jTs)
            end do
          end do
        else
          do iTs=1,nTs
            do jTs=1,nTs
              if (def_solv==1 .or. def_solv==3) then
                GeoGrd = GeoGrd+(Q1(1,iTs)*Q2(1,jTs) &
                               +Two*Q1(2,iTs)*Q2(2,jTs) &
                               +Two*Q1(1,iTs)*Q2(2,jTs) &
                               -Q1(2,iTs)*Q1(2,jTs))*DerDM(iTs,jTs)
              else if (def_solv==4) then
                GeoGrd = GeoGrd+(Q1(1,iTs)*Q2(1,jTs) &
                                +Two*Q1(2,iTs)*Q2(2,jTs) &
                                +Two*Q1(1,iTs)*Q2(2,jTs) &
                                -Q1(2,iTs)*PCM_SQ_ori(2,jTs))*DerDM(iTs,jTs)
              end if
            end do
          end do
        end if
      else
        !! otherwise, Q1 = Q2
        do iTs=1,nTs
          Qi = Q1(1,iTs)+Q1(2,iTs)
          do jTs=1,nTs
            Qj = Q2(1,jTs)+Q2(2,jTs)
            GeoGrd = GeoGrd+Qi*DerDM(iTs,jTs)*Qj
          end do
        end do
      end if
      GeoGrd = GeoGrd*Eps/(Eps-One)
    end if
    Grd(IXYZ,IAtom) = GeoGrd*Half
  end do
end do

!write (6,'(" Gradient contributions in GeoDer")')
!do iatom = 1, nat
!  write (6,'(i3,3f20.10)') iatom,(Grd(IXYZ,iAtom),ixyz=1,3)
!end do

return

end subroutine GeoDer
