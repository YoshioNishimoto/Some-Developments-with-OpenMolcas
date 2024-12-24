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

subroutine MatPCM(NTs,Eps,Conductor,ISphe,Coor_Sph,Tessera,DMat,SMat,SDMat,TMat,RMat)
! Compute PCM matrix with the formalism in
! M. Cossi, N. Rega, G. Scalmani, V. Barone JCP in press;
! D.M. Chipman JCP 112, 5558 (2000).
! Solvation charges are defined through
! Tq = RV
! where V is the solute electrostatic potential. Here T^-1*R is computed
! and finally returned in DMat.

use PCM_Arrays, only: DiagScale
use Constants, only: Zero, Half, One, Two, Four, Pi
use Definitions, only: wp, iwp
use rctfld_module, only: GauASC, ZetaASC, ISlPar, NS, RSlPar
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer(kind=iwp), intent(in) :: NTs, ISPhe(NTs)
real(kind=wp), intent(in) :: Eps, Coor_Sph(4,*), Tessera(4,NTs)
real(kind=wp), intent(out) :: DMat(NTs,NTs), SMat(NTs,NTs), SDMat(NTs,NTs), TMat(NTs,NTs), RMat(NTs,NTs)
logical(kind=iwp), intent(in) :: Conductor
integer(kind=iwp) :: ITs, JTs, KTs, LI
real(kind=wp) :: EpsFac, Fac, Prod, RIJ, XI, XJ, XNI, YI, YJ, YNI, ZI, ZJ, ZNI
real(kind=wp), parameter :: FPI = Four*Pi, TPI = Two*Pi

integer(kind=iwp) :: LJ,NAT
real(kind=wp) :: alpha,gamma_S,Rhat,RI,RIKJ,Rin,Rsw,Sik,weight,zetai,zetaj,zetap,XIK,YIK,ZIK
real(kind=wp),allocatable :: nTSpSph(:)

if (Conductor) then

  ! Conductor model


  ! S matrix
  EpsFac = Eps/(Eps-One)
  SMat(:,:) = Zero
  ZetaASC = RSlPar(53)
  GauASC = .false.
  if (ZetaASC/=Zero) GauASC = .true.
  if (GauASC) then
    write (6,*) "GauASC = .True."
    ZetaASC = 4.9d+00
!   ZetaASC = 1.0d+00
    weight = sqrt(1.0d+00/NTs) ! we do not consider the weight
!   NAT = ISlPar(42) ! number of atoms
    gamma_S = 1.0d+00 ! degree of switching

! write (6,'(3f20.10)') erf(zetaasc*1.0d-5)/1.0d-05,zetaasc*sqrt(2.0d+00/pi),abs(erf(zetaasc*1.0d-5)/1.0d-05-zetaasc*sqrt(2.0d+00/pi))
! write (6,'(3f20.10)') erf(zetaasc*1.0d-6)/1.0d-06,zetaasc*sqrt(2.0d+00/pi),abs(erf(zetaasc*1.0d-6)/1.0d-06-zetaasc*sqrt(2.0d+00/pi))
! write (6,'(3f20.10)') erf(zetaasc*1.0d-7)/1.0d-07,zetaasc*sqrt(2.0d+00/pi),abs(erf(zetaasc*1.0d-7)/1.0d-07-zetaasc*sqrt(2.0d+00/pi))
! write (6,'(3f20.10)') erf(zetaasc*1.0d-8)/1.0d-08,zetaasc*sqrt(2.0d+00/pi),abs(erf(zetaasc*1.0d-8)/1.0d-08-zetaasc*sqrt(2.0d+00/pi))
! write (6,'(3f20.10)') erf(zetaasc*1.0d-9)/1.0d-09,zetaasc*sqrt(2.0d+00/pi),abs(erf(zetaasc*1.0d-9)/1.0d-09-zetaasc*sqrt(2.0d+00/pi))
! write (6,'(3f20.10)') erf(zetaasc*1.0d-10)/1.0d-10,zetaasc*sqrt(2.0d+00/pi),abs(erf(zetaasc*1.0d-10)/1.0d-10-zetaasc*sqrt(2.0d+00/pi))
! write (6,*) zetaasc/sqrt(2.0d+0*pi)
! write (6,*) 2.0d+00*sqrt(zetaasc/pi)
! write (6,*) "diagscale = ", diagscale

    call mma_allocate(nTSpSph,NS,Label='nTSpSph') ! number of tessera per sphere
    nTSpSph = 0.0d+00
    do ITs = 1, nTS
      LI = ISphe(ITs)
      nTSpSph(LI) = nTSpSph(LI) + 1
    end do
  end if
  do ITs=1,NTs
    XI = Tessera(1,iTs)
    YI = Tessera(2,iTs)
    ZI = Tessera(3,iTs)
    if (GauASC) then
      !! J. Phys. Chem. A 1999, 103, 11069. and possibly J. Chem. Phys. 2010, 133, 244111.
      LI = ISphe(iTs)
      XIK = XI - Coor_Sph(1,LI)
      YIK = YI - Coor_Sph(2,LI)
      ZIK = ZI - Coor_Sph(3,LI)
      Sik = 1.0d+00
      do LI = 1, NS ! number of spheres; some atoms do not have tesserae
      ! RIKJ = SQRT((XIK-Coor_sph(1,LI))**2.0d+00 + (YIK-Coor_Sph(2,LI))**2.0d+00 + (ZIK-Coor_Sph(3,LI))**2.0d+00)
        RIKJ = SQRT((XI-Coor_sph(1,LI))**2.0d+00 + (YI-Coor_Sph(2,LI))**2.0d+00 + (ZI-Coor_Sph(3,LI))**2.0d+00)
        RI = Coor_Sph(4,LI) ! radius of Ith atomic sphere
        Rsw  = gamma_S*sqrt(14.0d+00/nTSpSph(LI))*RI ! Eq.(86)
      ! Rsw  = 1.0d+00*RI
        alpha= 0.5d+00 + RI/Rsw - sqrt((RI/Rsw)**2.0d+00-0.03571428571d+00) ! Eq.(87)
        Rin  = RI - alpha*Rsw ! Eq.(62)
        Rhat = (RIKJ - Rin)/Rsw
!       write (*,*) li,rhat
        if (Rhat < 0.0d+00) then
          Sik = 0.0d+00
          exit
        else if (Rhat > 1.0d+00) then
          cycle
        else
          Sik = Sik*Rhat**3*(10.0d+00 - 15.0d+00*Rhat + 6.0d+00*Rhat*Rhat) ! Eq.(64)
        end if
      end do
      LI = ISphe(iTs)
      zetai = ZetaASC/(weight*Coor_Sph(4,LI))
      sik = sik*4.0d+00
      !! Note that lim_{R->0} S_{ij} = 2*zeta/sqrt(pi)
      SMat(ITs,ITs) = -EpsFac*zetai/(sqrt(TPI)*Sik) ! Eq.(67)
!     SMat(ITs,ITs) = -EpsFac*zetai/(sqrt(PI)*Sik) ! Eq.(67)
! if (its.le.10) then
! write(6,'(1i3,3f20.10)') its,zetai/(sqrt(TPI)*Sik),diagscale*sqrt(FPI/Tessera(4,ITs)),abs(zetai/(sqrt(TPI)*Sik)-diagscale*sqrt(FPI/Tessera(4,ITs)))
! write (6,*) "area = ", tessera(4,its)
! write(6,'(1i3,3f20.10)') its,zetai*2.0d+00/(sqrt(PI)*Sik),diagscale*sqrt(FPI/Tessera(4,ITs)),2.0d+00*abs(zetai/(sqrt(PI)*Sik)-diagscale*sqrt(FPI/Tessera(4,ITs)))
! write(6,'(1i3,3f20.10)') its,zetai/(sqrt(PI)*Sik),diagscale*sqrt(FPI/Tessera(4,ITs)),abs(zetai/(sqrt(PI)*Sik)-diagscale*sqrt(FPI/Tessera(4,ITs)))
! write(6,'(1i3,3f20.10)')&
! its,zetai/(sqrt(PI*0.5d+00)*Sik),sqrt(FPI/Tessera(4,ITs)),abs(zetai/(sqrt(PI*0.5d+00)*Sik)-sqrt(FPI/Tessera(4,ITs)))
! write (6,'(3f20.10)') xi,yi,zi
! end if
!     SMat(ITs,ITs) = -DiagScale*EpsFac*zetai/(sqrt(PI*0.5d+00)*Sik) ! Eq.(3.7)
!     SMat(ITs,ITs) = -DiagScale*EpsFac*sqrt(FPI/Tessera(4,ITs))
    else
      SMat(ITs,ITs) = -DiagScale*EpsFac*sqrt(FPI/Tessera(4,ITs))
    end if
    do JTs=1,ITs-1
      XJ = Tessera(1,jTs)
      YJ = Tessera(2,jTs)
      ZJ = Tessera(3,jTs)
      RIJ = sqrt((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
      if (GauASC) then
        LI = ISphe(iTs)
        zetai = ZetaASC/(weight*Coor_Sph(4,LI))
        LJ = ISphe(jTs)
        zetaj = ZetaASC/(weight*Coor_Sph(4,LJ))
        zetap = zetai*zetaj/sqrt(zetai*zetai+zetaj*zetaj) ! Eq.(68)
        SMat(ITs,JTs) = -EpsFac*erf(zetap*RIJ)/RIJ ! Eq.(66)
      else
        SMat(ITs,JTs) = -EpsFac*One/RIJ
      end if
      SMat(JTs,ITs) = SMat(ITs,JTs)
    end do
  end do

  if (GauASC) then
    call mma_deallocate(nTSpSph)
  end if

  ! Invert S matrix and store it in D

  if (Eps > One) then
    call MatInvert(SMat,nTs)
    DMat(:,:) = SMat(:,:)
  else
    DMat(:,:) = Zero
  end if


else

  ! Dielectric model:

  ! S and D* matrices
  DMat(:,:) = Zero
  do ITs=1,NTs
    XI = Tessera(1,iTs)
    YI = Tessera(2,iTs)
    ZI = Tessera(3,iTs)
    LI = ISphe(ITs)
    XNI = (XI-Coor_Sph(1,LI))/Coor_Sph(4,LI)
    YNI = (YI-Coor_Sph(2,LI))/Coor_Sph(4,LI)
    ZNI = (ZI-Coor_Sph(3,LI))/Coor_Sph(4,LI)
    SMat(ITs,ITs) = DiagScale*sqrt(FPI/Tessera(4,ITs))
    DMat(ITs,ITs) = DMat(ITs,ITs)-TPI/Tessera(4,ITs)
    do JTs=1,NTs
      if (JTs == ITs) cycle
      XJ = Tessera(1,jTs)
      YJ = Tessera(2,jTs)
      ZJ = Tessera(3,jTs)
      RIJ = sqrt((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
      SMat(ITs,JTs) = One/RIJ
      Prod = (XI-XJ)*XNI+(YI-YJ)*YNI+(ZI-ZJ)*ZNI
      DMat(ITs,JTs) = -Prod/RIJ**3
      DMat(JTs,JTs) = DMat(JTs,JTs)-DMat(ITs,JTs)*Tessera(4,ITs)/Tessera(4,JTs)
    end do
  end do

  ! S*A*D matrix
  SDMat(:,:) = Zero
  do ITs=1,NTs
    do JTs=1,NTs
      do KTs=1,NTs
        SDMat(ITs,JTs) = SDMat(ITs,JTs)+SMat(ITs,KTs)*Tessera(4,KTs)*DMat(KTs,JTs)
      end do
    end do
  end do

  ! The charges are defined as
  ! q = T-1 R V,         T = f(e)*S - SAD / 2p

  ! T and R matrices
  Fac = (Eps+One)/(Eps-One)
  TMat(:,:) = Fac*SMat(:,:)-SDMat(:,:)/TPI
  do ITs=1,NTs
    do JTs=1,NTs
      RMat(ITs,JTs) = DMat(JTs,ITs)*Tessera(4,JTs)/TPI
    end do
    RMat(ITs,ITs) = RMat(ITs,ITs)-One
  end do

  ! Invert T matrix

  if (Eps > One) then
    call MatInvert(TMat,nTs)
  else
    TMat(:,:) = Zero
  end if

  ! Form T^-1 * R and store it in D

  call DGEMM_('N','N',nTs,nTs,nTs,One,TMat,nTs,RMat,nTs,Zero,DMat,nTs)

end if

!! Moved from pcm_deriver
do iTs=1,nTs
  do jTs=1,iTs-1
    Fac = Half*(DMat(iTs,jTs)+DMat(jTs,iTs))
    DMat(iTs,jTs) = Fac
    DMat(jTs,iTs) = Fac
  end do
end do

return

end subroutine MatPCM
