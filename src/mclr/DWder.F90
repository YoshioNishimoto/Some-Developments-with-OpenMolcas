!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2024, Yoshio Nishimoto                                 *
!***********************************************************************
  Subroutine DWder(mode,idsym,der1,nder1,der2,nder2,DWOut, &
                   ntash,nRoots,ERASSCF,itoc,LUJOB, &
                   nSym,nAsh,nIsh,nBas,ipCM)

  use definitions, only: iwp,wp
  use stdalloc, only: mma_allocate, mma_deallocate
  use ISRotation, only: weight
  use DWSol, only: DWSol_der

  implicit none

  integer(kind=iwp), intent(in) :: mode,idSym,nder1,nder2
  real(kind=wp), intent(in) :: der1(nder1),der2(nder2)
  real(kind=wp), intent(inout) :: DWOut(nRoots,nRoots)

  integer(kind=iwp),intent(in) :: ntash,nRoots,itoc(*),LUJOB,nSym,nAsh(*),nIsh(*),nBas(*),ipCM(*)
  real(kind=wp),intent(in) :: ERASSCF(1:nRoots)

  real(kind=wp), allocatable :: G1q(:),G2q(:),DERHII(:),DEROMG(:)

  integer(kind=iwp) :: i,iA,ip1,iS,j,jA,jS,jDisk,k,ng1,ng2, &
                       kA,lA,ij1,ij2,kl1,kl2,nnA,ip2
  real(kind=wp) :: OMGDER,rdum(1),Scal
  integer(kind=iwp) :: itri

  itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
!
! derivative of the dynamical weight
! I think this implementation uses the following constraint condition
! w_I^Z * (<I|H|I> - E_I)
! where E is used as a parameter, and w_I^Z is the undetermined multiplier (not weight itself)
!
! For DW-MCSCF, partial derivative of the Lagrangian wrt H_{IJ} can be evaluated with either orbital or CI rotations.
! At present, it is evaluated in ISR_TimesE2.
! For DW solvation, it is only with orbital rotations (here)?
!
  ng1=ntash*(ntash+1)/2
  ng2=ng1*(ng1+1)/2
  Call mma_allocate(G1q,ng1,Label='G1q')
  Call mma_allocate(G2q,ng2,Label='G2q')
  call mma_allocate(DEROMG,nRoots,Label='DEROMG')
  DEROMG = 0.0D+00

  jdisk=itoc(3)
  nnA = ntash
  Do j = 1, nRoots
    Call dDaFile(LUJOB ,2,G1q,ng1,jDisk)
    Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
    Call dDaFile(LUJOB ,2,G2q,Ng2,jDisk)
    Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)

    OMGDER = 0.0d+00
    Do iS=1,nSym
      If (nAsh(iS).gt.0) Then
        jS=iEOr(is-1,iDSym-1)+1
        Do iA=1,nAsh(is)
          Do jA=1,nAsh(js)
            ij1 = itri(ia,ja)
            ip1=nIsh(iS)+iA-1 + nBas(jS)*(nIsh(js)+jA-1)+ipCM(is)
!           write (*,*) dercont(ip1),G1q(max(iA,jA)*(max(iA,jA)-1)/2+min(iA,jA))
            OMGDER = OMGDER + der1(ip1)*G1q(ij1)
            if (mode==1) then
              !! DW-MCSCF
              Do kA=1,nAsh(is)
                Do lA=1,nAsh(js)
                  kl1 = itri(ka,la)
                  !! dimension of rmoaa is itri(nnA**2,nnA**2), but the actual array is similar to G2q
                  !! dimension of G2q  : itri(itri(nnA,nnA),itri(nnA,nnA))
                  scal=1.0d+00
                  if(ij1.ge.kl1 .and. kA.eq.lA) scal=2.0d+00
                  if(ij1.lt.kl1 .and. iA.eq.jA) scal=2.0d+00
                  OMGDER = OMGDER + scal*der2(itri(ij1,kl1))*G2q(itri(ij1,kl1))*0.5d+00
                end do
              end do
            end if
          End Do
        End Do
      End If
    End Do
    DEROMG(j) = OMGDER
!   write (*,*) j,omgder
  End Do

  call mma_allocate(DERHII,nRoots,Label='DERHII')
  DERHII = 0.0D+00

  !! weight is optional
  call DWSol_Der(mode,DEROMG,DERHII,ERASSCF,weight)

  !! this contributions should not be treated as multipliers, I guess
  do i = 1, nRoots
    DWOut(i,i) = DWOut(i,i) + DERHII(i)*0.5d+00
  end do
! write (*,*) "dwout"
! call sqprt(dwout,nroots)

  call mma_deallocate(G1q)
  call mma_deallocate(G2q)
  call mma_deallocate(DERHII)
  call mma_deallocate(DEROMG)

  return

  End Subroutine DWder
