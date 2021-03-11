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
* Copyright (C) 2011, Jose Manuel Hermida Ramon                        *
************************************************************************
      Subroutine OneOverR_Sl(iFil,Ax,Ay,Az,BoMaH,BoMaO,EEDisp,iCNum
     &                  ,Eint,iQ_Atoms,outxyz,Eint_Nuc)
      Implicit Real*8 (a-h,o-z)
*----------------------------------------------------------------------*
*     Jose. Modification of the OneOverR subroutine to include the     *
*     electrostatic penetration of the charge density in the           *
*     electrostatic operator of the Hamiltonian. 2011-05-30            *
*----------------------------------------------------------------------*

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "WrkSpc.fh"

      Parameter (MxK=(MxMltp*(MxMltp**2+6*MxMltp+11)+6)/6)

      Dimension iFil(MxQCen,10)
      Dimension Eint(MxQCen,10),BoMaH(MxAt),BoMaO(MxAt),outxyz(MxQCen,3)
      Dimension Eint_Nuc(MxAt),EintSl(MxK)
      Dimension CTemp(3,MxCen),DTemp(MxCen),DInvTemp(MxCen)
      Logical   lAtom

      EEdisp=0.0d0
*----------------------------------------------------------------------*
* Compute some distances and inverted distances etc. The potential,    *
* the field and etc. and when we already have the numbers, we also do  *
* the dispersion interaction.                                          *
*----------------------------------------------------------------------*
      iCi=(iQ_Atoms*(iQ_Atoms+1))/2
      Do 601, k=1,iCi
        lAtom=.false.
        If(k.le.iQ_atoms) lAtom=.true.
        Gx=outxyz(k,1)+Ax
        Gy=outxyz(k,2)+Ay
        Gz=outxyz(k,3)+Az
        Do 602, j=iCnum+1,nPart
          i=1+(j-1)*nCent
          ip=1+(j-1)*nPol
          Rabx1=Cordst(i,1)-Gx  !Below follows a lot of
          Raby1=Cordst(i,2)-Gy  !distances to and fro.
          Rabz1=Cordst(i,3)-Gz
          Rabx2=Cordst(i+1,1)-Gx
          Raby2=Cordst(i+1,2)-Gy
          Rabz2=Cordst(i+1,3)-Gz
          Rabx3=Cordst(i+2,1)-Gx
          Raby3=Cordst(i+2,2)-Gy
          Rabz3=Cordst(i+2,3)-Gz
          Rabx4=Cordst(i+3,1)-Gx
          Raby4=Cordst(i+3,2)-Gy
          Rabz4=Cordst(i+3,3)-Gz
          Rabx5=Cordst(i+4,1)-Gx
          Raby5=Cordst(i+4,2)-Gy
          Rabz5=Cordst(i+4,3)-Gz
          R21=Rabx1**2+Raby1**2+Rabz1**2
          R22=Rabx2**2+Raby2**2+Rabz2**2
          R23=Rabx3**2+Raby3**2+Rabz3**2
          R24=Rabx4**2+Raby4**2+Rabz4**2
          R25=Rabx5**2+Raby5**2+Rabz5**2
          Rg1=Sqrt(r21)
          Rg2=Sqrt(r22)
          Rg3=Sqrt(r23)
          Rg4=Sqrt(r24)
          Rg5=Sqrt(r25)
          S1i=1/Rg1
          S2i=1/Rg2
          S3i=1/Rg3
          S4i=1/Rg4
          S5i=1/Rg5
          S1e=S1i
          S2e=S2i
          S3e=S3i
          S4e=S4i
          S5e=S5i
          Rab13i=S1e/R21
          Rab23i=S2e/R22
          Rab33i=S3e/R23
          Rab43i=S4e/R24
          Rab53i=S5e/R25
*----------------------------------------------------------------------*
* The dispersion interaction between QM-atoms and solvent is computed, *
* with or without damping. The initial if-clause sees to that only     *
* atom-centers are included, while bonds and virtual centers are       *
* ignored.                                                             *
*----------------------------------------------------------------------*
          If(lAtom) then
            Call DispEnergy(EEDisp,BoMah,BoMaO,rg1,rg2,rg3
     &                     ,Rab13i,Rab23i,Rab33i,k)

          Endif
          Ux1=RabX1*s1i
          Uy1=RabY1*s1i
          Uz1=RabZ1*s1i
          Ux2=RabX2*s2i
          Uy2=RabY2*s2i
          Uz2=RabZ2*s2i
          Ux3=RabX3*s3i
          Uy3=RabY3*s3i
          Uz3=RabZ3*s3i
          Ux4=RabX4*s4i
          Uy4=RabY4*s4i
          Uz4=RabZ4*s4i
          Ux5=RabX5*s5i
          Uy5=RabY5*s5i
          Uz5=RabZ5*s5i
          Rg1=Sqrt(r21)
          Rg2=Sqrt(r22)
          Rg3=Sqrt(r23)
          Rg4=Sqrt(r24)
          Rg5=Sqrt(r25)
          S1i=1/Rg1
          S2i=1/Rg2
          S3i=1/Rg3
          S4i=1/Rg4
          S5i=1/Rg5
*----------------------------------------------------------------------*
* Now we wrap up the electrostatics.                                   *
*----------------------------------------------------------------------*
* Jose. First we check if distance with at least one of the centers of *
* the clasical molecule is inside the Cut-off                          *
*----------------------------------------------------------------------*
          distMin=min(rg1,rg2)
          distMin=min(distMin,rg3)
          distMin=min(distMin,rg4)
          distMin=min(distMin,rg5)
          If(distMin.le.Cut_Elc) then
            CTemp(1,1)=Rabx1
            CTemp(2,1)=Raby1
            CTemp(3,1)=Rabz1
            CTemp(1,2)=Rabx2
            CTemp(2,2)=Raby2
            CTemp(3,2)=Rabz2
            CTemp(1,3)=Rabx3
            CTemp(2,3)=Raby3
            CTemp(3,3)=Rabz3
            CTemp(1,4)=Rabx4
            CTemp(2,4)=Raby4
            CTemp(3,4)=Rabz4
            CTemp(1,5)=Rabx5
            CTemp(2,5)=Raby5
            CTemp(3,5)=Rabz5
            DTemp(1)=Rg1
            DTemp(2)=Rg2
            DTemp(3)=Rg3
            DTemp(4)=Rg4
            DTemp(5)=Rg5
            DInvTemp(1)=S1i
            DInvTemp(2)=S2i
            DInvTemp(3)=S3i
            DInvTemp(4)=S4i
            DInvTemp(5)=S5i
            If(lQuad) then
              nMltTemp=nMlt-1 ! this is done because for Sl_Grad
                              ! charge is L=0, dipole is L=1 and
                              ! quadrupole is L=2. In QmStat they
                              ! are 1, 2 and 3, respectively.

              Call Sl_Grad(nSlSiteC,lMltSlC,CTemp,DTemp,DInvTemp
     &               ,SlExpC,SlFactC,SlPC,nMltTemp,SlExpQ(1,k),DifSlExp
     &               ,EintSl,EintSl_Nuc,lAtom)

*-----------------------------------------------------------------------*
* Change in the order of field gradients beacuse subroutine Sl_Grad
* have the same order than Molcas xx=5, xy=6, xz=7, yy=8, yz=9 and zz=10
* and we need the QmStat order: xx=5, xy=6, yy=7, xz=8, yz=9 and zz=10
*-----------------------------------------------------------------------*
              EintSlTemp=EintSl(7)
              EintSl(7)=EintSl(8)
              EintSl(8)=EintSlTemp

              Do jhr=1,10
                Eint(k,jhr)=Eint(k,jhr)-EintSl(jhr) !Check bellow why
                                               ! it is a subtraction
                                               ! and not a sum
              End do
              If(lAtom) then
                Eint_Nuc(k)=Eint_Nuc(k)-EintSl_Nuc
              End if
              go to 3466
            Else
              nMltTemp=nMlt-1 ! this is done because for Sl_Grad
                              ! charge is L=0, dipole is L=1 and
                              ! quadrupole is L=2. In QmStat they
                              ! are 1, 2 and 3, respectively.

              ijhr=min(nMltTemp,1)
              Call Sl_Grad(nSlSiteC,lMltSlC,CTemp,DTemp,DInvTemp
     &                 ,SlExpC,SlFactC,SlPC,ijhr,SlExpQ(1,k),DifSlExp
     &                 ,EintSl,EintSl_Nuc,lAtom)

              Do jhr=1,4
                Eint(k,jhr)=Eint(k,jhr)-EintSl(jhr) ! Check bellow why
                                               ! it is a subtraction
                                               ! and not a sum
              End do
              If(lAtom) then
                Eint_Nuc(k)=Eint_Nuc(k)-EintSl_Nuc
              End if
              go to 3456
            Endif
          Endif
*----------------------------------------------------------------------*
* The Eint(k,1) term will below turn into the interaction between
* charges on water and the MME-charges on the QM-molecule.
* Plus the interaction between the charge densities in Water and the
* charge densities in QM-Molecule when the Penetration is evaluated
*----------------------------------------------------------------------*
          Eint(k,1)=-Qsta(1)*S2e-Qsta(2)*S3e-Qsta(3)*S4e
     &             -Qsta(4)*S5e+Eint(k,1)
          If(lAtom) then
            Eint_Nuc(k)=-Qsta(1)*S2e-Qsta(2)*S3e-Qsta(3)*S4e
     &             -Qsta(4)*S5e+Eint_Nuc(k)
          Endif

         !These three terms will below turn into the interaction
         !between water charges and the MME-dipoles on the QM-mol.
         !Change sign of charge, change sign of vector and then we
         !should also change sign since when a dipole interacts with
         !a field we have a minus sign, but this minus sign we have
         !omitted in hel; therefore this calculation gives the right
         !number eventually.
          Eint(k,2)=-Qsta(1)*RabX2*Rab23i-Qsta(2)*RabX3*Rab33i
     &-Qsta(3)*RabX4*Rab43i-Qsta(4)*RabX5*Rab53i+Eint(k,2)
          Eint(k,3)=-Qsta(1)*RabY2*Rab23i-Qsta(2)*RabY3*Rab33i
     &-Qsta(3)*RabY4*Rab43i-Qsta(4)*RabY5*Rab53i+Eint(k,3)
          Eint(k,4)=-Qsta(1)*RabZ2*Rab23i-Qsta(2)*RabZ3*Rab33i
     &-Qsta(3)*RabZ4*Rab43i-Qsta(4)*RabZ5*Rab53i+Eint(k,4)
          !And here it is the MME-quadrupoles that are prepared.
          !Change sign of charges, change sign two times of the
          !vector (in effect, zero times then) and then in the
          !energy expression for the interaction between the field
          !vector from a charge and a quarupole there is a plus
          !sign, so a minus is the right sign below.

3456      Continue

          Eint(k,5)=Eint(k,5)-Qsta(1)*Ux2**2*Rab23i-Qsta(2)*Ux3**2
     &*Rab33i-Qsta(3)*Ux4**2*Rab43i-Qsta(4)*Ux5**2*Rab53i
          Eint(k,7)=Eint(k,7)-Qsta(1)*Uy2**2*Rab23i-Qsta(2)*Uy3**2
     &*Rab33i-Qsta(3)*Uy4**2*Rab43i-Qsta(4)*Uy5**2*Rab53i
          Eint(k,10)=Eint(k,10)-Qsta(1)*Uz2**2*Rab23i-Qsta(2)*Uz3**2
     &*Rab33i-Qsta(3)*Uz4**2*Rab43i-Qsta(4)*Uz5**2*Rab53i
          Eint(k,6)=Eint(k,6)-Qsta(1)*Ux2*Uy2*Rab23i-Qsta(2)*Ux3
     &*Uy3*Rab33i-Qsta(3)*Ux4*Rab43i*Uy4-Qsta(4)*Ux5*Rab53i*Uy5
          Eint(k,8)=Eint(k,8)-Qsta(1)*Ux2*Uz2*Rab23i-Qsta(2)*Ux3
     &*Uz3*Rab33i-Qsta(3)*Ux4*Rab43i*Uz4-Qsta(4)*Ux5*Rab53i*Uz5
          Eint(k,9)=Eint(k,9)-Qsta(1)*Uz2*Uy2*Rab23i-Qsta(2)*Uz3
     &*Uy3*Rab33i-Qsta(3)*Uz4*Rab43i*Uy4-Qsta(4)*Uz5*Rab53i*Uy5

3466      Continue

*----------------------------------------------------------------------*
* And now a whole lot of grad(1/r) and higher...                       *
*----------------------------------------------------------------------*
      !Unipoles.
          Work(iFil(k,1)-1+ip)=Rabx1*Rab13i
          Work(iFil(k,1)-1+ip+1)=Rabx2*Rab23i
          Work(iFil(k,1)-1+ip+2)=Rabx3*Rab33i
          Work(iFil(k,1)-1+nPart*nPol+ip)=Raby1*Rab13i
          Work(iFil(k,1)-1+nPart*nPol+ip+1)=Raby2*Rab23i
          Work(iFil(k,1)-1+nPart*nPol+ip+2)=Raby3*Rab33i
          Work(iFil(k,1)-1+2*nPart*nPol+ip)=Rabz1*Rab13i
          Work(iFil(k,1)-1+2*nPart*nPol+ip+1)=Rabz2*Rab23i
          Work(iFil(k,1)-1+2*nPart*nPol+ip+2)=Rabz3*Rab33i
      !Dipole -- x-component.
          Work(iFil(k,2)-1+ip)=-(1-3*Ux1**2)*Rab13i
          Work(iFil(k,2)-1+ip+1)=-(1-3*Ux2**2)*Rab23i
          Work(iFil(k,2)-1+ip+2)=-(1-3*Ux3**2)*Rab33i
          Work(iFil(k,2)-1+nPart*nPol+ip)=Uy1*Ux1*Rab13i*3
          Work(iFil(k,2)-1+nPart*nPol+ip+1)=Uy2*Ux2*Rab23i*3
          Work(iFil(k,2)-1+nPart*nPol+ip+2)=Uy3*Ux3*Rab33i*3
          Work(iFil(k,2)-1+2*nPart*nPol+ip)=Uz1*Ux1*Rab13i*3
          Work(iFil(k,2)-1+2*nPart*nPol+ip+1)=Uz2*Ux2*Rab23i*3
          Work(iFil(k,2)-1+2*nPart*nPol+ip+2)=Uz3*Ux3*Rab33i*3
      !Dipole -- y-component.
          Work(iFil(k,3)-1+ip)=Uy1*Ux1*Rab13i*3
          Work(iFil(k,3)-1+ip+1)=Uy2*Ux2*Rab23i*3
          Work(iFil(k,3)-1+ip+2)=Uy3*Ux3*Rab33i*3
          Work(iFil(k,3)-1+nPart*nPol+ip)=-(1-3*Uy1**2)*Rab13i
          Work(iFil(k,3)-1+nPart*nPol+ip+1)=-(1-3*Uy2**2)*Rab23i
          Work(iFil(k,3)-1+nPart*nPol+ip+2)=-(1-3*Uy3**2)*Rab33i
          Work(iFil(k,3)-1+2*nPart*nPol+ip)=Uz1*Uy1*Rab13i*3
          Work(iFil(k,3)-1+2*nPart*nPol+ip+1)=Uz2*Uy2*Rab23i*3
          Work(iFil(k,3)-1+2*nPart*nPol+ip+2)=Uz3*Uy3*Rab33i*3
      !Dipole -- z-component.
          Work(iFil(k,4)-1+ip)=Uz1*Ux1*Rab13i*3
          Work(iFil(k,4)-1+ip+1)=Uz2*Ux2*Rab23i*3
          Work(iFil(k,4)-1+ip+2)=Uz3*Ux3*Rab33i*3
          Work(iFil(k,4)-1+nPart*nPol+ip)=Uz1*Uy1*Rab13i*3
          Work(iFil(k,4)-1+nPart*nPol+ip+1)=Uz2*Uy2*Rab23i*3
          Work(iFil(k,4)-1+nPart*nPol+ip+2)=Uz3*Uy3*Rab33i*3
          Work(iFil(k,4)-1+2*nPart*nPol+ip)=-(1-3*Uz1**2)*Rab13i
          Work(iFil(k,4)-1+2*nPart*nPol+ip+1)=-(1-3*Uz2**2)*Rab23i
          Work(iFil(k,4)-1+2*nPart*nPol+ip+2)=-(1-3*Uz3**2)*Rab33i
      !Quadrupole -- xx-component.
          Work(iFil(k,5)-1+ip)=(5*Ux1*(Ux1*Ux1-.4))*Rab13i*S1e
          Work(iFil(k,5)-1+ip+1)=(5*Ux2*(Ux2*Ux2-.4))*Rab23i*S2e
          Work(iFil(k,5)-1+ip+2)=(5*Ux3*(Ux3*Ux3-.4))*Rab33i*S3e
          Work(iFil(k,5)-1+nPart*nPol+ip)=5*Uy1*Ux1*Ux1*Rab13i*S1e
          Work(iFil(k,5)-1+nPart*nPol+ip+1)=5*Uy2*Ux2*Ux2*Rab23i*S2e
          Work(iFil(k,5)-1+nPart*nPol+ip+2)=5*Uy3*Ux3*Ux3*Rab33i*S3e
          Work(iFil(k,5)-1+2*nPart*nPol+ip)=5*Uz1*Ux1*Ux1*Rab13i*S1e
          Work(iFil(k,5)-1+2*nPart*nPol+ip+1)=5*Uz2*Ux2*Ux2*Rab23i*S2e
          Work(iFil(k,5)-1+2*nPart*nPol+ip+2)=5*Uz3*Ux3*Ux3*Rab33i*S3e
      !Quadrupole -- yy-component.
          Work(iFil(k,7)-1+ip)=5*Uy1*Uy1*Ux1*Rab13i*S1e
          Work(iFil(k,7)-1+ip+1)=5*Uy2*Uy2*Ux2*Rab23i*S2e
          Work(iFil(k,7)-1+ip+2)=5*Uy3*Uy3*Ux3*Rab33i*S3e
          Work(iFil(k,7)-1+nPart*nPol+ip)=(5*Uy1*(Uy1*Uy1-.4))
     &                                   *Rab13i*S1e
          Work(iFil(k,7)-1+nPart*nPol+ip+1)=(5*Uy2*(Uy2*Uy2-.4))
     &                                   *Rab23i*S2e
          Work(iFil(k,7)-1+nPart*nPol+ip+2)=(5*Uy3*(Uy3*Uy3-.4))
     &                                   *Rab33i*S3e
          Work(iFil(k,7)-1+2*nPart*nPol+ip)=5*Uz1*Uy1*Uy1*Rab13i*S1e
          Work(iFil(k,7)-1+2*nPart*nPol+ip+1)=5*Uz2*Uy2*Uy2*Rab23i*S2e
          Work(iFil(k,7)-1+2*nPart*nPol+ip+2)=5*Uz3*Uy3*Uy3*Rab33i*S3e
      !Quadrupole -- zz-component.
          Work(iFil(k,10)-1+ip)=5*Uz1*Uz1*Ux1*Rab13i*S1e
          Work(iFil(k,10)-1+ip+1)=5*Uz2*Uz2*Ux2*Rab23i*S2e
          Work(iFil(k,10)-1+ip+2)=5*Uz3*Uz3*Ux3*Rab33i*S3e
          Work(iFil(k,10)-1+nPart*nPol+ip)=5*Uz1*Uz1*Uy1*Rab13i*S1e
          Work(iFil(k,10)-1+nPart*nPol+ip+1)=5*Uz2*Uz2*Uy2*Rab23i*S2e
          Work(iFil(k,10)-1+nPart*nPol+ip+2)=5*Uz3*Uz3*Uy3*Rab33i*S3e
          Work(iFil(k,10)-1+2*nPart*nPol+ip)=(5*Uz1*(Uz1*Uz1-.4))
     &                                      *Rab13i*S1e
          Work(iFil(k,10)-1+2*nPart*nPol+ip+1)=(5*Uz2*(Uz2*Uz2-.4))
     &                                      *Rab23i*S2e
          Work(iFil(k,10)-1+2*nPart*nPol+ip+2)=(5*Uz3*(Uz3*Uz3-.4))
     &                                      *Rab33i*S3e
      !Quadrupole -- xy-component.
          Work(iFil(k,6)-1+ip)=(5*Uy1*(Ux1*Ux1-.2))*Rab13i*S1e
          Work(iFil(k,6)-1+ip+1)=(5*Uy2*(Ux2*Ux2-.2))*Rab23i*S2e
          Work(iFil(k,6)-1+ip+2)=(5*Uy3*(Ux3*Ux3-.2))*Rab33i*S3e
          Work(iFil(k,6)-1+nPart*nPol+ip)=(5*Ux1*(Uy1*Uy1-.2))
     &                                   *Rab13i*S1e
          Work(iFil(k,6)-1+nPart*nPol+ip+1)=(5*Ux2*(Uy2*Uy2-.2))
     &                                   *Rab23i*S2e
          Work(iFil(k,6)-1+nPart*nPol+ip+2)=(5*Ux3*(Uy3*Uy3-.2))
     &                                   *Rab33i*S3e
          Work(iFil(k,6)-1+2*nPart*nPol+ip)=5*Uz1*Uy1*Ux1*Rab13i*S1e
          Work(iFil(k,6)-1+2*nPart*nPol+ip+1)=5*Uz2*Uy2*Ux2*Rab23i*S2e
          Work(iFil(k,6)-1+2*nPart*nPol+ip+2)=5*Uz3*Uy3*Ux3*Rab33i*S3e
      !Quadrupole -- xz-component.
          Work(iFil(k,8)-1+ip)=(5*Uz1*(Ux1*Ux1-.2))*Rab13i*S1e
          Work(iFil(k,8)-1+ip+1)=(5*Uz2*(Ux2*Ux2-.2))*Rab23i*S2e
          Work(iFil(k,8)-1+ip+2)=(5*Uz3*(Ux3*Ux3-.2))*Rab33i*S3e
          Work(iFil(k,8)-1+nPart*nPol+ip)=5*Uz1*Uy1*Ux1*Rab13i*S1e
          Work(iFil(k,8)-1+nPart*nPol+ip+1)=5*Uz2*Uy2*Ux2*Rab23i*S2e
          Work(iFil(k,8)-1+nPart*nPol+ip+2)=5*Uz3*Uy3*Ux3*Rab33i*S3e
          Work(iFil(k,8)-1+2*nPart*nPol+ip)=(5*Ux1*(Uz1*Uz1-.2))
     &                                     *Rab13i*S1e
          Work(iFil(k,8)-1+2*nPart*nPol+ip+1)=(5*Ux2*(Uz2*Uz2-.2))
     &                                     *Rab23i*S2e
          Work(iFil(k,8)-1+2*nPart*nPol+ip+2)=(5*Ux3*(Uz3*Uz3-.2))
     &                                     *Rab33i*S3e
      !Quadrupole -- yz-component.
          Work(iFil(k,9)-1+ip)=5*Uz1*Uy1*Ux1*Rab13i*S1e
          Work(iFil(k,9)-1+ip+1)=5*Uz2*Uy2*Ux2*Rab23i*S2e
          Work(iFil(k,9)-1+ip+2)=5*Uz3*Uy3*Ux3*Rab33i*S3e
          Work(iFil(k,9)-1+nPart*nPol+ip)=(5*Uz1*(Uy1*Uy1-.2))
     &                                   *Rab13i*S1e
          Work(iFil(k,9)-1+nPart*nPol+ip+1)=(5*Uz2*(Uy2*Uy2-.2))
     &                                   *Rab23i*S2e
          Work(iFil(k,9)-1+nPart*nPol+ip+2)=(5*Uz3*(Uy3*Uy3-.2))
     &                                   *Rab33i*S3e
          Work(iFil(k,9)-1+2*nPart*nPol+ip)=(5*Uy1*(Uz1*Uz1-.2))
     &                                   *Rab13i*S1e
          Work(iFil(k,9)-1+2*nPart*nPol+ip+1)=(5*Uy2*(Uz2*Uz2-.2))
     &                                   *Rab23i*S2e
          Work(iFil(k,9)-1+2*nPart*nPol+ip+2)=(5*Uy3*(Uz3*Uz3-.2))
     &                                   *Rab33i*S3e
*----------------------------------------------------------------------*
* If damping of the field is requested, then do it.                    *
*----------------------------------------------------------------------*
          If(FieldDamp) then
            Do 620, ijhr=1,10
              Do 621, jjhr=0,2
                Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip)=
     &          Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip)
     &          *(1-exp(CAFieldG*rg1))**CFexp
                Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip+1)=
     &          Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip+1)
     &          *(1-exp(CBFieldG*rg2))**CFexp
                Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip+2)=
     &          Work(iFil(k,ijhr)-1+jjhr*nPart*nPol+ip+2)
     &          *(1-exp(CBFieldG*rg3))**CFexp
621           Continue
620         Continue
          Endif
602     Continue
601   Continue

      Return
      End
