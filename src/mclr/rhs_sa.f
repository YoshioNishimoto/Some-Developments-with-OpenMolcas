************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       Subroutine rhs_sa(Fock,SLag,ipS2)
       use Arrays, only: Int1, CMO, fimo,famo, g1t,g2t
       use Data_Structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
       use ipPage, only: W
       use rctfld_module, only: lRf
       use PCM_grad, only: PCM_grad_dens,PCM_grad_dens2,DSSMO,DSSAO,
     *                     PCMSSMO,dscfmo,dscfao,pcmscfmo,potnuc_pcm,
     *                     preppcm,PrepPCM2,PCMSSAO,def_solv,DZMO,
     *                     PCM_grad_CLag,ISRot,do_RF
       use ISRotation, only: InvSCF,ISR,InvEne
       Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
#include "Files_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"

       Real*8 Fock(*)
       real*8, optional :: SLag(*)
       integer, optional :: ipS2
       Dimension rdum(1)
       Real*8, Allocatable:: T(:), F(:), G1q(:), G2q(:), G1r(:), G2r(:)
       logical :: first, Dff, NonEq
       Real*8, Allocatable:: h1(:), TwoHam(:), D(:)

*                                                                      *
************************************************************************
*                                                                      *
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
C      calL compute_lag
C      call abend
*                                                                      *
************************************************************************
*                                                                      *
*
       if(doDMRG)then ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),ntash,0)
       end if

       ng1=itri(ntash,ntash)
       ng2=itri(ng1,ng1)
*
       Call mma_allocate(T,ndens2,Label='T')
       Call mma_allocate(F,ndens2,Label='F')
       Call mma_allocate(G1q,ng1,Label='G1q')
       Call mma_allocate(G2q,ng2,Label='G2q')
       Call mma_allocate(G1r,ntash**2,Label='G1r')
       Call mma_allocate(G2r,itri(ntash**2,ntash**2),Label='G2r')
*
**     Pick up densities from JobIph file
*
       iR=iroot(istate)
C      write (6,*) "iR = ", ir, istate
       jdisk=itoc(3)
       Do ii=1,iR-1
         Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
         Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
         Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)
       End Do
       Call dDaFile(LUJOB ,2,G1q,ng1,jDisk)
       Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
       Call dDaFile(LUJOB ,2,G2q,Ng2,jDisk)
       Call dDaFile(LUJOB ,0,rdum,Ng2,jDisk)

       !! test with the state-averaged density
C      call dcopy_(ng1,g1t,1,g1q,1)
C      call dcopy_(ng2,g2t,1,g2q,1)
*
       !! Add SLag (rotations of states) contributions from the partial
       !! derivative of the CASPT2 energy. G1q and G2q are modified.
       !! The modified density will be used in out_pt2.f and ptrans_sa.f
       !! If the target energy is not invariant, rotations that are in
       !! the SLag array is evaluated here (even if everything zero).
C      If ((PT2.and.nRoots.gt.1) .or.
C    *     (lRF.and.def_solv.ne.3.and.InvSCF)) then
C      If (.not.InvEne .and. InvSCF) then
       If (.not.InvEne) then
         write (6,*) "state-specific before slag"
         Do iB=1,ntash
          Do jB=1,ntash
          G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
          End Do
         End Do
         call sqprt(g1r,nash(1))

        !g1q = 0.0d+00
        !g2q = 0.0d+00
C        call daxpy_(nRoots**2,1.0d+00,ISR%p,1,SLag,1)
         Call PT2_SLag()
       End If
*
       Call Put_dArray('P2mo',G2q,ng2)
       Call Put_dArray('D1mo',G1q,ng1)
*
       Do iB=1,ntash
        Do jB=1,ntash
        G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
        End Do
       End Do
       write (6,*) "state-specific G1r after slag rotation"
       call sqprt(g1r,nash(1))

C      if (lRF.and.def_solv.ne.3.and.InvSCF) then
C      if (lRF .and. .not.InvEne .and. InvSCF) then
C      If (lRF .and. ((PT2.and.nRoots.gt.1) .or.
C    *     (def_solv.ne.3.and.InvSCF))) then
       If (do_RF .and. ((PT2.and.nRoots.gt.1)
     *             .or. (.not.InvEne .and. InvSCF))) then
C        write (6,*) "doing PCM-relevant after G1r rotation"
         !! Update state-specific quantities using rotated G1r
         call dcopy_(ntash**2,g1r,1,dssmo,1)
C        write (*,*) "dssmo before"
C        call sqprt(dssmo,ntash)
C        write (*,*) "dssao before"
C        call sqprt(dssao,nbas(1))
C        write (*,*) "pcmssao before"
C        call sqprt(pcmssao,nbas(1))
C        write (*,*) "pcmssmo before"
C        call sqprt(pcmssmo,nbas(1))
         call PrepPCM2(1,dssmo,dssao,pcmssao,pcmssmo)
         ! compute CLag again with the rotated SS density
         ! There is no need to construct the derivative twice, but
         ! we anyway have to evaluate the derivatives of energy and
         ! Lagrangian separately, so there is no much benefit to change
         Call PCM_grad_CLag(2,ipCI,ipS2,SLag)


C        do i = 1, 5
C          write (*,*) "slag = ", slag(1:ntash**2)
C        write (*,*) "before"
C        call sqprt(g1q,ntash)
C        Call PT2_SLag()
C        write (*,*) "after"
C        call sqprt(g1q,ntash)
C      Do iB=1,ntash
C       Do jB=1,ntash
C       G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
C       End Do
C      End Do
C        call dcopy_(ntash**2,g1r,1,dssmo,1)
C        call PrepPCM2(dssmo,dssao,pcmssao,pcmssmo)
C        Call PCM_grad_CLag(1,ipCI,ipS2,SLag)
C        end do

C        write (*,*) "dssmo after"
C        call sqprt(dssmo,ntash)
C        write (*,*) "dssao after"
C        call sqprt(dssao,nbas(1))
C        write (*,*) "pcmssao after"
C        call sqprt(pcmssao,nbas(1))
C        write (*,*) "pcmssmo after"
C        call sqprt(pcmssmo,nbas(1))
       end if

       if (.false.) then
       !! tentative
       if (lrf) call dcopy_(ntash**2,g1r,1,dssmo,1)
       e1 = 0.0d+00
       e1fimo = 0.0d+00
       write (6,*) "fimosize = ", size(fimo)
       write (6,*) "int1"
       call sqprt(int1,nbas(1))
       write (6,*) "fimo"
       call sqprt(fimo,nbas(1))
       write (6,*) "nbas(1) = ", nbas(1)
       e12 = 0.0d+00
       do i = 1, nish(1)
         nseq = i*(i-1)/2 + i
C        e1 = e1 + 2.0d+00*int1(nseq)
         e1 = e1 + 2.0d+00*int1(i+nbas(1)*(i-1))
         e1fimo = e1fimo + 2.0d+00*fimo(i+nbas(1)*(i-1))
         if (lrf) then
           e1fimo = e1fimo - 2.0d+00*pcmscfmo(i+nbas(1)*(i-1),3)
         end if
       end do
       write (6,'("one-electron energy1= ",f20.10)') e1
       write (6,'("FIMO         energy1= ",f20.10)') e1fimo
       e12 = e1 + (e1fimo-e1)*0.5d+00
       write (6,'("e12 = rcore/2       = ",f20.10)') e12
       e1fimoact = e1fimo
       e1sav = e1
       erfx = 0.0d+00
       ecorr = 0.0d+00
       do i = 1, nash(1)
         do j = 1, nash(1)
           nseq = max(nish(1)+i,nish(1)+j)
     *     *(max(nish(1)+i,nish(1)+j)-1)/2+min(nish(1)+i,nish(1)+j)
C          e1 = e1 + int1(nseq)*g1r(i+ntash*(j-1))
           e1 = e1 + int1(nish(1)+i+nbas(1)*(nish(1)+j-1))
     *             *g1r(i+ntash*(j-1))
           e1fimo = e1fimo
     *     + g1r(i+ntash*(j-1))*fimo(nish(1)+i+nbas(1)*(nish(1)+j-1))
       if (lrf) then
       erfx = erfx
     *   + dscfmo(i,j)*pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
      ecorr =  ecorr
     *   + 0.5d+00*(dscfmo(i,j)-dssmo(i,j))
     *     *pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
        end if
         end do
       end do
       e1fimoact = e1fimo-e1fimoact
       enuc = 28.65127677d+00
       write (6,'("nuclear repulsion   = ",f20.10)') enuc
       write (6,'("one-electron energy = ",f20.10)') e1
       write (6,'("FIMO         energy = ",f20.10)') e1fimo
       write (6,'("one-2        energy = ",f20.10)') e12 + enuc
       write (6,'("one-3        energy = ",f20.10)')
     *  e12 + enuc + e1fimoact
       write (6,'("inact-act    energy = ",f20.10)') e1fimoact
       write (6,'("Erfx                = ",f20.10)') -0.5d+00*erfx
       write (6,'("EMY                 = ",f20.10)')
     *  e12 + enuc + potnuc_pcm  -0.5d+00*erfx
      end if

       Do iB=1,ntash
        Do jB=1,ntash
         iDij=iTri(ib,jB)
         iRij=jb+(ib-1)*ntash
         Do kB=1,ntash
          Do lB=1,ntash
           iDkl=iTri(kB,lB)
           iRkl=lb+(kb-1)*ntash
           fact=One
           if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
           if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
           iijkl=itri(iDij,iDkl)
           iRijkl=itri(iRij,iRkl)
           G2r(iRijkl)=Fact*G2q(iijkl)
          End Do
         End Do
        End Do
       End Do

C      write (*,*) "g2r"
C      do i = 1, itri(ntash**2,ntash**2)
C      write (*,'(i4,f20.10)') i,g2r(i)
C      end do
C      call abend

      if (do_RF) then
      write (6,*) "calling compute_energy"
C     if (lRF) call preppcm(iSpin)
      call compute_energy(eee,G1r,G2r,rdum)
      end if

      if (.false.) then
      write (6,*) "computing the energy"
      call compute_energy(eee,G1r,G2r,rdum)
      write (6,*) "eee = ", eee
      call abend

       e2 = 0.0d+00
       Do iB=1,ntash
        Do jB=1,ntash
         iDij=iTri(ib,jB)
         iRij=jb+(ib-1)*ntash
         Call Coul(1,1,1,1,nish(1)+iB,nish(1)+jB,T,F)
C        write (6,*) "active ib,jb= ", ib,jb
C        call sqprt(t,nbas(1))
         Do kB=1,ntash
          Do lB=1,ntash
           iDkl=iTri(kB,lB)
           iRkl=lb+(kb-1)*ntash
           fact=One
           if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
           if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
           iijkl=itri(iDij,iDkl)
           iRijkl=itri(iRij,iRkl)
           G2r(iRijkl)=Fact*G2q(iijkl)
C          e2 = e2 + t(nish(1)+kb-1+nbas(1)*(nish(1)+lb-1))*g2r(irijkl)
C          iDik = itri(ib,kb)
C          iDjl = itri(jb,lb)
C          iijkl=itri(iDik,iDjl)
        !e2 = e2 + t(nish(1)+kb+nbas(1)*(nish(1)+lb-1))*g2q(iijkl)
         e2 = e2 + t(nish(1)+kb+nbas(1)*(nish(1)+lb-1))*g2r(irijkl)
          End Do
         End Do
        End Do
       End Do
       write (6,'("two-electron energy = ",f20.10)') e2
       write (6,*)
C      write (6,'("total (?) energy    = ",f20.10)') enuc + e1fimo + e2
       write (6,'("total (?) energy    = ",f20.10)') e12+enuc+e1fimoact
     *   + e2*0.5d+00
      etot =  e12+enuc+potnuc_pcm-0.5d+00*erfx+e1fimoact + e2*0.5d+00
       write (6,'("total (?) energy    = ",f20.10)') etot
       eint = e12 + enuc + e1fimoact + e2*0.5d+00
       write (6,'("delta variational   = ",f20.10)') ecorr
       write (6,'("total (?) energy+cor= ",f20.10)') etot + ecorr
       write (6,'("reference energy    = ",f20.10)') -112.14570255169
       write (6,'("electronic energy   = ",f20.10)')
     *  -112.14570255169-enuc
       write (6,*)
C      write (6,'("core energy fmat    = ",f20.10)') -103.360940266019
C      write (6,'("inactive-active     = ",f20.10)') -14.6797467804585
C      write (6,'("CAS energy (core+int) ",f20.10)') -118.04068706478
       write (6,'("core energy fmat    = ",f20.10)') -103.369714485012
       write (6,'("inactive-active     = ",f20.10)') -14.6807793028748
       write (6,'("CAS energy (core+int) ",f20.10)') -118.050493787886
       end if
*
       if(doDMRG)then ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),nna,0)
       end if

      !! PCM
      if (lRf .and. .false.) then
        nh1 = 0 ! number of AOs
        do iSym = 1, nSym
          nh1 = nh1 + nBas(iSym)*(nBas(iSym)+1)/2
        end do
        RepNuc=Zero ! not used actually
        First = .true.
        Dff = .false.
        NonEq = .false.
        iCharge=Int(Tot_Charge)
*
        Call mma_allocate(h1,ndens2,Label='h1')
        Call mma_allocate(TwoHam,ndens2,Label='TwoHam')
        Call mma_allocate(D,ndens2,Label='D')
*
        TwoHam(:) = Zero
        h1(:) = Zero
C       call PCM_grad_dens(2) !! DSSMO
        call PCM_grad_dens2(1,DSSMO,DSSAO)
        call PCM_grad_dens(1) !! DSCFMO
        call PCM_grad_dens2(1,DScfMO,DScfAO)
        !! use the density
      ! call SADENS()
      ! !! transform MO to AO
      ! h1(:) = Zero
      ! Do iS=1,nSym
      !   iOrb = nOrb(iS)
      !   Do iB=1,nFro(is)+nIsh(is)
      !     h1(ipmat(iS,iS)+iB-1+iOrb*(iB-1)) = Two
      !   End Do
      !
      !   Do iB=1,nAsh(iS)
      !     Do jB=1,nAsh(iS)
      !       iiB=nA(iS)+ib
      !       ijB=nA(iS)+jb
      !       iij=iTri(iib,ijb)
      !       iiB=nFro(iS)+nIsh(iS)+ib
      !       ijB=nFro(iS)+nIsh(iS)+jb
      !       h1(ipmat(iS,iS)+iiB-1+iOrb*(ijB-1)) = G1q(iij)
      !     End Do
      !   End Do
      ! End Do
C     ! write (6,*) "D in MO"
C     ! call sqprt(h1,nbas(1))
      ! call tcmo(h1,1,-2)
      ! call fold(nSym,nBas,h1,d)
        call fold(nSym,nBas,DSSAO,d)
      ! h1(:) = Zero
        !! D(AO) -> V(AO)
        Call DrvRF(h1,TwoHam,D,RepNuc,nh1,First,Dff,NonEq,
     &             iCharge)
        !! V(AO) -> V(MO)
        !! First, add nuclear contributions to electrostatic
        call daxpy_(nh1,1.0d+00,h1,1,TwoHam,1)
C       call tcmo(TwoHam,1,1)
        !! temporarily
        call square(TwoHam,PCMSSMO,1,nBas(1),nBas(1))
        call tcmo(PCMSSMO,1,1)
       !write (6,*) "pcmssmo"
       !call sqprt(pcmssmo,nbas(1))
        !! Add to the initial Lagrangian

C       call PCM_grad_dens(2) !! DSSMO
C       call PCM_grad_dens2(1,DSSMO,DSSAO)
C       call PCM_grad_dens(1) !! DSCFMO
C       call PCM_grad_dens2(1,DScfMO,DScfAO)
*
        Call mma_deallocate(h1)
        Call mma_deallocate(TwoHam)
        Call mma_deallocate(D)
      end if

       Call FockGen(One,G1r,G2r,T,Fock,1)
*       Do iS=1,nsym
*        Call RecPrt(' ',' ',fock(ipMat(is,is)),nbas(is),nbas(is))
*       End Do

      If (.not.debug) Then !yma debug ??
       renergy=Zero
       Do ii=1,nsym
        Do jj=1,nbas(ii)
         renergy = renergy + T(ipmat(ii,ii)+jj-1+nbas(ii)*(jj-1))
        End DO
       End DO

      rcorei=Zero
      rcorea=Zero
      Do iS=1,nSym
       Do iB=1,nIsh(is)
       rcorei=rcorei+Two*Int1(ipCM(is)+nOrb(iS)*(ib-1)+ib-1)
       End Do

       Do iB=1,nAsh(iS)
        Do jB=1,nAsh(iS)
         iiB=nA(iS)+ib
         ijB=nA(iS)+jb
         iij=iTri(iib,ijb)
         iiB=nIsh(iS)+ib
         ijB=nIsh(iS)+jb
         rcorea=rcorea+G1q(iij)*Int1(ipCM(is)+nOrb(is)*(iib-1)+ijB-1)
        End Do
       End Do
      End Do
!      rcore=rCorei+rcoreA
!      write(6,*) 'In rhs_sa'
!      Write(6,*) 'Checking energy',0.5d0*renergy+potnuc+half*rcore !yma
!      Write(6,*) 'Checking energy',0.5d0*renergy,potnuc,rcore      !yma
!      write(6,*)
      End if
!      Do iS=1,nsym
!       Call RecPrt(' ',' ',fock(ipMat(is,is)),nbas(is),nbas(is))
!      End Do
*
*
       Call mma_deallocate(G1q)
       Call mma_deallocate(G2q)
*
       Call TCMO(T,1,-2)
       ijb=0
       Do is=1,nsym
        Do ib=1,nbas(is)
         Do jb=1,ib-1
          ijb=ijb+1
          F(ijb)=T(ipmat(is,is)+nbas(is)*(JB-1)+IB-1)
     &          +T(ipmat(is,is)+nbas(is)*(IB-1)+JB-1)
         End Do
         ijb=ijb+1
         F(ijb)=T(ipmat(is,is)+nbas(is)*(iB-1)+IB-1)
        End Do
       End Do
C      write (*,*) "fockocc in ..."
C      do i = 1, ndens2
C      write (*,'(i3,f20.10)') i,f(i)
C      end do
       Call Put_dArray('FockOcc',F,nDens2)

!       call recprt('RHS',' ',fock,ndens2,1)
*
       Call mma_deallocate(G1r)
       Call mma_deallocate(G2r)
       Call mma_deallocate(T)
       Call mma_deallocate(F)

*
       Return

       Contains

      Subroutine PT2_SLag

      Implicit Real*8 (A-H,O-Z)
      ! integer opout
      Real*8, Allocatable:: CIL(:), CIR(:)
      integer :: i,j

!     At present, Molcas accepts equally-weighted MCSCF reference,
!     so all SLag values are employed in the following computation.
!     For unequally-weighted reference as in GAMESS-US, some more
!     operations are required, but the CP-MCSCF part has to be
!     modified, so this may not be realized easily.

      nConfL=Max(nconf1,nint(xispsm(1,1)))
      nConfR=Max(nconf1,nint(xispsm(1,1)))
      call mma_allocate(CIL, nConfL, Label='CIL')
      call mma_allocate(CIR, nConfR, Label='CIR')
      !! iR = iRLXRoot
      Do jR = 1, nRoots
        Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,1)
        Do kR = 1, jR !! jR-1
          iSLag = jR + nRoots*(kR-1)
          vSLag = +SLag(iSLag)
          write (6,*) "jr,kr,slag = ", jr,kr,vslag
          If (abs(vSLag).le.1.0d-10) Cycle
C
          Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,1)
          ! iRC=opout(ipCI)
          Call CSF2SD(W(ipCI)%Vec(1+(kR-1)*nconf1),CIR,1)
          ! iRC=opout(ipCI)
          ! iRC=ipnout(-1)
          ! icsm=1
          ! issm=1
          Call Densi2(2,G1r,G2r,CIL,CIR,0,0,0,n1dens,n2dens)
          !! For RDM1
          ij=0
          Do i=0,ntAsh-1
            Do j=0,i-1
              ij=ij+1
              G1q(ij)=G1q(ij)+
     *          (G1r(1+i*ntAsh+j)+
     *           G1r(1+j*ntAsh+i))*Half*vSLag
            End Do
            ij=ij+1
            G1q(ij)=G1q(ij) + G1r(1+i*ntAsh+i)*vSLag
          End Do
          !! For RDM2
          Do i=1,ntAsh**2
            j=itri(i,i)
            G2r(j)=Half*G2r(j)
          End Do
          Do i=0,ntAsh-1
            Do j=0,i-1
              ij=i*(i+1)/2+j
              Do k=0,ntAsh-1
                Do l=0,k
                  kl=k*(k+1)/2+l
                  If (ij.ge.kl) Then
                    factor=Quart*vSLag
                    If (ij.eq.kl) factor=Half*vSLag
                    ijkl=ij*(ij+1)/2+kl
                    ij2=i*ntAsh+j
                    kl2=k*ntAsh+l
                    G2q(1+ijkl)=G2q(1+ijkl)
     *                + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                    ij2=Max(j*ntAsh+i,l*ntAsh+k)
                    kl2=Min(j*ntAsh+i,l*ntAsh+k)
                    G2q(1+ijkl)=G2q(1+ijkl)
     &                + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                    If (k.ne.l) Then
                      ij2=i*ntAsh+j
                      kl2=l*ntAsh+k
                      G2q(1+ijkl)=G2q(1+ijkl)
     &                  + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                      If (ij.ne.kl) Then
                        ij2=Max(j*ntAsh+i,k*ntAsh+l)
                        kl2=Min(j*ntAsh+i,k*ntAsh+l)
                        G2q(1+ijkl)=G2q(1+ijkl)
     &                    + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                      End If
                    End If
                  End If
                End Do
              End Do
            End Do
            ij=i*(i+1)/2+i
            Do k=0,ntAsh-1
              Do l=0,k
                kl=k*(k+1)/2+l
                If (ij.ge.kl) Then
                  factor=Half*vSLag
                  If (ij.eq.kl) factor=One*vSLag
                  ijkl=ij*(ij+1)/2+kl
                  ij2=i*ntAsh+i
                  kl2=k*ntAsh+l
                  G2q(1+ijkl)=G2q(1+ijkl)
     *              + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                  If (k.ne.l) Then
                    kl2=l*ntAsh+k
                    G2q(1+ijkl)=G2q(1+ijkl)
     &                + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                  End If
                End If
              End Do
            End Do
          End Do
        End Do
      End Do
      call mma_deallocate(CIL)
      call mma_deallocate(CIR)
      nConf=ncsf(1)
C
      Return
C
      End Subroutine PT2_SLag

      Subroutine SADENS()

      Implicit Real*8 (A-H,O-Z)
      ! integer opout
      Real*8, Allocatable:: CIL(:), CIR(:)
      integer :: i,j

      nConfL=Max(nconf1,nint(xispsm(1,1)))
      call mma_allocate(CIL, nConfL, Label='CIL')
      call mma_allocate(CIR, nConfR, Label='CIR')
      G1q(:) = Zero
      Do jR = 1, nRoots
C       if (jR .ne. iRlxRoot) cycle
        Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,1)
        Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIR,1)
        Call Densi2(1,G1r,G2r,CIL,CIR,0,0,0,n1dens,n2dens)
        write (6,*) "jR = ", jR
        call sqprt(g1r,nash(1))
        !! For RDM1
        ij=0
        Do i=0,ntAsh-1
          Do j=0,i-1
            ij=ij+1
            G1q(ij)=G1q(ij)+
     *        (G1r(1+i*ntAsh+j)+
     *         G1r(1+j*ntAsh+i))*Half
          End Do
          ij=ij+1
          G1q(ij)=G1q(ij) + G1r(1+i*ntAsh+i)
        End Do
      End Do
      call mma_deallocate(CIL)
      call mma_deallocate(CIR)
      nConf=ncsf(1)

      scal = 1.0D+00/DBLE(nRoots)
      call dscal_(ng1,scal,G1q,1)
C
      Return
C
      End Subroutine SADENS
C
       End

      subroutine compute_lag
      use Arrays, only: CMO,CMO_inv,int1,fimo,f0sqmo
      use OneDat, only: sNoNuc, sNoOri
      use ipPage, only: W
      use rctfld_module, only: lRf
      use PCM_grad, only: Preppcm
      implicit real*8 (a-h,o-z)

#include "Input.fh"
#include "stdalloc.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "real.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"

      dimension WRK(nBas(1),nBas(1)),
     * wrk2(nbas(1),nbas(1)),WRK3(nBas(1),nBas(1))
      Real*8, Allocatable:: WRKDER(:,:),G1q(:),G2q(:),G1r(:),G2r(:)
      Real*8, Allocatable:: STmat(:), Smat(:)
      character(len=8) :: Label

      dimension tao(nbas(1),nbas(1),nbas(1),nbas(1)),
     *          tmo(nbas(1)**4)
      integer keep(8)
      logical isq,olag,clag

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      olag = .true.
      clag = .false.

      ng1=itri(ntash,ntash)
      ng2=itri(ng1,ng1)

      Call mma_allocate(G1q,ng1,Label='G1q')
      Call mma_allocate(G2q,ng2,Label='G2q')
      Call mma_allocate(G1r,ntash**2,Label='G1r')
      Call mma_allocate(G2r,itri(ntash**2,ntash**2),Label='G2r')

      call compute_G(1,0,0.0d+00,G1q,G2q,G1r,G2r)
      Call DaName_MF_wa(LuTri1,FnTri1)
      CALL DANAME_wa(LUTRI2,FNTRI2)
      Call DAName_wa(LUHLF2,FNHLF2)
      Call DAName_wa(LUHLF3,FNHLF3)

      Call mma_allocate(STmat,nbas(1)**2,Label='STmat')
      Call mma_allocate(Smat,nbas(1)**2,Label='Smat')
*
      iSymlbl=1
      iOpt=ibset(ibset(0,sNoOri),sNoNuc)
      Label='Mltpl  0'
      iComp=1
      Call RdOne(irc,iOpt,Label,iComp,STmat,iSymlbl)
*
      index = 1
      iOff = 0
      Do iSym = 1, nSym
         Do i = 1, nBas(iSym)
            Do j = 1, i-1
               Smat(j+nBas(iSym)*(i-1)+iOff) =STmat(index)
               Smat(i+nBas(iSym)*(j-1)+iOff) =STmat(index)
               index = index + 1
            End Do
            Smat(i+nBas(iSym)*(i-1)+iOff) =STmat(index)
            index = index + 1
         End Do
         ioff=ioff+nBas(iSym)**2
      End Do
      Call mma_deallocate(STmat)

      idum = 1
      npq = 0
      iopt = 1
      ipq = 0
      lbuf = nbas(1)**2
      Call GetOrd(IRC,isq,nSym,nBas,KEEP)
      tao = 0.0d+00
      nbast = nbas(1)
      nseq = 0
      do ip = 1, nbas(1)
        do iq = 1, ip
          IPQ=IPQ+1
          nseq = nseq + 1
            IF (IPQ.GT.NPQ) THEN
          CALL RDORD(IRC,IOPT,idum,idum,idum,idum,tmo,
     *     ntbtri*ntbsqr,NPQ)
          IF(IRC.GT.1) then
            write (6,*) "wrong"
            call abend
          end if
          IPQ=2
          iopt = 1
            ENDIF
C         CALL SQUARE (WRK,tao(1,1,ip,iq),1,nBas(1),nBas(1))
C         CALL SQUARE (WRK,tao(1,1,iq,ip),1,nBas(1),nBas(1))
      CALL SQUARE (tmo(1+ntbtri*(nseq-1)),tao(1,1,ip,iq),
     *  1,nBas(1),nBas(1))
      CALL SQUARE (tmo(1+ntbtri*(nseq-1)),tao(1,1,iq,ip),
     *  1,nBas(1),nBas(1))
C         write (6,*) "ip,iq = ", ip,iq
C         call sqprt(tao(1,1,ip,iq),nbas(1))
        end do
      end do
C     call abend

      tmo = 0.0d+00
      call transformation(CMO,tao,tmo)

      !! initialize orbital Lagrangian
      WRK = 0.0d+00

      !call daxpy_(nish(1)*norb(1),1.0d+00,int1,1,wrk,1)
      !call daxpy_(nish(1)*norb(1),1.25d+00,fimo,1,wrk,1)
C     call dgemm_('T','N',norb(1),nash(1),nash(1),
C    *            1.0d+00,fimo(1+nish(1)),norb(1),G1r,nash(1),
C    *            0.0d+00,WRK(1,nish(1)+1),norb(1))
C     call dgemm_('N','T',norb(1),nash(1),nash(1),
C    *            1.0d+00,fimo(1+norb(1)*nish(1)),norb(1),
C    *                    G1r,nash(1),
C    *            1.0d+00,WRK(1,nish(1)+1),norb(1))

      e2 = 0.0d+00
      Do iB=1,ntash
       Do jB=1,ntash
        iDij=iTri(ib,jB)
        iRij=jb+(ib-1)*ntash
        Call Coul(1,1,1,1,nish(1)+iB,nish(1)+jB,wrk2,wrk3)
        Do kB=1,ntash
         Do lB=1,ntash
          iDkl=iTri(kB,lB)
          iRkl=lb+(kb-1)*ntash
          fact=One
          if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
          if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
          iijkl=itri(iDij,iDkl)
          iRijkl=itri(iRij,iRkl)
       !  G2r(iRijkl)=Fact*G2q(iijkl)
      ! e2 = e2 + t(nish(1)+kb+nbas(1)*(nish(1)+lb-1))*g2r(irijkl)
      ! Call Coul(1,1,1,1,nish(1)+iB,nish(1)+jB,wrk2,wrkder)
      !    do iq = 1, nbas(1)
      !    wrk(iq,nish(1)+kb) = wrk(iq,nish(1)+kb)
     *!      + wrk2(iq,nish(1)+lb)*g2r(irijkl)*0.25d+00
      !    wrk(iq,nish(1)+lb) = wrk(iq,nish(1)+lb)
     *!      + wrk2(nish(1)+kb,iq)*g2r(irijkl)*0.25d+00
      !    end do
      ! Call Coul(1,1,1,1,nish(1)+kB,nish(1)+lB,wrk2,wrkder)
      !    do iq = 1, nbas(1)
      !    wrk(iq,nish(1)+ib) = wrk(iq,nish(1)+ib)
     *!      + wrk2(iq,nish(1)+jb)*g2r(irijkl)*0.25d+00
      !    wrk(iq,nish(1)+jb) = wrk(iq,nish(1)+jb)
     *!      + wrk2(nish(1)+ib,iq)*g2r(irijkl)*0.25d+00
      !    end do
         End Do
        End Do
       End Do
      End Do

      delta0 = 1.0d-06
      if (olag) then
        Call mma_allocate(WRKDER,nbas(1),nbas(1),Label='WRKDER')
        WRKDER = 0.0d+00

        do mu = 1, nbas(1)
        do nu = 1, nbas(1)
        do ivib = 1, 2
          if (ivib.eq.1) delta =  delta0
          if (ivib.eq.2) delta = -delta0

          CMO(mu+nbas(1)*(nu-1)) = CMO(mu+nbas(1)*(nu-1)) + delta

          !! transform two-electron integrals
          Call mma_allocate(CMO_inv,nbas(1)**2,Label='CMO_Inv')
          Call dGemm_('T','N', nOrb(1),nBas(1),nBas(1),
     &               1.0d0,CMO(1),nBas(1),Smat(1),nBas(1),
     &               0.0d0,CMO_Inv(1),nOrb(1))
          Call SetUp_CASPT2_Tra(nSym,nBas,nOrb,nIsh,nAsh,
     &                          nFro,nDel,CMO,nDens2,
     &                          LuTri1,LuTri2,LuHlf2,LuHlf3)
          call mma_deallocate(cmo_inv)
          iType=3  ! Means that TraCtl is called by MCLR
          Call TraCtl_Drv(iType,.True.,1)
          Call DaClos(LuTri2)
          Call DaClos(LuHlf2)
          Call DaClos(LuHlf3)
C       call trctl_mclr
          tmo = 0.0d+00
          call transformation(CMO,tao,tmo)

          if (lRF) call preppcm()

          !! compute energies
          call compute_energy(eee,G1r,G2r,tmo)
C         Call DaClos(LuTri2)
C         Call DaClos(LuTri3)
C         Call DaClos(LuTri4)
C         Call DaClos(LuTri5)

          CMO(mu+nbas(1)*(nu-1)) = CMO(mu+nbas(1)*(nu-1)) - delta

          if (ivib.eq.1) e1 = eee
          if (ivib.eq.2) e2 = eee

          if (ivib.eq.2) then
            eder = (e1-e2)/(2.0d+00*delta)
            WRKDER(mu,nu) = eder
          end if
        end do
        end do
        end do

        Call DGEMM_('T','N',nBas(1),nBas(1),nBas(1),
     *              1.0D+00,CMO,nBas(1),WRKDER,nBas(1),
     *              1.0D+00,WRK,nBas(1))
        write (6,*) "orbital lagrangian"
        call sqprt(WRK,nBas(1))

        wrkder = 0.0d+00
C       Call DGeSub(WRK,nBas(1),'N',WRK,nBas(1),'T',
C    *              WRKDER,nBas(1),nBas(1),nBas(1))
        do i = 1, nbas(1)
          do j = 1, nbas(1)
            tmp = wrk(i,j) - wrk(j,i)
            wrkder(i,j) = tmp
            wrkder(j,i) =-tmp
          end do
        end do
        write (6,*) "Anti-Symmetrized Orbital Lagrangian"
        Call SqPrt(WRKDER,nBas(1))

        write (6,*) "after addgrad2"
        call addgrad2(wrkder,1,1.0d+00)
        Call SqPrt(WRKDER,nBas(1))

        Call DGeSub(F0SQMO(1),nOrb(1),'N',
     &              F0SQMO(1),nOrb(1),'T',
     &              WRKDER,nOrb(1),
     &              nOrb(1),nOrb(1))
        Call SqPrt(WRKDER,nBas(1))

        call mma_deallocate(WRKDER)
      else if (clag) then
        nConfL=Max(nconf1,nint(xispsm(1,1)))
        Call mma_allocate(WRKDER,nconfl,nroots,Label='WRKDER')
        WRKDER = 0.0d+00

        do nu = 1, nroots
        do mu = 1, ncsf(1) ! nconfl
        do ivib = 1, 2
          if (ivib.eq.1) delta =  delta0
          if (ivib.eq.2) delta = -delta0

          W(ipCI)%Vec(mu+(nu-1)*ncsf(1))
     *      = W(ipCI)%Vec(mu+(nu-1)*ncsf(1)) + delta
C       Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,1)

          if (lRF) call PrepPCM()

          call compute_G(1,mu,delta,G1q,G2q,G1r,G2r)

          !! compute energies
          call compute_energy(eee,G1r,G2r,tmo)

          W(ipCI)%Vec(mu+(nu-1)*ncsf(1))
     *      = W(ipCI)%Vec(mu+(nu-1)*ncsf(1)) - delta

          if (ivib.eq.1) e1 = eee
          if (ivib.eq.2) e2 = eee

          if (ivib.eq.2) then
            eder = (e1-e2)/(2.0d+00*delta)
            WRKDER(mu,nu) = eder
          end if
        end do
        end do
        end do

        write (6,*) "clag after the loop"
        do i = 1, nroots
        write (6,'("roots = ",i3)') i
        do j = 1, ncsf(1) ! nconfl
          write (6,'(i4,f20.10)') j,wrkder(j,i)
        end do
        end do

        write (6,*) "Numerical State Lagrangian and projection"
        call proj(wrkder)
C       ijRoots = 0
C       Do iRoots = 1, nRoots
C         Do jRoots = 1, iRoots-1
C           ijRoots = ijRoots + 1
C           SLagV = DDot_(nConf,Work(ipWRK1+nConf*(iRoots-1)),1,
C    *                          Work(ipWRK2+nConf*(jRoots-1)),1)
C    *            - DDot_(nConf,Work(ipWRK1+nConf*(jRoots-1)),1,
C    *                          Work(ipWRK2+nConf*(iRoots-1)),1)
C           write (6,'(2i3,f20.10)') iRoots,jRoots,SLagV
C         End Do
C       End Do
C       End If

        ! projection
        write (6,*) "clag after projection"
        do i = 1, nroots
        write (6,'("roots = ",i3)') i
        do j = 1, ncsf(1) ! nconfl
          write (6,'(i4,f20.10)') j,wrkder(j,i)
        end do
        end do
C       write (6,*) "After Projection"
C       Do iRoots = 1, nRoots
C         write (6,*) "iRoots = ", iRoots
C         Do jRoots = 1, nRoots
C           Scal = DDot_(nConf,Work(ipWRK1+nConf*(jRoots-1)),1,
C    *                         Work(ipWRK2+nConf*(iRoots-1)),1)
C           Call DaXpY_(nConf,-Scal,Work(ipWRK1+nConf*(jRoots-1)),1,
C    *                              Work(ipWRK2+nConf*(iRoots-1)),1)
C         End Do
C         Do indCI = 1, nConf
C           write (6,'(i5,f20.10)') indCI,Work(ipWRK2+indCI-1)
C         End Do
C       End Do

        Call mma_deallocate(WRKDER)
      end if

      Call mma_deallocate(G1q)
      Call mma_deallocate(G2q)
      Call mma_deallocate(G1r)
      Call mma_deallocate(G2r)
      Call mma_deallocate(smat)

      end subroutine compute_lag

      subroutine proj(derci)

      use ipPage, only: W

      implicit real*8 (a-h,o-z)

#include "Input.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "Pointers.fh"
#include "real.fh"
#include "sa.fh"
#include "stdalloc.fh"

      dimension :: derci(max(nconf1,nint(xispsm(1,1))),nRoots)

      Real*8, Allocatable:: CI(:,:),DERtmp(:,:)

      nConfL=Max(nconf1,nint(xispsm(1,1)))
      call mma_allocate(CI, ncsf(1), nRoots, Label='CI')
      call mma_allocate(DERtmp, ncsf(1), nRoots, Label='DERtmp')
C     write (6,*) "nconf etc. = ", nconf,nconf1,nconfl,ncsf(1)

      Do iRoots = 1, nRoots
C       Call CSF2SD(W(ipCI)%Vec(1+(iRoots-1)*ncsf(1)),CI(1,iRoots),1)
        call dcopy_(ncsf(1),W(ipCI)%Vec(1+(iRoots-1)*ncsf(1)),1,
     *  CI(1,iRoots),1)
      End Do

      !! Numerical state lagrangian (should be zero)
      write (6,'(" Numerical State Lagrangian (should be zero)")')
      ijRoots = 0
      Do iRoots = 1, nRoots
        Do jRoots = 1, iRoots-1
          ijRoots = ijRoots + 1
          SLagV = DDot_(ncsf(1),CI(1,iRoots),1,DERCI(1,jRoots),1)
     *          - DDot_(ncsf(1),CI(1,jRoots),1,DERCI(1,iRoots),1)
          write (6,'(2i3,f20.10)') iRoots,jRoots,SLagV
        End Do
      End Do

      !! Projection
      DERtmp = DERCI
      Do iRoots = 1, nRoots
        Do jRoots = 1, nRoots
          Scal = DDot_(ncsf(1),DERtmp(1,iRoots),1,CI(1,jRoots),1)
          Call DaXpY_(ncsf(1),-Scal,CI(1,jRoots),1,DERCI(1,iRoots),1)
        End Do
      End Do

      call mma_deallocate(CI)
      call mma_deallocate(DERtmp)
      nConf=ncsf(1)

      end subroutine proj

      subroutine compute_G(imode,idelta,delta,G1q,G2q,G1r,G2r)

      use ipPage, only: W

      implicit real*8 (a-h,o-z)

#include "Input.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "Pointers.fh"
#include "real.fh"
#include "sa.fh"
#include "stdalloc.fh"

      dimension :: G1q(*),G2q(*),G1r(*),G2r(*)

      Real*8, Allocatable:: CIL(:), CIR(:)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      nConfL=Max(nconf1,nint(xispsm(1,1)))
      nConfR=Max(nconf1,nint(xispsm(1,1)))
      call mma_allocate(CIL, nConfL, Label='CIL')
      call mma_allocate(CIR, nConfR, Label='CIR')

      ng1 = ntash*(ntash-1)/2 + ntash
      ng2 = ng1*(ng1-1)/2 + ng1
      ng3 = ntash**2
      ng4 = ng3*(ng3-1)/2 + ng3

      G1q(1:ng1) = 0.0d+00
      G2q(1:ng2) = 0.0d+00
      G1r(1:ng3) = 0.0d+00
      G2r(1:ng4) = 0.0d+00

      !! iR = iRLXRoot
      Do jR0 = 1, nRoots
        if (imode.eq.1) then
          !! state-specific
          if (jR0.ne.iRlxRoot) cycle
          vSLag = 1.0d+00
          if (isNAC)  then
            jR = NSSA(2)
            kR = NSSA(1)
          else
            jR = jR0
            kR = jR0
          end if
        else if (imode.eq.0) then
          vSLag = 1.0d+00/dble(nRoots)
          jR = jR0
        else
          write (6,*) "unknown imode in compute_G"
          call abend
        end if
        Call CSF2SD(W(ipCI)%Vec(1+(jR-1)*nconf1),CIL,1)
        ! iRC=opout(ipCI)
        Call CSF2SD(W(ipCI)%Vec(1+(kR-1)*nconf1),CIR,1)


C       if (idelta.eq.0) then
C         !! unperturbed density
C       else
C         !! perturbed density
C         istate = idelta/nconfl + 1
C         ivec   = mod(idelta,nconfl)
C         if (ivec.eq.0) then
C           istate = istate-1
C           ivec = nconfl
C         end if
C         if (istate.eq.jR0) then
C           CIL(ivec) = CIL(ivec) + delta
C           CIR(ivec) = CIR(ivec) + delta
C         end if
C       end if

        ! iRC=opout(ipCI)
        ! iRC=ipnout(-1)
        ! icsm=1
        ! issm=1
        Call Densi2(2,G1r,G2r,CIL,CIR,0,0,0,n1dens,n2dens)
        !! For RDM1
        ij=0
        Do i=0,ntAsh-1
          Do j=0,i-1
            ij=ij+1
            G1q(ij)=G1q(ij)+
     *        (G1r(1+i*ntAsh+j)+
     *         G1r(1+j*ntAsh+i))*Half*vSLag
          End Do
          ij=ij+1
          G1q(ij)=G1q(ij) + G1r(1+i*ntAsh+i)*vSLag
        End Do
        !! For RDM2
        Do i=1,ntAsh**2
          j=itri(i,i)
          G2r(j)=Half*G2r(j)
        End Do
        Do i=0,ntAsh-1
          Do j=0,i-1
            ij=i*(i+1)/2+j
            Do k=0,ntAsh-1
              Do l=0,k
                kl=k*(k+1)/2+l
                If (ij.ge.kl) Then
                  factor=Quart*vSLag
                  If (ij.eq.kl) factor=Half*vSLag
                  ijkl=ij*(ij+1)/2+kl
                  ij2=i*ntAsh+j
                  kl2=k*ntAsh+l
                  G2q(1+ijkl)=G2q(1+ijkl)
     *              + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                  ij2=Max(j*ntAsh+i,l*ntAsh+k)
                  kl2=Min(j*ntAsh+i,l*ntAsh+k)
                  G2q(1+ijkl)=G2q(1+ijkl)
     &              + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                  If (k.ne.l) Then
                    ij2=i*ntAsh+j
                    kl2=l*ntAsh+k
                    G2q(1+ijkl)=G2q(1+ijkl)
     &                + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                    If (ij.ne.kl) Then
                      ij2=Max(j*ntAsh+i,k*ntAsh+l)
                      kl2=Min(j*ntAsh+i,k*ntAsh+l)
                      G2q(1+ijkl)=G2q(1+ijkl)
     &                  + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                    End If
                  End If
                End If
              End Do
            End Do
          End Do
          ij=i*(i+1)/2+i
          Do k=0,ntAsh-1
            Do l=0,k
              kl=k*(k+1)/2+l
              If (ij.ge.kl) Then
                factor=Half*vSLag
                If (ij.eq.kl) factor=One*vSLag
                ijkl=ij*(ij+1)/2+kl
                ij2=i*ntAsh+i
                kl2=k*ntAsh+l
                G2q(1+ijkl)=G2q(1+ijkl)
     *            + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                If (k.ne.l) Then
                  kl2=l*ntAsh+k
                  G2q(1+ijkl)=G2q(1+ijkl)
     &              + factor*G2r(1+ij2*(ij2+1)/2+kl2)
                End If
              End If
            End Do
          End Do
        End Do
      End Do
      call mma_deallocate(CIL)
      call mma_deallocate(CIR)
      nConf=ncsf(1)

       Do iB=1,ntash
        Do jB=1,ntash
        G1r(ib+(jb-1)*ntash) = G1q(itri(ib,jb))
        End Do
       End Do

       Do iB=1,ntash
        Do jB=1,ntash
         iDij=iTri(ib,jB)
         iRij=jb+(ib-1)*ntash
         Do kB=1,ntash
          Do lB=1,ntash
           iDkl=iTri(kB,lB)
           iRkl=lb+(kb-1)*ntash
           fact=One
           if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
           if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
           iijkl=itri(iDij,iDkl)
           iRijkl=itri(iRij,iRkl)
           G2r(iRijkl)=Fact*G2q(iijkl)
          End Do
         End Do
        End Do
       End Do

      end subroutine compute_G

      subroutine transformation(CMO,TAO,TMO)

      implicit real*8 (a-h,O-z)

#include "Input.fh"

      dimension wrk1(nbas(1),nbas(1),nbas(1),nbas(1)),
     *          wrk2(nbas(1),nbas(1),nbas(1),nbas(1)),
     *          tao(nbas(1),nbas(1),nbas(1),nbas(1)),
     *          tmo(nbas(1),nbas(1),nbas(1),nbas(1))
      dimension cmo(*)

      ! 4th
      nbast=nbas(1)
      call dgemm_('N','N',nOrb(1)**3,nBasT,nOrb(1),
     *            1.0D+00,TAO,nOrb(1)**3,CMO,nBasT,
     *            0.0D+00,WRK1,nOrb(1)**3)
      ! 3rd
      do i = 1, nBasT
        call dgemm_('N','N',nOrb(1)**2,nBasT,nOrb(1),
     *            1.0D+00,WRK1(1,1,1,i),nOrb(1)**2,CMO,nBasT,
     *            0.0D+00,WRK2(1,1,1,i),nOrb(1)**2)
      end do

      do j = 1, nBasT
      do i = 1, nBasT
        ! 1st
        call dgemm_('T','N',nBasT,nOrb(1),nOrb(1),
     *            1.0D+00,CMO,nBasT,WRK2(1,1,i,j),nOrb(1),
     *            0.0D+00,WRK1,nBasT)
        ! 2nd
        call dgemm_('N','N',nBasT,nBasT,nOrb(1),
     *              1.0D+00,WRK1,nBasT,CMO,nBasT,
     *              0.0D+00,TMO(1,1,i,j),nBasT)
      end do
      end do

      end subroutine transformation

      subroutine compute_energy(eee,G1r,G2r,tmo)

      use Arrays, only: Int1, fimo, FAMO, F0SQMO, INT2,CMO
      use OneDat, only: sOpSiz
      use rctfld_module, only: lRf
      use PCM_grad, only: PCM_grad_dens,PCM_grad_dens2,DSSMO,DSSAO,
     *                    PCMSSMO,dscfmo,dscfao,pcmscfmo,potnuc_pcm,
     *                    pcmscfao,pcmssmo,pcmssao,def_solv

      implicit real*8 (a-h,o-z)

#include "Input.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "Pointers.fh"
#include "sa.fh"

      Character*8 Label

      real*8, intent(in) :: G1r(*),G2r(*)
      Real*8, Allocatable:: T(:), F(:)
      Real*8, Allocatable:: Q(:), Tmp2(:,:), T3(:)
      dimension temp1(nbas(1)**2),tmo(nbas(1),nbas(1),nbas(1),nbas(1))

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      imode = 0

      Call mma_allocate(T,nbas(1)**2,Label='T')
      Call mma_allocate(F,nbas(1)**2,Label='F')

      ! transform H and FI

      iRc=-1
      iOpt=ibset(0,sOpSiz)
      iisym=2**0
      Label='ONEHAM'
      iComp=1
      iisym=2**0
      iRc=-1
      iOpt=0
      Call RdOne(iRc,iOpt,Label,iComp,Temp1,iisym)
      Call Square(Temp1,t,1,nBas(1),nBas(1))
      if (lRF) then
        call daxpy_(nbas(1)**2,1.0d+00,pcmscfao(1,1),1,t,1)
      end if
      Call DGEMM_('T','N',nOrb(1),nBas(1),nBas(1),
     &            1.0d0,CMO(1),nBas(1),t,nBas(1),
     &            0.0d0,f,nOrb(1))
      Call DGEMM_('N','N',nOrb(1),nOrb(1),nBas(1),
     &            1.0d0,f,nOrb(1),CMO(1),nBas(1),
     &            0.0d0,Int1(1),nOrb(1))
      ! call dcopy_(nbas(1)**2,int1,1,t,1)
      ! call dcopy_(nbas(1)**2,[0.0d+00],0,int1,1)
C     call FckMat
      Call mma_allocate(Q,nDens2,Label='Q')
      Call mma_allocate(Tmp2,ndens2,2,Label='Tmp2')
      Call mma_allocate(T3,ndens2,Label='T3')
      Call Read22_2(Int2,F0SQMO,Q,FIMO,FAMO,Tmp2(:,1),Tmp2(:,2),T3)
      Call mma_deallocate(T3)
      Call mma_deallocate(Tmp2)
      Call mma_deallocate(Q)
      ! call dcopy_(nbas(1)**2,t,1,int1,1)
      ! call daxpy_(nbas(1)**2,1.0d+00,int1,1,fimo,1)

      enuc = potnuc - potnuc_pcm !! 28.65127677d+00
      if (isNAC) enuc = -potnuc_pcm

C     write (6,*) "int1"
C     call sqprt(int1,nbas(1))
C     write (6,*) "fimo"
C     call sqprt(fimo,nbas(1))
      e1 = 0.0d+00
      e1fimo = 0.0d+00
      e12 = 0.0d+00
      e12vac = 0.0d+00
      !! At present, the PCM contribution to the energy is
      !! D^SS * (V(N,SA) + V(e,SA))
      !! but the correct expression is
      !! D^SS * (V(N,SA) + V(e,SA)) / 2 + D^SA * (V(N,SS) + V(e,SS)) / 2
      !! so we need some corrections
      do i = 1, nish(1)
        nseq = i*(i-1)/2 + i
        !! D^SS * (V(N,SA) + V(e,SA))
        if (.not.isNAC) then
          e1 = e1 + 2.0d+00*int1(i+nbas(1)*(i-1))
          e1fimo = e1fimo + 2.0d+00*fimo(i+nbas(1)*(i-1))
          e12vac = e12vac
     *      + 2.0d+00*(int1(i+nbas(1)*(i-1))+fimo(i+nbas(1)*(i-1)))
        end if
        if (lRF) then
          !! - D^SS * V(e) / 2
C         e1fimo = e1fimo - 2.0d+00*pcmscfmo(i+nbas(1)*(i-1),3)

          if (def_solv.eq.1) then
            !! - D^SS * V(e,SS) / 2
            e1fimo = e1fimo - 2.0d+00*pcmscfmo(i+nbas(1)*(i-1),3)
          else if (def_solv.eq.2) then
            !! - D^SS * (V(N,SA)+V(e,SA)) / 2
            e1fimo = e1fimo - 2.0d+00*pcmscfmo(i+nbas(1)*(i-1),1)
            !! + D^SA * V(N,SS) / 2
            e1fimo = e1fimo + 2.0d+00*pcmssmo(i+nbas(1)*(i-1),2)
          else if (def_solv.eq.3) then
            !! - D^SS * V(e,SA) / 2
            if (.not.isNAC) then
              e1fimo = e1fimo - 2.0d+00*pcmscfmo(i+nbas(1)*(i-1),3)
            end if
          else if (def_solv.eq.4) then
            !! - D^SS * V(e,SA) / 2
            e1fimo = e1fimo - 2.0d+00*pcmscfmo(i+nbas(1)*(i-1),3)
          else if (def_solv.eq.5) then
            !! - D^SS * V(e,SA) / 2
            e1fimo = e1fimo - 2.0d+00*pcmscfmo(i+nbas(1)*(i-1),3)
          else if (def_solv.eq.6) then
            !! - D^SS * V(e,SA) / 2
            e1fimo = e1fimo - 2.0d+00*pcmscfmo(i+nbas(1)*(i-1),3)
          end if
        end if
      end do
      e12 = (e1fimo+e1)*0.5d+00
      e12vac = e12vac * 0.5d+00
C     write (6,'("one-electron energy1= ",f20.10)') e1
C     write (6,'("FIMO         energy1= ",f20.10)') e1fimo
C     write (6,'("e12 = rcore/2       = ",f20.10)') e12
C     aaa = 0.0d+00
C     bbb = 0.0d+00
C     do i = 1, nish(1)
C       aaa = aaa + pcmscfmo(i+nbas(1)*(i-1),3)
C       bbb = bbb + pcmssmo(i+nbas(1)*(i-1),3)
C     end do
C     write (6,*) aaa,bbb
C     aaa = ddot_(nbas(1)**2,pcmscfao,1,dssao,1)
C     bbb = ddot_(nbas(1)**2,pcmssao,1,dscfao,1)
C     do i = 1, 10
C     write (6,'(i4,2f20.10)') i,pcmscfao(i,1),pcmssao(i,1)
C     end do
C     write (6,*) aaa,bbb
C     aaa = ddot_(nbas(1)**2,pcmscfao(1,2),1,dssao,1)
C     bbb = ddot_(nbas(1)**2,pcmssao(1,2),1,dscfao,1)
C     do i = 1, 10
C     write (6,'(i4,2f20.10)') i,pcmscfao(i,2),pcmssao(i,2)
C     end do
C     write (6,*) aaa,bbb
C     aen = aaa
C     ben = bbb
C     aaa = ddot_(nbas(1)**2,pcmscfao(1,3),1,dssao,1)
C     bbb = ddot_(nbas(1)**2,pcmssao(1,3),1,dscfao,1)
C     do i = 1, 10
C     write (6,'(i4,2f20.10)') i,pcmscfao(i,3),pcmssao(i,3)
C     end do
C     write (6,*) aaa,bbb
C     aee = aaa
C     bee = bbb

C     write (6,*) "asdf"
C     write (6,*) aen+aee, ben+bee

C     call abend

      e1fimoact = e1fimo
      e1sav = e1
      erfx = 0.0d+00
      ecorr = 0.0d+00
C     write (6,*) "g1r in compute"
C     call sqprt(g1r,ntash)
C     write (6,*) "FIMO"
C     call sqprt(fimo,nbas(1))
      do i = 1, nash(1)
        do j = 1, nash(1)
          nseq = max(nish(1)+i,nish(1)+j)
     *    *(max(nish(1)+i,nish(1)+j)-1)/2+min(nish(1)+i,nish(1)+j)
          e1 = e1 + int1(nish(1)+i+nbas(1)*(nish(1)+j-1))
     *            *g1r(i+ntash*(j-1))
          !! D^SS * (V(N) + V(e))
          e1fimo = e1fimo
     *    + g1r(i+ntash*(j-1))*fimo(nish(1)+i+nbas(1)*(nish(1)+j-1))
          e12vac =  e12vac +
     *    + g1r(i+ntash*(j-1))*fimo(nish(1)+i+nbas(1)*(nish(1)+j-1))

          if (lRF) then
          if (def_solv.eq.1) then
            !! - D^SS * V(e,SS) / 2
C           e1fimo = e1fimo
C    *        - 1.0d+00*g1r(i+ntash*(j-1))
C    *                 *pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
            erfx = erfx
     *      + dscfmo(i,j)*pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
          else if (def_solv.eq.2) then
            !! - D^SS * (V(N,SA)+V(e,SA)) / 2
            e1fimo = e1fimo
     *        - 0.5d+00*g1r(i+ntash*(j-1))
     *                 *pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),1)
            !! + D^SA * V(N,SS) / 2
     *        + 0.5d+00*dscfmo(i,j)
     *                 *pcmssmo(nish(1)+i+nbas(1)*(nish(1)+j-1),2)
          else if (def_solv.eq.3) then
            if (.not.isNAC) then
              erfx = erfx
     *        + dscfmo(i,j)*pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
            end if
          else if (def_solv.eq.4) then
            erfx = erfx
     *      + dscfmo(i,j)*pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
            ecorr =  ecorr
     *        + 0.5d+00*(dscfmo(i,j)-dssmo(i,j))
     *        *pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
          else if (def_solv.eq.5) then
            erfx = erfx
     *      + dscfmo(i,j)*pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
            ecorr =  ecorr
     *        + (dscfmo(i,j)-dssmo(i,j))
     *        *pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),1)
          else if (def_solv.eq.6) then
            erfx = erfx
     *      + dscfmo(i,j)*pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
            ecorr =  ecorr
     *        + (dscfmo(i,j)-dssmo(i,j))
     *        *pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),2)
          end if


C       erfx = erfx
C    *  + dscfmo(i,j)*pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
C       ecorr =  ecorr
C    *  + 0.5d+00*(dscfmo(i,j)-dssmo(i,j))
C    *    *pcmscfmo(nish(1)+i+nbas(1)*(nish(1)+j-1),3)
          end if
        end do
      end do
      e1fimoact = e1fimo-e1fimoact
C     write (6,*) "e1fimoact = ", e1fimoact

      write (6,*) "e12 = ", e12
      write (6,*) "e1fimoact = ", e1fimoact
      write (6,*) "enuc = ", enuc
      write (6,*) "potnuc_pcm = ", potnuc_pcm
      write (6,*) "erfx etc. = ", -0.5d+00*erfx
      write (6,*) "ecorr = ", ecorr

      e1cont = e12 + enuc + potnuc_pcm - 0.5d+00*erfx
     *       + e1fimoact + ecorr
      write (6,*) "e1cont = ", e1cont

      e2 = 0.0d+00
      Do iB=1,ntash
       Do jB=1,ntash
        iDij=iTri(ib,jB)
        iRij=jb+(ib-1)*ntash
        Call Coul(1,1,1,1,nish(1)+iB,nish(1)+jB,T,F)
        Do kB=1,ntash
         Do lB=1,ntash
          iDkl=iTri(kB,lB)
          iRkl=lb+(kb-1)*ntash
          fact=One
          if(iDij.ge.iDkl .and. kB.eq.lB) fact=Two
          if(iDij.lt.iDkl .and. iB.eq.jB) fact=Two
          iijkl=itri(iDij,iDkl)
          iRijkl=itri(iRij,iRkl)
C         G2r(iRijkl)=Fact*G2q(iijkl)
        e2 = e2 + t(nish(1)+kb+nbas(1)*(nish(1)+lb-1))*g2r(irijkl)
C       e2 = e2 + tmo(nish(1)+ib,nish(1)+jb,nish(1)+kb,nish(1)+lb)
C    *   *g2r(irijkl)
C     write (6,*)  t(nish(1)+kb+nbas(1)*(nish(1)+lb-1)),
C    *   tmo(nish(1)+ib,nish(1)+jb,nish(1)+kb,nish(1)+lb)
         End Do
        End Do
       End Do
      End Do
C     call abend

      e2 = e2 * 0.5d+00
      write (6,*) "two-ele = ", e2

      eee = e1cont + e2

      write (6,*) "active energy = ", e1fimoact + e2 -0.5d+00*erfx
     *  + ecorr
      write (6,*) "eigen energy = ", e12 + enuc + potnuc_pcm + e1fimoact
     *  + e2

      write (6,*) "energy = ", eee
      ! wrong
!     write (6,*) "eigen energy = ", e12vac + enuc + potnuc_pcm
!    * + e1fimoact + e2
!     write (6,*) "eigen energy = ", e12vac + e1fimoact + e2
C     call abend()

      Call mma_deallocate(T)
      Call mma_deallocate(F)

      end subroutine compute_energy
