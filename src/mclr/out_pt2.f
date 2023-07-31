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
      SubRoutine Out_Pt2(iKapDisp,iCIDisp)
********************************************************************
*                                                                  *
********************************************************************
      use Arrays, only: CMO
      use ipPage, only: W
      Implicit Real*8 (a-h,o-z)
#include "detdim.fh"
#include "Input.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "disp_mclr.fh"
#include "cicisp_mclr.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"
#include "SysDef.fh"
      Character*8 Method
      Logical CI, Is_Roots_Set
      Character(LEN=80) Note
! Added for DMRG calculation
      real*8,allocatable::tmpDe(:,:),tmpP(:),tmpDeM(:,:),tmpPM(:,:,:,:)
      Integer iKapDisp(nDisp),iCiDisp(nDisp)
      Character(Len=16) mstate
      Dimension rdum(1),idum(7,8)
      Real*8, Allocatable:: D_K(:), Tmp(:), K1(:), K2(:), DAO(:),
     &                      D_CI(:), D1(:), P_CI(:), P1(:), Conn(:),
     &                      OCCU(:), CMON(:), DTmp(:), G1q(:), G1m(:),
     &                      Temp(:), tTmp(:), DM(:), DMs(:)
*                                                                      *
************************************************************************
*                                                                      *
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
       isym=1
       CI=.true.
       Call Setup_MCLR(iSym)
       nbas_tot=0
       ntot1=0
       nDLMO=0
       nLCMO=0
       Do is=1,nsym
          nbas_tot=nbas_tot+nbas(is)
          ntot1=ntot1+nbas(is)*(nbas(is)+1)/2
          nDLMO=nDLMO+nash(is)
          nLCMO=nLCMO+nbas(is)*nbas(is)
       End Do
       nNAC=(nDLMO+nDLMO**2)/2
       nDLMO=nDLMO*(nDLMO+1)/2
       nPLMO=nDLMO*(nDLMO+1)/2
*
       Call mma_allocate(K1,  nDens2,Label='K1')
       Call mma_allocate(K2,  nDens2,Label='K2')
       Call mma_allocate(DAO,nDens2,Label='DAO')
       Call mma_allocate(D_CI,n1Dens,Label='D_CI')
       Call mma_allocate(D1,n1Dens,Label='D1')
       Call mma_allocate(P_CI,n2Dens,Label='P_CI')
       Call mma_allocate(P1,n2Dens,Label='P1')
       Call mma_allocate(Conn,nDens2,Label='Conn')
       Call mma_allocate(OCCU,nbas_tot,Label='OCCU')
       Call mma_allocate(CMON,ndens2,Label='CMON')
*      OBS nBuf might not be def.
       Call mma_MaxDBLE(nBuf)
       Call mma_allocate(Dtmp,nDens2,Label='DTmp')

*
*
*
* 1)   CI Part
*
*      All multipliers are introduced as densities
*

       If (CI) Then
         nconf1=ncsf(State_sym)
         ilen=nconf1*nroots ! nroot = # of roots in SA
         ipcip=ipget(nconf1*nroots)
         iDisk=iCIDisp(1)
         irc=ipin(ipCIp)
         Call dDaFile(LuTemp,2,W(ipCIp)%Vec,iLen,iDisk)
*
*-------Calculate the densities that correct the nonvariational CI stuff
*
         Call CIDens_sa(.true.,ipCIp,ipCI,
     &                  State_sym,State_sym,
     &                  P_CI,D_CI) ! \bar{d} and \bar{D}

! ======================================================================
         if(doDMRG)then  ! yma
           call dmrg_dim_change_mclr(LRras2(1:8),ntash,0)
           call dmrg_dim_change_mclr(RGras2(1:8),ndim,0)

           call mma_allocate(tmpDe,ndim,ndim,Label='TmpDe')
           call mma_allocate(tmpP,ndim**2*(ndim**2+1)/2,Label='tmpP')
           call mma_allocate(tmpDeM,ntash,ntash,Label='tmpDeM')
           call mma_allocate(tmpPM,ntash,ntash,ntash,ntash,
     &                       Label='tmpPM')
           tmpDe=0.0d0
           tmpP=0.0d0
           tmpDeM=0.0d0
           tmpPM=0.0d0

           ij=0
           do i=1,ntash
             do j=1,ntash
               ij=ij+1
               if(abs(D_CI(ij)).lt.1.0e-12)then
                 D_CI(ij)=0.0d0
               end if
               tmpDeM(i,j)=D_CI(ij)
             end do
           end do

           ij=0
           do i=1,ndim
             do j=1,ndim
               ij=ij+1
               if(i.gt.ntash.or.j.gt.ntash)then
                 tmpDe(i,j)=0.0d0
               else
                 tmpDe(i,j)=tmpDeM(i,j)
               end if
             end do
           end do

           Do i=1,ntash
             Do j=1,ntash
               Do k=1,ntash
                 Do l=1,ntash
                   ij1=ntash*(i-1)+j
                   ij2=ntash*(j-1)+i
                   kl1=ntash*(k-1)+l
                   kl2=ntash*(l-1)+k
                   if(ij1.ge.kl1)then
                   if(abs(P_CI(itri(ij1,kl1))).lt.1.0e-12)then
                     P_CI(itri(ij1,kl1))=0.0d0
                   end if
                   tmpPM(i,j,k,l)=P_CI(itri(ij1,kl1))
                   end if
                 End Do
               End Do
             End Do
           End Do

           do i=1,ndim
             do j=1,ndim
               do k=1,ndim
                 do l=1,ndim
                   ij1=ndim*(i-1)+j
                   ij2=ndim*(j-1)+i
                   kl1=ndim*(k-1)+l
                   kl2=ndim*(l-1)+k
                   if(ij1.ge.kl1)then
         if(i.gt.ntash.or.j.gt.ntash.or.k.gt.ntash.or.l.gt.ntash)then
                       tmpP(itri(ij1,kl1))=0.0d0
                     else
                       tmpP(itri(ij1,kl1))=tmpPM(i,j,k,l)
                     end if
                   end if
                 end do
               end do
             end do
           end do

           ij=0
           do i1=1,ndim
             do j1=1,ndim
               ij=ij+1
               D_CI(ij)=tmpDe(i1,j1)
             end do
           end do
           do i=1,n2dens
             P_CI(i)=tmpP(i)
           end do
           call mma_deallocate(tmpDe)
           call mma_deallocate(tmpDeM)
           call mma_deallocate(tmpP)
           call mma_deallocate(tmpPM)
           call dmrg_dim_change_mclr(RGras2(1:8),
     &                               ntash,0)
         end if
! ===================================================================

*
*-------Some administrative shit
*
*       Store densities in triangular form
*
         Do i=1,ntAsh
          Do j=1,i
           D1(itri(i,j))=D_CI((i-1)*ntash+j)
          End Do
         End Do

         Do i=1,ntAsh
          Do j=1,i
           ij=itri(i,j)
           ij2=i+(j-1)*ntash
           ji2=j+(i-1)*ntash
           Do k=1,i
            Do l=1,k
             kl=itri(k,l)
             kl2=k+(l-1)*ntash
             lk2=l+(k-1)*ntash
             ijkl=itri(ij2,kl2)
             jikl=itri(ji2,kl2)
             ijlk=itri(ij2,lk2)
             jilk=itri(ji2,lk2)
             P1(itri(ij,kl))=Quart*(P_CI(ijkl)+P_CI(jikl)+
     &                              P_CI(ijlk)+P_CI(jilk))
            End Do
           End Do
          End Do
         End Do

        DO  K=1,NTASH
         DO L = 1, K
          KL = K*(K-1)/2 + L
          KLROW = KL*(KL-1)/2
          IF( L .EQ. K ) THEN
           IMAX = K
          ELSE
           IMAX = K-1
          END IF
          DO I = 1,IMAX
           II= I*(I+1)/2
           IIKL= KLROW + II
           P1(IIKL) = P1(IIKL)*Half
          End Do
         End Do
         End Do


C         Do i=1,ntAsh
C         Do j=1,i
C         ij=itri(i,j)
C         ij2=i+(j-1)*ntash
C         ji2=j+(i-1)*ntash
C         Do k=1,ntAsh
C         Do l=1,k
C          kl=itri(k,l)
C          kl2=k+(l-1)*ntash
C          ijkl=itri(ij2,kl2)
C          jikl=itri(ji2,kl2)
C          fact=Half
C          if(ij.ge.kl .and. k.eq.l) fact=Quart
C          if(ij.lt.kl .and. i.eq.j) fact=Quart
C          P1(itri(ij,kl))=
C     &        fact*(P_CI(ijkl)+P_CI(jikl))
C         End Do
C         End Do
C         End Do
C         End Do
C         If (debug) Call triprt('P1',' ',P1,(ntash**2+ntash)/2)
c
c Write the 'bar' densities to disk,  not symmetry blocked.
c

!         Call Put_dArray('DLMO',D1,ndim1) ! \bar{D} triangular  ! yma
!         Call Put_dArray('PLMO',P1,ndim2) ! \bar{d} triangular  ! yma

         Call Put_dArray('DLMO',D1,nDLMO) ! \bar{D} triangular
         Call Put_dArray('PLMO',P1,nPLMO) ! \bar{d} triangular
*
       End If
*
*      2) Orbital response
*         ================
*
*       Read in from disk
*
       iDisk=iKapDisp(1)
       Call dDaFile(LuTemp,2,K1,nDensC,iDisk) ! Read \bar{kappa}
       Call Uncompress(K1,K2,1)
c
c If we want to estimate the error
c
       If (esterr) Then
*        Do iestate=1,lroots
*          Call calcerr(K2,iestate)
*        End do
         Call calcerr(K2,istate)
       End If
*
*----- First we fix the renormalization contribution
*
       Call mma_allocate(D_K,nLCMO,Label='D_K')
       Call Get_dArray_chk('FockOcc',D_K,nLCMO)
*      Calculates the effective Fock matrix
       Call Make_Conn(Conn,K2,P_CI,D_CI)   !D_CI not changed
       Call DaxPy_(ndens2,One,D_K,1,Conn,1)
*      call dcopy_(ndens2,D_K,1,Conn,1)
       If (PT2) Then
         !! Add the WLag term (will be contracted with overlap
         !! derivative) computed in CASPT2
         Do i = 1, nTot1
           Read(LuPT2,*) Val
           Conn(i) = Conn(i) + Val
         End Do
       End If
       Call Put_dArray('FockOcc',Conn,nTot1)
*
*      Transposed one index transformation of the density
*      (only the inactive density to store it separately)
*
       Call OITD(K2,1,DAO,Dtmp,.False.)
*
       If (PT2) Then
         !! For gradient calculation. D^var couples with inactive
         !! orbitals only, so D0(1,4) (in integral_util/prepp.f), which
         !! couples with active orbitals has to be modified,
         Do iSym = 1, nSym
           nBasI = nBas(iSym)
           Do iI = 1, nBasI
             Do iJ = 1, nBasI
               Read(LuPT2,*) Val
               DAO(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *       = DAO(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *       + Val
             End Do
           End Do
         End Do
         !! The PT2 density will be used later again.
         Do iSym = 1, nSym
           nBasI = nBas(iSym)
           Do iI = 1, nBasI
             Do iJ = 1, nBasI
               Backspace LuPT2
             End Do
           End Do
         End Do
       End If
*
*      Transformation to AO basis (covariant)
*
c
c Transforms to AO differently dep on last arg.
c
       Call TCMO(DAO,1,-2)
*
*      Fold AO density and write to disk
c Mult all terms that are not diag by 2
*
       Call FOLD2(nsym,nbas,DAO,K1)
*
       Call Put_dArray('DLAO',K1,ntot1)
*
*      Now with active density too, to form the variational density
*
!      gives \tilde{D}
       Call OITD(K2,1,D_K,Dtmp,.True.)
*
       Do iS=1,nsym
c
c C*\tilde{\kappa} --> ipDAO
c
          If (nBas(is).ge.1)
     &       CALL DGEMM_('N','N',
     &                   NBAS(is),NBAS(is),NBAS(is),
     &                   One,CMO(ipCM(is)),NBAS(is),
     &                   K2(ipmat(is,is)),NBAS(is),
     &                   Zero,DAO(ipCM(is)),NBAS(is))
       End Do
*
       Call Put_dArray('LCMO',DAO,nLCMO)
*
       if(doDMRG)then  ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),ntash,0)
         call dmrg_spc_change_mclr(RGras2(1:8),nash)
       end if
*
       If (isNAC) Then
         ng1=nNAC
         Call mma_allocate(G1q,ng1,Label='G1q')
         Call Get_dArray_chk('D1mo',G1q,ng1)
         iR = 0 ! set to dummy value.
       Else
         iR=iroot(istate)
         jdisk=itoc(3)
         ng1=itri(ntash,ntash)
         ng2=itri(ng1,ng1)
         Call mma_allocate(G1q,n1dens,Label='G1q')
c
c Read active one el dens for state j from JOBIPH and store in G1q
c
         Call Get_cArray('Relax Method',Method,8)
         if(Method.eq.'MSPDFT  ') then
          Call Get_DArray('D1MOt           ',G1q,ng1)
         else
           Do i=1,iR-1  ! Dummy read until state j
             Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
             Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
             Call dDaFile(LUJOB ,0,rdum,ng2,jDisk)
             Call dDaFile(LUJOB ,0,rdum,ng2,jDisk)
           End Do
           Call dDaFile(LUJOB ,2,G1q,ng1,jDisk)

           If (PT2.and.nRoots.gt.1) Call Get_dArray("D1mo",G1q,ng1)

         end if
       EndIf
*
*    Construct a variationally stable density matrix. In MO
c
c D_eff = D^j + \tilde{D} +\bar{D}
c D_K = (G1q + inact) + D_K + D_CI
*
C
       If (isNAC) Then
*
** For NAC, first build DAO and then DAO_var
*
         Do is=1,nSym
c Note: no inactive part for transition densities
          Do iA=1,nash(is)
           Do jA=1,nash(is)
            i=iA+nish(is)
            j=jA+nish(is)
            iAA=iA+na(is)
            jAA=jA+na(is)
            D_K(ipmat(is,is)+i-1+(j-1)*nbas(is))=
     &       D_K(ipmat(is,is)+i-1+(j-1)*nbas(is))
     &      +D_CI(iAA+(jAA-1)*ntash)
     &      +G1q(itri(iAA,jAA))
           End Do
          End Do
         End Do
C
         If (PT2) Then
           !! PT2 density (in MO)
           Do iSym = 1, nSym
             nBasI = nBas(iSym)
             Do iI = 1, nBasI
               Do iJ = 1, nBasI
                 Read(LuPT2,*) Val
                 D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *         = D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *         + Val
               End Do
             End Do
           End Do
           !! PT2C density (in MO)
           Do iSym = 1, nSym
             nBasI = nBas(iSym)
             Do iI = 1, nBasI
               Do iJ = 1, nBasI
                 Read(LuPT2,*) Val
                 D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *         = D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *         + Val*0.25d+00
                 D_K(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI)
     *         = D_K(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI)
     *         + Val*0.25d+00
               End Do
             End Do
           End Do
         End If
         Call mma_allocate(Temp,nBuf/2,Label='Temp')
         Call NatOrb(D_K,CMO,CMON,OCCU)
         Call dmat_MCLR(CMON,OCCU,Temp)
         Call Put_dArray('D1aoVar',Temp,nTot1)
         Call mma_deallocate(Temp)
*
** Transform the antisymmetric transition density matrix to AO
**  (there is no guarantee the symmetry will work here)
*
         iDisk=0
         LuDens=20
         Call DaName(LuDens,'MCLRDENS')
         Call dDaFile(LuDens,2,G1q,ng1,iDisk)
         Call DaClos(LuDens)
         Call mma_allocate(G1m,ndens2,Label='G1m')
         G1m(:)=Zero
* Reconstruct the square matrix
         Do is=1,nSym
          Do iA=1,nash(is)
           i=iA+nish(is)
           iAA=iA+na(is)
           Do jA=1,iA-1
            j=jA+nish(is)
            jAA=jA+na(is)
            G1m(ipmat(is,is)+i-1+(j-1)*nbas(is))=
     &           G1q(itri(iAA,jAA))
            G1m(ipmat(is,is)+j-1+(i-1)*nbas(is))=
     &          -G1q(itri(iAA,jAA))
           End Do
           G1m(ipmat(is,is)+i-1+(i-1)*nbas(is))=Zero
          End Do
         End Do
*
         If (PT2) Then
           !! PT2C density (in MO) for CSF derivative
           Do iSym = 1, nSym
             nBasI = nBas(iSym)
             Do iI = 1, nBasI
               Do iJ = 1, nBasI
                 Read(LuPT2,*) Val
                 G1m(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI)
     *         = G1m(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI) + Val
               End Do
             End Do
           End Do
         End If
* Transform
         Call TCMO(G1m,1,-2)
* Save the triangular form
         iOff=0
         Do is=1,nSym
          ibas=nbas(is)
          Do i=1,ibas
           Do j=1,i
            G1m(iOff+itri(i,j))=
     &        G1m(ipmat(is,is)+j-1+(i-1)*nbas(is))
           End Do
          End Do
          iOff=iOff+(ibas*ibas+ibas)/2
         End Do
         Call Put_dArray('D1ao-',G1m,nTot1)
         Call mma_deallocate(G1m)
*
       Else
*
** Normal SA gradient (no NAC)
*
         Do is=1,nSym
          Do i=1,nish(is)
c
c The inactive density
c
           D_K(ipmat(is,is)+i-1+(i-1)*nbas(is))=
     &     D_K(ipmat(is,is)+i-1+(i-1)*nbas(is))+Two
          End DO
          Do iA=1,nash(is)
           Do jA=1,nash(is)
            i=iA+nish(is)
            j=jA+nish(is)
            iAA=iA+na(is)
            jAA=jA+na(is)
c
c The active density G1q and \bar{D}
c
            D_K(ipmat(is,is)+i-1+(j-1)*nbas(is))=
     &       D_K(ipmat(is,is)+i-1+(j-1)*nbas(is))
     &      +D_CI(iAA+(jAA-1)*ntash)
     &      +G1q(itri(iAA,jAA))
           End Do
          End Do
         End Do
C
         If (PT2) Then
           !! Add PT2 density (in MO)
           Do iSym = 1, nSym
             nBasI = nBas(iSym)
             Do iI = 1, nBasI
               Do iJ = 1, nBasI
                 Read(LuPT2,*) Val
                 D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *         = D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *         + Val
               End Do
             End Do
           End Do
           !! Also, PT2C density (in MO)
           !! This density couples with inactive orbitals only,
           !! while the above PT2 density couples with inactive+active
           !! orbitals.
           Do iSym = 1, nSym
             nBasI = nBas(iSym)
             Do iI = 1, nBasI
               Do iJ = 1, nBasI
                 Read(LuPT2,*) Val
                 D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *         = D_K(ipMat(iSym,iSym)+iI-1+(iJ-1)*nBasI)
     *         + Val*0.25d+00
                 D_K(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI)
     *         = D_K(ipMat(iSym,iSym)+iJ-1+(iI-1)*nBasI)
     *         + Val*0.25d+00
               End Do
             End Do
           End Do
         End If
c
c Diagonalize the effective density to be able to use Prpt
c OCCU eigenvalues of eff dens
c CMON eigenvectors (new orb coef)
c
         Call NatOrb(D_K,CMO,CMON,OCCU)
         Call mma_Allocate(Tmp,nBuf/2,Label='Tmp')
         Call dmat_MCLR(CMON,OCCU,Tmp)
         Call Put_dArray('D1aoVar',Tmp,nTot1)
         Call mma_deallocate(Tmp)

         Call mma_allocate(TEMP,nNac,Label='TEMP')
         Call mma_allocate(tTmp,nNac,Label='tTmp')
         Call get_dArray_chk('D1mo',TEMP,nNac)
         Call get_dArray_chk('DLMO',tTmp,nNac)
         Call DaxPy_(nNac,1.0d0,tTmp,1,TEMP,1)
         Call mma_deallocate(TEMP)
         Call mma_deallocate(tTmp)

         Note='var'
         LuTmp=50
         LuTmp=IsFreeUnit(LuTmp)
         Call WrVec('TMPORB',LuTmp,'O',nSym,nBas,nBas,
     &            rDum,OCCU,rDum,iDum,Note)
         Call Prpt()
*                                                                      *
************************************************************************
*        There should now be dipole moments on the runfile which
*        corresponds to the gradient of the energy w.r.t. the
*        electric field. Let's update the list of values stored
*        on the runfile.
*
         Is_Roots_Set = .False.
         Call Qpg_iScalar('Number of roots',Is_Roots_Set)
         nRoots = 1
         If (Is_Roots_Set) Then
            Call Get_iScalar('Number of roots',nRoots)
         End If
*
         If (nRoots.ne.1) Then
*           Write (*,*) 'iR=',iR
            Call mma_allocate(DM,3,Label='DM')
            Call mma_allocate(DMs,3*nROOTS,Label='DMs')
            Call Get_dArray('Last Dipole Moments',DMs,3*nRoots)
*           Call RecPrt('Last Dipole Moments',' ',DMS,3,nRoots)
            Call Get_dArray('Dipole Moment',DM,3)
*           Call RecPrt('Dipole Moment',' ',DM,1,3)
            Call DCopy_(3,DM,1,DMS(1+(iR-1)*3),1)
*           Call RecPrt('Last Dipole Moments',' ',DMS,3,nRoots)
            Call Put_dArray('Last Dipole Moments',DMs,3*nRoots)
            Call mma_deallocate(DMs)
            Call mma_deallocate(DM)
         End If
************************************************************************
*                                                                      *
       End If
       Call mma_deallocate(G1q)
C
c--------------------------  debug -----
c
       if(doDMRG)then ! yma
         call dmrg_dim_change_mclr(LRras2(1:8),ntash,0)
         call dmrg_spc_change_mclr(LRras2(1:8),nash)
       end if

c
c  Write the effective active one el density to disk in the same format as g1q
c
c       Call mma_allocate(Deff_act,ndens2,Label='Deff_act')
c       call dcopy_(nDens2,D_K,1,Deff_act,1)
c       Do is=1,nSym
c        Do i=1,nish(is)
c
c Subtract the inactive density
c
c         Deff_act(ipmat(is,is)+i-1+(i-1)*nbas(is))=
c     &   D_K(ipmat(is,is)+i-1+(i-1)*nbas(is))-Two
c        End Do
c       End Do
c
c      Call Put_DEff(Deff_act,ndens2)
c
c       Call mma_deallocate(Deff_act)
c
c--------------------------------------------------
c
c Diagonalize the effective density to be able to use Prpt
c OCCU eigenvalues of eff dens
c CMON eigenvectors (new orb coef)
c
c      Call NatOrb(D_K,CMO,CMON,OCCU)
c      Call mma_allocate(Temp,nBuf/2,Label='Temp')
c      Call dmat_MCLR(CMON,OCCU,Temp)
c      Call Put_dArray('D1aoVar',Temp,nTot1)
c      Note='var'
c      LuTmp=50
c      LuTmp=IsFreeUnit(LuTmp)
c      Call WrVec('TMPORB',LuTmp,'O',nSym,nBas,nBas,
c    &            Dum,OCCU,Dum,iDum,Note)
c      Call Prpt()

c
c Standard routine, Temp effective dens in AO
c
*       Call dmat_MCLR(CMON,OCCU,Temp)
c
*       Call Put_dArray('D1aoVar',Temp,nTot1)
c      Call mma_deallocate(Temp)

       Call Put_iScalar('SA ready',1)
       If (isNAC) Then
         Write(mstate,'(1X,I7,",",I7)') NACStates(1),NACStates(2)
       Else
         Write(mstate,'(I16)') irlxroot
       End If
       If (override) mstate(1:1)='+'
       Call Put_cArray('MCLR Root',mstate,16)
*
       Call mma_deallocate(K1)
       Call mma_deallocate(K2)
       Call mma_deallocate(DAO)
       Call mma_deallocate(D_CI)
       Call mma_deallocate(D1)
       Call mma_deallocate(P_CI)
       Call mma_deallocate(P1)
       Call mma_deallocate(Conn)
       Call mma_deallocate(OCCU)
       Call mma_deallocate(CMON)
       Call mma_deallocate(Dtmp)
       Call mma_deallocate(D_K)

       irc=ipclose(-1)
*
       Return
#ifdef _WARNING_WORKAROUND_
       If (.False.) Call Unused_integer(irc)
#endif
       End

c --------------------------------------------------------------------------
c
      Subroutine OITD(rK,isym,D,Dtmp,act)
*
      use Arrays, only: G1t
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "real.fh"
#include "sa.fh"
      Real*8 rK(*),D(*),Dtmp(*)
      Logical act
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
      call dcopy_(ndens2,[Zero],0,Dtmp,1)
*
*     Note: even with NAC we set the inactive block,
*     because this is the SA density, not the transition density
      Do iS=1,nSym
        Do iB=1,nIsh(iS)
          Dtmp(1+(ipCM(iS)+(ib-1)*nOrb(iS)+ib-1)-1) = Two
        End Do
      End Do
      If (act) Then
       Do iS=1,nSym
        Do iB=1,nAsh(iS)
         Do jB=1,nAsh(iS)
          Dtmp(1+(ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nOrb(is)-1)-1)=
     &    G1t((itri((nA(is)+ib),(nA(is)+jb))))
         End Do
        End Do
       End Do
      End If
*
      Do iS=1,nsym
         jS=ieor(iS-1,isym-1)+1
         If (nOrb(iS)*nOrb(jS).ge.1) Then
            Call DGEMM_('N','T',nOrb(iS),nOrb(jS),nOrb(iS),One,
     &                 Dtmp(1+ipCM(iS)-1),nOrb(iS),
     &                 rK(ipMat(jS,iS)),nOrb(jS),
     &                 Zero,D(ipMat(iS,jS)),nOrb(iS))
            Call DGEMM_('T','N',nOrb(iS),nOrb(jS),nOrb(jS),-One,
     &                 rK(ipMat(jS,iS)),nOrb(jS),
     &                 Dtmp(1+ipCM(jS)-1),nOrb(jS),
     &                 One,D(ipMat(iS,jS)),nOrb(iS))
         End If
      End Do
      Return
      End
      Subroutine NatOrb(Dens,CMOO,CMON,OCCN)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "stdalloc.fh"
#include "real.fh"
      Real*8 Dens(*),CMOO(*),CMON(*),OCCN(*)
      Real*8, Allocatable:: EVal(:), EVec(:)

      Call mma_allocate(EVec,ndens2,Label='EVec')
      Call mma_allocate(EVal,ndens2,Label='EVal')
C
C         Diagonalize the density matrix and transform orbitals
C
      If (iAnd(kprint,8).eq.8) Then
         Write(6,*)
         Write(6,*) '           Effective natural population '
         Write(6,*) '           ============================ '
         Write(6,*)
      End If
      io=0
      Do is=1,nsym
         ij=0
         Do i=0,nbas(is)-1
            Do j=0,i
               ij=ij+1
               Eval(ij)=Dens(ipMat(is,is)+i+j*nbas(is))
            End DO
         End DO
         EVec(:)=Zero
         Call dCopy_(nBas(is),[One],0,EVec,nbas(is)+1)
         CALL JACOB(EVal,EVec,nbas(is),nbas(is))
         ii=0
         DO i=1,nbas(is)
            ii=ii+i
            OCCN(io+i)=Eval(ii)
         END DO
         IST=IO+1
         IEND=IO+NBAS(is)
         If (iAnd(kprint,2).eq.2)
     &      Write (6,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))')
     &             'sym',iS,':',(OCCN(I),I=IST,IEND)
         If (nBas(is).ge.1)
     &      CALL DGEMM_('N','N',
     &                  NBAS(is),NBAS(is),NBAS(is),
     &                  One,CMOO(ipCM(is)),NBAS(is),
     &                  EVec,NBAS(is),
     &                  Zero,CMON(ipCM(is)),NBAS(is))
         io=io+nbas(is)
      End DO

      Call mma_deallocate(EVec)
      Call mma_deallocate(Eval)

      Return
      End
