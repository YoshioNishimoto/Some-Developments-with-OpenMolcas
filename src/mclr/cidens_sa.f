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
      SubRoutine CIDens_sa(RSP,iLS,iRS,iL,iR,rP,rD)
      use ipPage, only: W
      use pcm_grad, only: do_RF,DZACTMO,iStpPCM
      use ISRotation, only: InvSCF,ISR
      Implicit Real*8(a-h,o-z)
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "stdalloc.fh"
#include "crun_mclr.fh"

#include "Input.fh"
#include "Pointers.fh"
#include "spinfo_mclr.fh"
#include "cands.fh"
#include "dmrginfo_mclr.fh"
#include "sa.fh"
      Real*8 rP(*),rD(*)
      Logical RSP
      integer opout
      Real*8, Allocatable:: De(:), Pe(:), CIL(:), CIR(:)

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

* LS = CI
*
*     Ok, we once more want to hide Jeppe's routines from
*     the eyes of the world, so everyone believes that I have done
*     all the work.
*     If we have spin dependent Hamiltonian we will work with
*     SD in all parts of the program, no CSF is necessary,
*     otherwise we will do the optimization in CSF's to increase
*     convergence.
*
*     Input:
*
*     RSP: true if the density should be used for response calculations
*     iLS: CI Coeff for left state
*     iRS: CI Coeff for right state
*     iL : Symmetry of left state
*     iR : Symmetry of right state
*
*
*               +       +
*     iS=1 E  =a  a  + a  a
*           pq  ap aq   Bp Bq
*
*                +       +
*     iS=-1 T  =a  a  - a a
*            pq  ap aq   Bp Bq
*
*     Output:
*
*      rP : Two Electron Density
*      rD : One Electron Density
*
*
      If (nconf1.eq.0) return

      if(doDMRG)then
        call dmrg_dim_change_mclr(LRras2(1:8),ndim,0)
        call dmrg_dim_change_mclr(LRras2(1:8),nna,0)
        n1dens=ndim**2
        n2dens=n1dens*(n1dens+1)/2
      end if

      Call mma_allocate(De,n1dens,Label='De')
      Call mma_allocate(Pe,n2dens,Label='Pe')
      call dcopy_(n1dens,[0.0d0],0,rD,1)
      call dcopy_(n2dens,[0.0d0],0,rP,1)
      if (do_RF) call dcopy_(n1dens,[0.0d+00],0,DZACTMO,1)

*
      nConfL=Max(ncsf(il),nint(xispsm(il,1)))
      nConfR=Max(ncsf(iR),nint(xispsm(iR,1)))

      Call mma_allocate(CIL,nConfL,Label='CIL')
      Call mma_allocate(CIR,nConfR,Label='CIR')
*
      Do i=0,nroots-1
        irc=ipin(iLS)
        irc=ipin(iRS)
        Call CSF2SD(W(iLS)%Vec(1+i*ncsf(il)),CIL,iL)
        irc=opout(iLS)
        Call CSF2SD(W(iRS)%Vec(1+i*ncsf(ir)),CIR,iR)
        irc=opout(iRS)
        irc=ipnout(-1)
        icsm=iR
        issm=iL
        Call Densi2(2,De,Pe,CIL,CIR,0,0,0,n1dens,n2dens)

        If (RSP) Then
           Do iA=1,nnA
             Do jA=1,nnA
               Do kA=1,nnA
                Do la=1,nnA
                 ij1=nnA*(iA-1)+ja
                 ij2=nna*(ja-1)+ia
                 kl1=nnA*(ka-1)+la
                 kl2=nna*(la-1)+ka
                 if (ij1.ge.kl1)
     &           rp(itri(ij1,kl1))=rp(itri(ij1,kl1))+weight(1+i)*(
     &            Pe(itri(ij1,kl1))+Pe(itri(ij2,kl2)))
                End Do
               End Do
             End Do
           End Do
           Do iA=1,nnA
              Do jA=1,nnA
                 ij1=nnA*(iA-1)+ja
                 ij2=nna*(ja-1)+ia
                 rD(ij1)=rD(ij1)+
     *                weight(1+i)*(De(ij1)+De(ij2))
              End Do
           End Do
        Else
           call daxpy_(n2dens,weight(i+1),Pe,1,rp,1)
           call daxpy_(n1dens,Weight(i+1),De,1,rD,1)
        End If
*
        !! rotations between internal states
        if (.not.InvSCF) then
          if (iStpPCM==2) then
            irc=ipin(iLS)
            Call CSF2SD(W(iLS)%Vec(1+i*ncsf(il)),CIL,iL)
            irc=opout(iLS)
          else if (iStpPCM==3) then
            irc=ipin(iRS)
            Call CSF2SD(W(iRS)%Vec(1+i*ncsf(ir)),CIL,iR)
            irc=opout(iRS)
          end if
          do j = 0, i
            scal = ISR%p(i+1,j+1)
C           write (6,*) i+1,j+1,scal
            if (abs(scal).le.1.0d-10) cycle
            if (iStpPCM==2) then
              irc=ipin(iLS)
              Call CSF2SD(W(iLS)%Vec(1+j*ncsf(il)),CIR,iL)
              irc=opout(iLS)
            else if (iStpPCM==3) then
              irc=ipin(iRS)
              Call CSF2SD(W(iRS)%Vec(1+j*ncsf(ir)),CIR,iR)
              irc=opout(iRS)
            end if
            irc=ipnout(-1)
            icsm=iR
            issm=iL
            Call Densi2(2,De,Pe,CIL,CIR,0,0,0,n1dens,n2dens)
            call dscal_(n1dens,scal,De,1)
            call dscal_(n2dens,scal,Pe,1)
            If (RSP) Then
               Do iA=1,nnA
                 Do jA=1,nnA
                   Do kA=1,nnA
                    Do la=1,nnA
                     ij1=nnA*(iA-1)+ja
                     ij2=nna*(ja-1)+ia
                     kl1=nnA*(ka-1)+la
                     kl2=nna*(la-1)+ka
                     if (ij1.ge.kl1)
     &               rp(itri(ij1,kl1))=rp(itri(ij1,kl1))+
     &                Pe(itri(ij1,kl1))+Pe(itri(ij2,kl2))
                    End Do
                   End Do
                 End Do
               End Do
               Do iA=1,nnA
                  Do jA=1,nnA
                     ij1=nnA*(iA-1)+ja
                     ij2=nna*(ja-1)+ia
                     rD(ij1)=rD(ij1) + De(ij1)+De(ij2)
                  End Do
               End Do
            Else
               call daxpy_(n2dens,1.0d+00,Pe,1,rp,1)
               call daxpy_(n1dens,1.0d+00,De,1,rD,1)
            End If
          end do
        end if
      End Do
*
      !! PCM: for implicit contributions
      if (do_RF) call dcopy_(n1dens,rD,1,DZACTMO,1)
C     write (*,*) "asdf"
C     call sqprt(rd,nna)
C      write (*,*) "g2r"
C      do i = 1, itri(ntash**2,ntash**2)
C      write (*,'(i4,f20.10)') i,rp(i)*2.0d+00
C      end do
C      call abend
C     call abend
*
      Call mma_deallocate(CIL)
      Call mma_deallocate(CIR)
      Call mma_deallocate(Pe)
      Call mma_deallocate(De)

      if(doDMRG)then
        call dmrg_dim_change_mclr(RGras2(1:8),ndim,0)
        call dmrg_dim_change_mclr(RGras2(1:8),nna,0)
        n1dens=ndim**2
        n2dens=n1dens*(n1dens+1)/2
      end if
*
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Call Unused_integer(irc)
#endif
      End
