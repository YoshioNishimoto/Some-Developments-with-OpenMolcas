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
      SubRoutine DMInvKap(rMFact,rIn,nrIn,rOut,nrOut,rtemp,nrtemp,
     &                    isym,iter)
************************************************************************
*                                                                      *
*     _____     -1                                                     *
*     Kappa  = M  Kappa                                                *
*          ip   pq     iq                                              *
*                                                                      *
*                                                                      *
*     In: rMFact        Factorized preconditioner (diagonal part       *
*                         of the electronic hessian that couple          *
*                         rotations with one common index)               *
*     In,Out rOut       Orbital rotaotion                              *
*                                                                      *
*     iSym              Symmetry of rotation                           *
*                                                                      *
************************************************************************

      Implicit Real*8 (a-h,o-z)
#include "Input.fh"
#include "standard_iounits.fh"
#include "Pointers.fh"
#include "dmrginfo_mclr.fh"
      Real*8 rOut(nrOut),rMFact(*),rIn(nrIn),rtemp(nrTemp)
*                                                                      *
************************************************************************
*                                                                      *
      Call DMInvKap_Internal(rMFact)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine DMInvKap_Internal(rMFact)
      Use Iso_C_Binding
      Real*8, Target :: rMFact(*)
      Integer, Pointer :: iMFact(:)
      ip1=1
*
      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nash)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if

C     If (ActRot) Call DCopy_(nDensC,rIn,1,rOut,1)
      Call DCopy_(nDensC,rIn,1,rOut,1)
      Call Uncompress2(rIn,rtemp,isym)
C     write (*,*) "norb= ", norb(1)
C     write (*,*) "uncompressed"
C     call sqprt(rtemp,norb(1))
C     call abend
*                                                                      *
************************************************************************
*                                                                      *
      Do jS=1,nSym
         iS=iEOr(js-1,iSym-1)+1
*                                                                      *
************************************************************************
*                                                                      *
*        kappa
*            ip
*                                                                      *
************************************************************************
*                                                                      *
         Do iI=1,nIsh(js)
            nD=nOrb(is)-nIsh(is)
            If (nd.ne.0) Then
               ip2=ipMat(is,js)+nOrb(is)*(iI-1)
               if (actrot) ip2 = ip2 + nIsh(iS) !! jS?
               irc=0
               call c_f_pointer(c_loc(rMFact(ip1+nd**2)),iMFact,[ND])
               call dgetrs_('N',ND,1,rMFact(ip1),nd,
     &                       iMFact,rtemp(ip2),nd,irc)
               nullify(iMFact)
               If (irc.ne.0) then
                   Write(6,*) 'Error in DGETRS called from dminvkap'
                   Call Abend
               endif
               ip1=ip1+nD*(nD+1)
            End If
         End Do

*                                                                      *
************************************************************************
*                                                                      *
*        kappa
*             ap
*                                                                      *
************************************************************************
*                                                                      *
         Do iI=1,nAsh(js)
            If (ActRot) Then
               nD=nOrb(is)
            Else
C              nD=nOrb(is)-nAsh(is)
               If (iI.le.nRs1(jS)) THen
                 nD=nOrb(is)-nRs1(js)
               Else If (iI.le.nRs1(jS)+nRs2(jS)) Then
                 nD=nOrb(is)-nRs2(js)
               Else If (iI.le.nRs1(jS)+nRs2(jS)+nRs3(jS)) Then
                 nD=nOrb(is)-nRs3(js)
               End If
            End If
            If (nd.ne.0) Then
               ip2=ipMat(is,js)+nOrb(is)*(iI-1+nIsh(js))
               irc=0
               call c_f_pointer(c_loc(rMFact(ip1+nd**2)),iMFact,[ND])
               call dgetrs_('N',ND,1,rMFact(ip1),nd,
     &                      iMFact,rtemp(ip2),nd,irc)
               nullify(iMFact)
               If (irc.ne.0) then
                   Write(6,*) 'Error in DGETRS called from dminvkap'
                   Call Abend
               endif
               ip1=ip1+nD*(nD+1)
            End If
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*                                                                      *
************************************************************************
*
C     write (*,*) "rtemp"
C     do i = 1, nrtemp
C       write (*,'(i3,f20.10)') i,rtemp(i)
C     end do
C     call sqprt(rtemp,12)
C     write (*,*) "before compess2"
C     do i = 1, nrout
C       write (*,'(i3,f20.10)') i,rout(i)
C     end do
C       do i = 1, nrtemp
C         if (abs(rtemp(i)).le.1.0d-10) rtemp(i)=0.0d+00
C       end do
      Call Compress2(rtemp,nrtemp,rOut,nrOut,isym)
C     write (*,*) "after compess2"
C     do i = 1, nrout
C       write (*,'(i3,f20.10)') i,rout(i)
C     end do
C     call abend

      if(doDMRG)then
        call dmrg_spc_change_mclr(LRras2(1:8),nash)
      end if
*                                                                      *
************************************************************************
*                                                                      *
*     Warn if the trail vector becomes large
*
      If ((ddot_(ndensc,rout,1,rout,1).gt.100.D0).and.(iter.eq.1)) Then
          write(LuWr,*)'****************************************'
          write(LuWr,*)'*                                      *'
          write(LuWr,*)'*           WARNING!!                  *'
          write(LuWr,*)'* Elements in the E^[2] matrix small!! *'
          write(LuWr,*)'* The calculation might diverge.       *'
          write(LuWr,*)'*                                      *'
          write(LuWr,*)'* Check your active space!!!!          *'
          write(LuWr,*)'*                                      *'
          write(LuWr,*)'* Make sure degenerate orbitals do not *'
          write(LuWr,*)'* belong to different spaces.          *'
          write(LuWr,*)'* Note that no LR code can handle      *'
          write(LuWr,*)'* 2.0d0 occupancy in active orbitals!! *'
          write(LuWr,*)'****************************************'
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine DMInvKap_Internal
*
      end
