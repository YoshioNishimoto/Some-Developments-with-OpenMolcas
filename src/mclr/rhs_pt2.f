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
* Copyright (C) 1998, Anders Bernhardsson                              *
************************************************************************
      Subroutine RHS_PT2(rkappa,iprci,CLag,SLag)
*
*    Calculate and reads in from CASPT2 the RHS needed
*    for the calculation of the Lagrangian multipliers
*    for CASPT2.
*
*    From the same Fock matrix the effective Fock matrix
*    needed for the renormalization term is  calculated.

*    Here is the structure, it is not debugged and one
*    needs to check the detail, but all "input" is there.
*
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
#include "glbbas_mclr.fh"
#include "Files_mclr.fh"
#ifndef _DEBUGED_
#include "detdim.fh"
#include "csfbas_mclr.fh"
#endif
C     #include "caspt2.fh"
      Real*8 rKappa(*),CLag(*),SLag(*)
      Half=0.5d0

*
*     Read in a and b part of effective gradient from CASPT2
*
C     write (*,*) "Lugamma = ", lugamma
      nOLag = 0
      nCLag = 0
      DO i = 1, nSym
        nOLag = nOLag + nOrb(i)*nOrb(i)
        nCLag = nCLag + nRoots*nCSF(i)
      END DO
      nSLag = nRoots*(nRoots-1)/2
C     Call Get_dArray('CLAG',Work(ipCLag),nClag)
C          LUTMP = 88
           Do i = 1, nCLag
             Read (LUPT2,*) CLag(i)
           End Do
           Do i = 1, nOLag
             Read (LUPT2,*) tmp ! rKappa(i)
             rKappa(i) = rKappa(i) + tmp
           End Do
           Do i = 1, nSLag
             Read (LUPT2,*) SLag(i)
           End Do
C       call sqprt(rkappa,12)
C       call dscal(nclag,-2.00d+00,clag,1)
C       call dcopy(nclag,0.0d+00,0,clag,1)
C       call dcopy(nolag,0.0d+00,0,rKappa,1)
C       call dcopy(nslag,0.0d+00,0,slag,1)
C       call dcopy(nclag,clag,1,clag2,1)
C     Call Get_dArray('SLAG',Work(ipSLag),nSlag)

      return
      write (*,*) "in rhs_pt2"
      write (*,*) "imethod =", imethod
      If (imethod.eq.2)Then
      i=0
      Call Getmem('TEMPCI','ALLO','REAL',ipT,nconf1)
      Call Getmem('TEMPCI','ALLO','REAL',ipT2,nconf1)
      Call Getmem('TEMPKAP','ALLO','REAL',ipK,ndens2)
      !! Call dDaFile(LuPT2,2,Work(ipK),ndens2,i)
#ifdef _DEBUGED_
*
*    CASPT2 mode
*
      !! Call dDaFile(LuPT2,2,Work(ipT),nconf1,i)
#else
*
* lucia mode
*
      !! n=1 !! nint(xispsm(State_SYM,1))
      !! Call dDaFile(LuPT2,2,Work(ipT),n,i)
#endif
*
*---- Transform from split GUGA to symmetric group
*
#ifdef _DEBUGED_
*
*    CASPT2 mode (split graph)
*
      !! Call Gugactl_MCLR(ipT,1)
#else
*
* lucia mode (Symmetric group)
*
      !!iprdia=0
      !!Call INCSFSD(STATE_SYM,STATE_SYM,.false.)
      !!CALL CSDTVC_MCLR(Work(ipT2),Work(ipT),2,
     &!!                 WORK(KDTOC),iWORK(KICTS(1)),
     &!!                 State_SYM,1,IPRDIA)
#endif
*
      Else
*
*     MP2
*
        i=0
        Call Getmem('TEMPKAP','ALLO','REAL',ipK,ndens2)
        Call dDaFile(LuPT2,2,Work(ipK),ndens2,i)
*       Call ThreeP(Kappa)
*       MP2
      End if
* ----
*
      Call GetMem('Temp1','ALLO','REAL',ipT1,ndens2)
      Call GetMem('Temp2','ALLO','REAL',ipT2,ndens2)
      Call GetMem('FockAO1','ALLO','REAL',ipFAO1,ndens2)
      Call GetMem('FockAO2','ALLO','REAL',ipFAO2,ndens2)
      Call GetMem('FockMO1','ALLO','REAL',ipFMO1,ndens2)
      Call GetMem('FockMO2','ALLO','REAL',ipFMO2,ndens2)
      Call GetMem('dens1','ALLO','REAL',ipDP,ndens2)
      Call GetMem('Dens2','ALLO','REAL',ipDP2,ndens2)
      Call GetMem('Dens4','ALLO','REAL',ipDCAS2,ndens2)
      Call GetMem('Kappa1','ALLO','REAL',ipK1,ndens2)
      Call GetMem('Kappa2','ALLO','REAL',ipK2,ndens2)
      Call GetMem('Dens5','ALLO','REAL',ipD,ndens2)
*
      !! Density matrices generated in caspt2 module
      !! They are rotated according to the (X)MS rotation matrix
      !! They are just state-specific (unrelaxed) for SS-CASPT2
C     Call Getmem('ONED','ALLO','REAL',ipG1q,ng1)
C     Call Getmem('TWOD','ALLO','REAL',ipG2q,ng2)
C     Call Getmem('ONED','ALLO','REAL',ipG1r,ntash**2)
C     Call Getmem('TWOD','ALLO','REAL',ipG2r,itri(ntash**2,ntash**2))
*
*---  Read in necessary densities.
*
      write (*,*) "calling get_d1ao"
      !! runfile_util/get_d1ao.f
      write (*,*) "ipdcas,ndens,ndens2 = ", ipdcas,ndens,ndens2,ndenslt
      Call Get_D1ao(ipDCAS,nDensLT) ! nDens2)
        write (*,*) "ipt1"
        do i = 1, ndensLT ! ndens2
          write (*,'(i4,f20.10)') i,work(ipdcas+i-1)
        end do
        write (*,*) "a"
C     irc=-1
C     iopt=0
C     Call RdRlx(irc,iopt,'D1PT22',Work(ipDP))
C     If (irc.ne.0) Goto 100
C     irc=-1
C     iopt=0
C     Call RdRlx(irc,iopt,'OVLP',rovlp)
C     If (irc.ne.0) Goto 100
*
*--- Squared density
*
        write (*,*) "b"
      Call UnFold_MCLR(Work(ipDP),Work(ipDP2))
      Call UnFold_MCLR(Work(ipDCAS),Work(ipDCAS2))
        write (*,*) "ipdcas2: ndens=",ndens
        do i = 1, ndens
          write (*,'(i4,f20.10)') i,work(ipDCAS2+i-1)
        end do
*
*==============================================================================*
*
*
*---  Make Fock matrixes
*
      ExFac=1.0D0
*
*  1) P(PT2)
*
      nFlt=0
      nBMX=0
      nBasT=0
      Do iSym = 1, nSym
         nFlt=nFlt+nBas(iSym)*(nBas(iSym)+1)/2
         nBMX=Max(nBMX,nBas(iSym))
         nBasT = nBasT + nBas(iSym)
      End Do
          write (*,*) "transformed to AO"
          call prtril(work(ipDCAS),nbast)
*
      Call FZero(Work(ipFAO1),nDens2)
      Call FockTwo_Drv(nSym,nBas,nBas,nSkip,
     &                 Work(ipDP),Work(ipDP2),Work(ipFAO1),nFlt,
     &                 ExFac,nDens2,nBMX)
      Call AO2MO(Work(ipFAO1),Work(ipFMO1))
C     write (*,*) "ipFMO1"
C     call sqprt(work(ipFMO1),nbast)
*
*  2) P(CAS)
*
      write (*,*) "ndens2 = ", ndens2
      write (*,*) "nbast  = ", nbast
      write (*,*) "ipdcas2"
      call sqprt(work(ipdcas2),nbast)
      Call FZero(Work(ipFAO2),nDens2)
      !! G_{\mu \nu} (D)
      !!   = \sum_{\rho \sigma} (\mu \nu || \rho \sigma) D_{\rho \sigma}
      !! The input of D in AO is
      !!   - Work(ipDCAS)  (triangular, but the off-diagonal is doubled)
      !!   - Work(ipDCAS2) (square)
      !! The output of G(D) is in the AO basis (triangular)
      Call FockTwo_Drv(nSym,nBas,nBas,nSkip,
     &                 Work(ipDCAS),Work(ipDCAS2),Work(ipFAO2),nFlt,
     &                 ExFac,nDens2,nBMX)
      !! ipFAO2: in triangular
      write (*,*) "ipFAO2"
      call prtril(work(ipFAO2),nbast)
      Call AO2MO(Work(ipFAO2),Work(ipFMO2))
      CALL DSCAL_(nDens2,2.0D+00,Work(ipFMO2),1) !! to be consistent with the GAMESS implementation
      !! ipFMO2: MO basis, square form
      write (*,*) "ipFMO2"
      call sqprt(work(ipFMO2),nbast)
      write (*,*) "finish after FockTwo_Drv"
C     call abend()
*
*-- Add one particle hamiltonian here ???
*
      !! Call DaXpY_(ndens2,-rovlp,Work(kint1),1,Work(ipFMO1),1)
      write (*,*) "finish after P(CAS)"
C     call abend()
*
*==============================================================================*
*
*---- CI Vector
*
*     <i|Sigma> = <i|F|0> - <0|F|0><i|0>+CI_a+CI_b
*
      Call CISigma(0,State_sym,State_sym,ipFMO1,0,0,ipci,iprci,'N')
      rE=ddot_(nconf1,Work(ipin(iprci)),1,Work(ipin(ipci)),1)
      Call Daxpy_(nconf1,1.0d0,Work(ipT),1,Work(ipin(iprci)),1)
      Call Daxpy_(nconf1,-rE,Work(ipin(ipCI)),1,Work(ipin(iprci)),1)
*==============================================================================*
*
*     D^CAS*P(PT2+CAS)
*
      Call DFock(Work(ipDCAS2),Work(ipFMO1),Work(ipK1))
*
*     D^(PT2+CAS)*P(CAS)
*
      Call DFock(Work(ipDTot2),Work(ipFMO2),Work(ipK2))
*
*==============================================================================*
*
*     Add all orbital terms together
*
      Call Daxpy_(nDens2,1.0d0,Work(ipK2),1,Work(ipK1),1)
      call daxpy_(ndens2,1.0d0,Work(ipK),1,Work(ipK1),1)
*
*==============================================================================*
*
*  OKIDOKI Two things are needed, first of all the symmetric fockmatrix for the
*  connection term:
*
*==============================================================================*
*
*---  Calculate efficent Fock matrix in AO basis (contravariant,folded)
*     and write to disk. Woba
*
      !! This is likely antisymmetrization of the orbital Lagrangian
      !! In fockgen.f, this is realized with DGeSub
      Do iS=1,nSym
        If (nbas(is).ne.0)
     *  Call DGEadd(-1.0d+00,Work(ipK1-1+ipMat(is,is)),nBas(is),'T',
     *              Work(ipK1-1+ipMat(is,is)),nBas(is),'N',
     *              Work(ipK2-1+ipMat(is,is)),nBas(is),
     *              nBas(is),nBas(is))
C    *  Call DGEadd(Work(ipK1-1+ipMat(is,is)),nBas(is),'N',
C    *              Work(ipK1-1+ipMat(is,is)),nBas(is),'T',
C    *              Work(ipK2-1+ipMat(is,is)),nBas(is),
C    *              nBas(is),nBas(is))
      End Do
      !! MO -> AO transformation ... why needed?
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call DGEMM_('N','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 Work(ipK2-1+ipMat(is,is)),nBas(iS),
     &                 0.0d0,Work(ipT2),nBas(iS))
           Call DGEMM_('N','T',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT2),nBas(iS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,Work(ipFAO1-1+ipMat(iS,iS)),nBas(is))
        End If
      End Do
      !! in integral_util/prepp.f
      !! maybe square (ipFAO1) -> triangular (ipFAO2) conversion
      Call Fold2(nsym,nbas,Work(ipFAO1),Work(ipFAO2))
      Call Put_Fock_Occ(Work(ipFAO2),nDens2)
*
*     And then I want the unsymmetric part to the MCLR
*
*
*---  Keep gradient terms woba!!!!!!!
*
      Do iS=1,nSym
        If (nbas(is).ne.0)
     *  Call DGESUB(Work(ipK1-1+ipMat(is,is)),nBas(is),'N',
     *              Work(ipK1-1+ipMat(is,is)),nBas(is),'T',
     *              rKappa(ipMat(is,is)),nBas(is),
     *              nBas(is),nBas(is))
      End Do
*
      !! zeroth- and first-order (CASSCF) part
      ! Call FockGen(One,Work(ipG1r),Work(ipG2r),Work(ipT),Fock,1)
*
*    OK Thats all folks!!
*
      Call GetMem('Temp1','FREE','REAL',ipT1,ndens2)
      Call GetMem('Temp2','FREE','REAL',ipT2,ndens2)
      Call GetMem('FockAO1','FREE','REAL',ipFAO1,ndens2)
      Call GetMem('FockAO2','FREE','REAL',ipFAO2,ndens2)
      Call GetMem('FockMO1','FREE','REAL',ipFMO1,ndens2)
      Call GetMem('FockMO2','FREE','REAL',ipFMO2,ndens2)
      Call GetMem('dens1','FREE','REAL',ipDTot,ndens2)
      Call GetMem('Dens2','FREE','REAL',ipDTOT2,ndens2)
      Call GetMem('dens3','FREE','REAL',ipDCAS,ndens2)
      Call GetMem('Dens4','FREE','REAL',ipDCAS2,ndens2)
      Call GetMem('Kappa1','FREE','REAL',ipK1,ndens2)
      Call GetMem('Kappa2','FREE','REAL',ipK2,ndens2)
      Call GetMem('Kappa','FREE','REAL',ipK,ndens2)
      Call GetMem('Dens5','FREE','REAL',ipD,ndens2)
      Call Getmem('TEMPCI','FREE','REAL',ipT,nconf1)
*
C     Call Getmem('ONED','FREE','REAL',ipG1r,ntash**2)
C     Call Getmem('TWOD','FREE','REAL',ipG2r,itri(ntash**2,ntash**2))
C     Call Getmem('TEMP','FREE','REAL',ipT,ndens2)
C     Call Getmem('TEMP','FREE','REAL',ipF,ndens2)
*
*............
      return
 100  Call SysHalt('rhs_pt2')
      end
      Subroutine DFock(DAO,FockMO,Fock)
      Implicit Real*8 (a-h,o-z)
#include "Pointers.fh"

#include "Input.fh"
#include "WrkSpc.fh"
      Real*8 DAO(*),FockMO(*),Fock(*)
      Call Getmem('Temp1','ALLO','REAL',ipT1,ndens2)
      Call Getmem('Temp2','ALLO','REAL',ipT2,ndens2)
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call DGEMM_('T','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 DAO(ipCM(is)),nBas(iS),
     &                 0.0d0,Work(ipT2),nBas(iS))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT2),nBas(iS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,Work(ipT1),nBas(is))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT1),nBas(iS),
     &                 FockMO(ipCM(iS)),nBas(iS),
     &                 0.0d0,Fock(ipCM(is)),nBas(is))
        End If
      End Do
      Call Getmem('Temp1','FREE','REAL',ipT1,ndens2)
      Call Getmem('Temp2','FREE','REAL',ipT2,ndens2)
      Return
      End
      Subroutine AO2MO(FAO ,FMO)
      Implicit Real*8 (a-h,o-z)
      Real*8 FAO(*),FMO(*)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Call GetMem('Temp','ALLO','REAL',ipT1,ndens2)
      Call GetMem('Temp','ALLO','REAL',ipT2,ndens2)
      ip=1
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call Square(FAO(ip),
     *                   Work(ipT1),
     *                   1,nBas(is),nBas(is))
           Call DGEMM_('T','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 Work(ipT1),nBas(iS),
     &                 0.0d0,Work(ipT2),nBas(iS))
           Call DGEMM_('N','N',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT2),nBas(iS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,FMO(ipMat(iS,iS)),nBas(is))
           ip=ip+nBas(is)*(nBas(iS)+1)/2
        End If
      End Do
      Call GetMem('Temp','FREE','REAL',ipT1,ndens2)
      Call GetMem('Temp','FREE','REAL',ipT2,ndens2)
      Return
      End
C
C-----------------------------------------------------------------------
C
      Subroutine MO2AO(FMO ,FAO)
      Implicit Real*8 (a-h,o-z)
      Real*8 FMO(*),FAO(*)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Call GetMem('Temp','ALLO','REAL',ipT1,ndens2)
      Call GetMem('Temp','ALLO','REAL',ipT2,ndens2)
      ip=1
      write (*,*) "Work(ipCMO)"
      call sqprt(work(ipCMO),nBas(is))
      Do iS=1,nSym
        If (nBas(is).ne.0) Then
           Call Square(FMO(ip),
     *                   Work(ipT1),
     *                   1,nBas(is),nBas(is))
      write (*,*) "Work(ipT1)"
      call sqprt(work(ipT1),nBas(is))
           Call DGEMM_('N','N',
     &                 nBas(iS),nBas(iS),nBas(iS),
     &                 1.0d0,Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 Work(ipT1),nBas(iS),
     &                 0.0d0,Work(ipT2),nBas(iS))
           Call DGEMM_('N','T',
     &                 nBas(is),nBas(iS),nBAs(iS),
     &                 1.0d0,Work(ipT2),nBas(iS),
     &                 Work(ipCMO+ipCM(iS)-1),nBas(iS),
     &                 0.0d0,FAO(ipMat(iS,iS)),nBas(is))
           ip=ip+nBas(is)*(nBas(iS)+1)/2
        End If
      End Do
      Call GetMem('Temp','FREE','REAL',ipT1,ndens2)
      Call GetMem('Temp','FREE','REAL',ipT2,ndens2)
      Return
      End
C
C-----------------------------------------------------------------------
C
      subroutine unfold_mclr(DTR,DSQ)
C
      implicit real*8 (a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
C
      dimension DTR(*),DSQ(*)
C
      !! DTR: triangular density matrix (in the AO basis)
      !! DSQ: square     density matrix (in the AO basis)
      CALL DCOPY_(nDens2,0.0D+00,0,DSQ,1)
      indT  = 0 ! index for the triangular matrix
      indS0 = 0 ! index for the square     matrix (not complete)
      Do iSym = 1, nSym
        iBas = nBas(iSym)
        Do iAO1 = 1, iBas
          Do iAO2 = 1, iAO1-1
            indT = indT + 1
            Val = DTR(indT)*0.5d+00
            DSQ(indS0+(iAO1-1)*iBas+iAO2) = Val
            DSQ(indS0+(iAO2-1)*iBas+iAO1) = Val
          End Do
          indT = indT + 1
          Val = DTR(indT)
          DSQ(indS0+(iAO1-1)*iBas+iAO2) = Val
        End Do
        indS0 = indS0 + iBas*iBas
      End Do
C
      return
C
      end subroutine unfold_mclr
C
C-----------------------------------------------------------------------
C
      subroutine fold_mclr(DSQ,DTR)
C
      implicit real*8 (a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
C
      dimension DSQ(*),DTR(*)
C
      !! DTR: triangular density matrix (in the AO basis)
      !! DSQ: square     density matrix (in the AO basis)
      CALL DCOPY_(nDensLT,0.0D+00,0,DTR,1)
      indT  = 0 ! index for the triangular matrix
      indS0 = 0 ! index for the square     matrix (not complete)
      Do iSym = 1, nSym
        iBas = nBas(iSym)
        Do iAO1 = 1, iBas
          Do iAO2 = 1, iAO1-1
            indT = indT + 1
            Val = DSQ(indS0+(iAO1-1)*iBas+iAO2)
            DTR(indT) = Val*2.0D+00
          End Do
          indT = indT + 1
          Val = DSQ(indS0+(iAO1-1)*iBas+iAO2)
          DTR(indT) = Val
        End Do
        indS0 = indS0 + iBas*iBas
      End Do
C
      return
C
      end subroutine fold_mclr
