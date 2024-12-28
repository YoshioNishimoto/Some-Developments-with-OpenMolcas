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
* Copyright (C) Anders Bernhardsson                                    *
************************************************************************
      Subroutine PutRlx(D,DS,P,DAO,C)
      use spin_correlation, only: tRootGrad
      use rctfld_module, only: lRF,RslPar
      use DWSol, only: DWSol_wgt,DWZeta,W_SOLV
      Implicit Real*8 (a-h,o-z)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
#include "output_ras.fh"
#include "WrkSpc.fh"
      Character*16 ROUTINE
      Parameter (ROUTINE='PUTRLX  ')
      Real*8 D(*),DS(*),P(*),DAO(*),C(*)
      Dimension rdum(1)

      IPRLEV=IPRLOC(3)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
*
*     Read in corresponding density matrixes
*

      ! empty line before root resolved orbital gradient
      NZ=NTOT2
      kDisk = IADR15(3)  ! we need to save this address to reset later
      jDisk = IADR15(3)
      Call GetMem('TEMP','ALLO','REAL',ipDA,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipDI,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipDS,NZ)
      if (tRootGrad) then
        write(lf,*)
        do i = 1, LRoots
          call ddafile(JOBIPH, 2, D, NACPAR, jDisk)
          call ddafile(JOBIPH, 2, DS, NACPAR, jDisk)
          call ddafile(JOBIPH, 2, P, NACPR2, jDisk)
          call ddafile(JOBIPH, 0, rdum, NACPR2, jDisk)
*
*         Construc D-ACTIVE AND D-INACTIVE IN AO BASIS
*
*
          Call Get_D1I_RASSCF(C,Work(ipDI))
*
          Call GetMem('TEMP','ALLO','REAL',ipD,NACPAR)
*
          CALL DCOPY_(NACPAR,DS,1,Work(ipD),1)
          Call DBLOCK(Work(ipD))
          CALL Get_D1A_RASSCF(C,Work(ipD),Work(ipDS))
*
          CALL DCOPY_(NACPAR,D,1,Work(ipD),1)
          Call DBLOCK(Work(ipD))
          CALL Get_D1A_RASSCF(C,Work(ipd),Work(ipDA))
*
*         Construct the Fock matrix used for the connection term.
*
          NFSIZE=MAX(NACPAR,NTOT4)
          Call GetMem('TEMP','ALLO','REAL',ipF,NFSIZE)
          Call GetMem('TEMP','ALLO','REAL',ipB,NZ)
          Call GetMem('TEMP','ALLO','REAL',ipQ,NZ)
          Call GetMem('TEMP','ALLO','REAL',ipFA,NZ)
          Call GetMem('TEMP','ALLO','REAL',ipFI,NZ)
          Call GetMem('TEMP','ALLO','REAL',ipPUVX,NFINT)
          Call GetMem('TEMP','ALLO','REAL',iptuvx,nacpr2)
          IPR=0
          IF(IPRLEV.EQ.3) IPR=1
          IF(IPRLEV.EQ.4) IPR=5
          IF(IPRLEV.EQ.5) IPR=10
          Call FZero(Work(ipPUVX),NFINT)
          Call TraCtl2(C,Work(ipPUVX),Work(ipTUVX),
     &                 Work(ipDI),Work(ipFI),
     &                 Work(ipDA),Work(ipFA),ipr,lsquare,ExFac)
          CALL SGFCIN(C,WORK(ipF),Work(ipFI),Work(ipDI),Work(ipDA),
     &                                                  Work(ipDS))
          call dcopy_(ntot4,[0.0d0],0,Work(ipF),1)
          call dcopy_(ntot4,[0.0d0],0,Work(ipB),1)
*
*         Prevent FMAT from changing Active fock matrix
*
          iTmp=newfock
          newFock=-99999
*
          Call Fmat(C,Work(ipPUVX),Work(ipD),Work(ipDA),Work(ipFI),
     &              Work(ipFA))
          newFock=itmp
          IFINAL=1
          rTmp=   CBLBM
          itmp= IBLBM
          jtmp= JBLBM
          istmp= ISYMBB

          IF(ISTORP(NSYM+1).GT.0) THEN
              CALL GETMEM('ISTRP','ALLO','REAL',ipP,ISTORP(NSYM+1))
              CALL PMAT_RASSCF(P,WORK(ipP))
          Else
              ipP = ip_Dummy
          END IF
*
          Call FOCK(Work(ipF),Work(ipB),Work(ipFI),
     &              Work(ipFA),Work(ipd),WORK(ipP),Work(ipQ),
     &              Work(ipPUVX), IFINAL,C)
          CBLBM=rtmp
          IBLBM=itmp
          JBLBM=jtmp
          ISYMBB=istmp
          IF(ISTORP(NSYM+1).GT.0)
     &        CALL GETMEM('ISTRP','FREE','REAL',ipP,ISTORP(NSYM+1))

          write(6,'(6x,a,i3,5x,f12.10)')
     &        "Norm of electronic gradient for root ", i,
     &        DNRM2_(NSXS,Work(ipB),1)
*
          Call GetMem('TEMP','FREE','REAL',ipD,NACPAR)
          Call GetMem('TEMP','FREE','REAL',iptuvx,idum)
          Call GetMem('TEMP','FREE','REAL',ippuvx,idum)
          Call GetMem('TEMP','FREE','REAL',ipFI,idum)
          Call GetMem('TEMP','FREE','REAL',ipFA,idum)
          Call GetMem('TEMP','FREE','REAL',ipQ,idum)
          Call GetMem('TEMP','FREE','REAL',ipB,idum)
          Call GetMem('TEMP','FREE','REAL',ipF,NFSIZE)
        end do
      end if

      Do i=1,iRlxRoot-1
        Call DDaFile(JOBIPH,0,rdum,NACPAR,kDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPAR,kDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
        Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
      End Do
      Call DDaFile(JOBIPH,2,D,NACPAR,kDisk)
      Call DDaFile(JOBIPH,2,DS,NACPAR,kDisk)
      Call DDaFile(JOBIPH,2,P,NACPR2,kDisk)
      Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
*
*     Construc D-ACTIVE AND D-INACTIVE IN AO BASIS
*
      Call Get_D1I_RASSCF(C,Work(ipDI))
*
      Call GetMem('TEMP','ALLO','REAL',ipD,NACPAR)
*
      CALL DCOPY_(NACPAR,DS,1,Work(ipD),1)
      Call DBLOCK(Work(ipD))
      CALL Get_D1A_RASSCF(C,Work(ipD),Work(ipDS))
*
      CALL DCOPY_(NACPAR,D,1,Work(ipD),1)
      Call DBLOCK(Work(ipD))
      CALL Get_D1A_RASSCF(C,Work(ipd),Work(ipDA))
*
*     Construct the Fock matrix used for the connection term.
*
      NFSIZE=MAX(NACPAR,NTOT4)
      Call GetMem('TEMP','ALLO','REAL',ipF,NFSIZE)
      Call GetMem('TEMP','ALLO','REAL',ipB,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipQ,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipFA,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipFI,NZ)
      Call GetMem('TEMP','ALLO','REAL',ipPUVX,NFINT)
      Call GetMem('TEMP','ALLO','REAL',iptuvx,nacpr2)
      IPR=0
      IF(IPRLEV.EQ.3) IPR=1
      IF(IPRLEV.EQ.4) IPR=5
      IF(IPRLEV.EQ.5) IPR=10
      Call FZero(Work(ipPUVX),NFINT)
      Call TraCtl2(C,Work(ipPUVX),Work(ipTUVX),
     &             Work(ipDI),Work(ipFI),
     &             Work(ipDA),Work(ipFA),ipr,lsquare,ExFac)
*
      if (IPCMROOT <=0 .or. DWZeta /= 0.0d+00) then
        !! Polarize PCM etc with weighted density
        kDisk = IADR15(3)
        Call GetMem('TEMP','ALLO','REAL',ipWRK1,NACPAR)
        Call GetMem('TEMP','ALLO','REAL',ipWRK2,NACPAR)
        Call GetMem('TEMP','ALLO','REAL',ipDA_ave,MAX(NACPAR,NZ))
        Call GetMem('TEMP','ALLO','REAL',ipDS_ave,MAX(NACPAR,NZ))
*
        call dcopy_(NACPAR,[0.0d+00],0,Work(ipDA_ave),1)
        call dcopy_(NACPAR,[0.0d+00],0,Work(ipDS_ave),1)
*
        call DWSol_wgt(ENER(:,ITER))
        Do i=1,lRoots
          wgt = W_SOLV(i)
          Call DDaFile(JOBIPH,2,Work(ipWRK1),NACPAR,kDisk)
          Call DDaFile(JOBIPH,2,Work(ipWRK2),NACPAR,kDisk)
          Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
          Call DDaFile(JOBIPH,0,rdum,NACPR2,kDisk)
          if (wgt == 0.0d+00) cycle
!          if (i.eq.2) wgt = 0.40d+00
!          if (i.eq.1) wgt = 0.40d+00
!          if (i.eq.3) wgt = 0.20d+00
          call daxpy_(NACPAR,wgt,Work(ipWRK1),1,Work(ipDA_ave),1)
          call daxpy_(NACPAR,wgt,Work(ipWRK2),1,Work(ipDS_ave),1)
        End Do
*
*       Construc D-ACTIVE AND D-INACTIVE IN AO BASIS
*
        CALL DCOPY_(NACPAR,Work(ipDS_ave),1,Work(ipD),1)
        call dcopy_(NZ,Work(ipDS),1,Work(ipDS_ave),1)
        Call DBLOCK(Work(ipD))
        CALL Get_D1A_RASSCF(C,Work(ipD),Work(ipDS))
*
        CALL DCOPY_(NACPAR,Work(ipDA_ave),1,Work(ipD),1)
        call dcopy_(NZ,Work(ipDA),1,Work(ipDA_ave),1)
        Call DBLOCK(Work(ipD))
        CALL Get_D1A_RASSCF(C,Work(ipd),Work(ipDA))
*
        Call GetMem('TEMP','FREE','REAL',ipWRK1,NACPAR)
        Call GetMem('TEMP','FREE','REAL',ipWRK2,NACPAR)
      end if
*
C     write (6,*) "sgfcin from putrlx doing SS density"
      CALL SGFCIN(C,WORK(ipF),Work(ipFI),Work(ipDI),Work(ipDA),
     &                                              Work(ipDS))
      call dcopy_(ntot4,[0.0d0],0,Work(ipF),1)
      call dcopy_(ntot4,[0.0d0],0,Work(ipB),1)
*
*     Prevent FMAT from changing Active fock matrix
*
      iTmp=newfock
      newFock=-99999
*
      if (IPCMROOT <=0 .or. DWZeta /= 0.0d+00) then
        !! The rest uses state-specific density, so restore them
        call dcopy_(NZ,Work(ipDA_ave),1,Work(ipDA),1)
        call dcopy_(NZ,Work(ipDS_ave),1,Work(ipDS),1)
        Call GetMem('TEMP','FREE','REAL',ipDA_ave,MAX(NACPAR,NZ))
        Call GetMem('TEMP','FREE','REAL',ipDS_ave,MAX(NACPAR,NZ))
*
        Call GetMem('TEMP','ALLO','REAL',ipWRK1,nTot1)
        Call GetMem('TEMP','ALLO','REAL',ipWRK2,nTot1)
        Call Fold(nSym,nBas,Work(ipD),Work(ipWRK1))
        Call Fold(nSym,nBas,Work(ipDA),Work(ipWRK2))
        Call Daxpy_(nTot1,1.0D+00,Work(ipWRK1),1,Work(ipWRK2),1)
        Call Put_dArray('D1ao',Work(ipWRK2),nTot1)
        Call Fold(nSym,nBas,Work(ipDS),Work(ipWRK1))
        Call Put_dArray('D1sao',Work(ipWRK1),nTot1)
        Call GetMem('TEMP','FREE','REAL',ipWRK1,nTot1)
        Call GetMem('TEMP','FREE','REAL',ipWRK2,nTot1)
      end if
      Call Fmat(C,Work(ipPUVX),Work(ipD),Work(ipDA),Work(ipFI),
     &          Work(ipFA))
      newFock=itmp
      IFINAL=1
      rTmp=   CBLBM
      itmp= IBLBM
      jtmp= JBLBM
      istmp= ISYMBB

      IF(ISTORP(NSYM+1).GT.0) THEN
          CALL GETMEM('ISTRP','ALLO','REAL',ipP,ISTORP(NSYM+1))
          CALL PMAT_RASSCF(P,WORK(ipP))
      Else
          ipP = ip_Dummy
      END IF
*
      Call FOCK(Work(ipF),Work(ipB),Work(ipFI),
     &          Work(ipFA),Work(ipd),WORK(ipP),Work(ipQ),Work(ipPUVX),
     &          IFINAL,C)
      CBLBM=rtmp
      IBLBM=itmp
      JBLBM=jtmp
      ISYMBB=istmp
      IF(ISTORP(NSYM+1).GT.0)
     &    CALL GETMEM('ISTRP','FREE','REAL',ipP,ISTORP(NSYM+1))

      RlxGrd=DNRM2_(NSXS,Work(ipB),1)
*
      Call GetMem('TEMP','FREE','REAL',ipD,NACPAR)
      Call GetMem('TEMP','FREE','REAL',iptuvx,idum)
      Call GetMem('TEMP','FREE','REAL',ippuvx,idum)
      Call GetMem('TEMP','FREE','REAL',ipFI,idum)
      Call GetMem('TEMP','FREE','REAL',ipFA,idum)
      Call GetMem('TEMP','FREE','REAL',ipQ,idum)
      Call GetMem('TEMP','FREE','REAL',ipB,idum)
      Call GetMem('TEMP','FREE','REAL',ipF,NFSIZE)
*
* Add up one electron densities
*
      call daxpy_(nZ,1.0d0,Work(ipDA),1,Work(ipDI),1)
      Call Fold(nSym,nBas,Work(ipDI),DAO)

      Call GetMem('TEMP','FREE','REAL',ipDS,NZ)
      Call GetMem('TEMP','FREE','REAL',ipDA,idum)
      Call GetMem('TEMP','FREE','REAL',ipDI,idum)

      return
      end
