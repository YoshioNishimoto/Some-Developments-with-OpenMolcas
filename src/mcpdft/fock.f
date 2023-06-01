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
      SUBROUTINE FOCK_m(F,FI,FP,D,P,FINT)
!
!     RASSCF program version IBM-3090: SX section
!
!     Calculation of the MCSCF fock matrix F(eq.(7) in I.J.Q.C.S14,175)
!     FP is the matrix FI+FA (FP is FA at entrance)
!     F is stored as a symmetry blocked square matrix, by columns.
!     Note that F contains all elements, also the zero elements
!     occurring when the first index is secondary.
!     F is used to construct the Brillouin elements and the first row
!     of the super-CI Hamiltonian, while FP is used as the effective
!     one-electron operator in the construction of the super-CI
!     interaction matrix.
!
!     Does not compute the Brillouin elements any more since they
!     are not needed!
!
!          ********** IBM-3090 MOLCAS Release: 90 02 22 **********
!
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: wp
      use mcpdft_output, only: debug, lf, iPrLoc
      implicit none

      real(kind=wp), dimension(*), intent(in) :: FI, D, P, FINT
      real(kind=wp), dimension(*), intent(inout) :: FP
      real(kind=wp), dimension(*), intent(out) :: F
      integer ISTSQ(8), ISTAV(8)
      real(kind=wp) ECAS0

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

      integer :: ipFMCSCF, IPRLEV
      integer :: ISTD, ISTFCK, ISTFP, ISTP, ISTZ, iSym
      integer :: ix1, jstf, N1, n2, nao, ni, nio
      integer :: nm, no, no2, nor, np, nt, ntm
      integer :: ntt, ntv, nuvx, nv, nvi, nvm
      real(kind=wp), dimension(:), allocatable :: Q
      real(kind=wp) :: casdft_en, qntm


      IPRLEV=IPRLOC(4)

      ISTSQ(1)=0
      ISTAV(1)=0
      DO iSym=2,nSym
         ISTSQ(iSym) = ISTSQ(iSym-1) + nBas(iSym-1)**2
         ISTAV(iSym) = ISTAV(iSym-1) + nBas(iSym-1)*nAsh(iSym-1)
      End Do

c     add FI to FA to obtain FP
      CALL DAXPY_(NTOT3,1.0D0,FI,1,FP,1)
C     LOOP OVER ALL SYMMETRY BLOCKS
C
      ISTFCK=0
      ISTFP=0
      ISTD=0
      IX1=0
      ISTZ=0
C
* A long loop over symmetry
      DO ISYM=1,NSYM
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NO=NORB(ISYM)
       NO2=(NO**2+NO)/2
       IF(NO == 0) GO TO 90
       CALL FZERO(F(ISTFCK+1),NO**2)
c
c      first index in F is inactive
c
       IF(NIO /= 0) THEN
        DO NP=1,NO
         DO NI=1,NIO
          N1=MAX(NP,NI)
          N2=MIN(NP,NI)
          F(ISTFCK+NO*(NP-1)+NI)=2*FP(ISTFP+(N1**2-N1)/2+N2)
         END DO
        END DO
       END IF
c
c      first index in F active
c
       IF(NAO /= 0) THEN
         ISTP=ISTORP(ISYM)+1
         JSTF=ISTORD(ISYM)+1
         NUVX=(ISTORP(ISYM+1)-ISTORP(ISYM))/NAO
c
c           first compute the Q-matrix (equation (19))
c              Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
c           P is packed in xy and pre-multiplied by 2
c                             and reordered
c
         call mma_allocate(q, no*nao, label="q")
         CALL DGEMM_('N','N',
     &              NO,NAO,NUVX,
     &              1.0d0,FINT(JSTF),NO,
     &              P(ISTP),NUVX,
     &              0.0d0,Q,NO)
c
c        active-active interaction term in the RASSCF energy
c
         ECAS0=ECAS
         DO NT=1,NAO
          NTT=(NT-1)*NO+NIO+NT
          ECAS=ECAS+0.5D0*Q(NTT)
         END DO
         IF(IPRLEV >= DEBUG) THEN
           write(lf,*) 'Two-electron contribution (Q term):', ECAS-ECAS0
         END IF
c
c        Fock matrix
c
         NTM=0
         DO NT=1,NAO
           DO NM=1,NO
             NTM=NTM+1
             QNTM=Q(NTM)
             DO NV=1,NAO
               NVI=NV+NIO
               NTV=ITRI(MAX(NT,NV))+MIN(NT,NV)+ISTD
               NVM=ITRI(MAX(NVI,NM))+MIN(NVI,NM)+ISTFP
               QNTM=QNTM+D(NTV)*FI(NVM)
             END DO
             F(ISTFCK+NO*(NM-1)+NT+NIO)=QNTM
           END DO
         END DO
         call mma_deallocate(q)
       ENDIF

90     CONTINUE
       ISTFCK=ISTFCK+NO**2
       ISTFP=ISTFP+NO2
       ISTD=ISTD+(NAO**2+NAO)/2
       IX1=IX1+NBAS(ISYM)
       ISTZ=ISTZ+(NAO**2-NAO)/2
* End of long loop over symmetry
      END DO
c
      If (iPrLev.ge.DEBUG ) then
        CASDFT_En=0.0d0
        If(KSDFT(1:3).ne.'SCF'.and.KSDFT(1:3).ne.'PAM')
     &   Call Get_dScalar('CASDFT energy',CASDFT_En)
        Write(LF,'(A,2F22.16)') ' RASSCF energy: ',
     &                  ECAS+CASDFT_En,VIA_DFT
      End If
      If(iPrLev.ge.DEBUG ) then
        Write(LF,'(A)')' MCSCF Fock-matrix in MO-basis'
        ipFMCSCF=1
        Do iSym=1,nSym
           nOr=nOrb(iSym)
           Call RecPrt(' ',' ',F(ipFMCSCF),nOr,nOr)
           ipFMCSCF=ipFMCSCF+nOr*nOr
        End Do
      end if

      If ( IPRLEV >= DEBUG ) then
         Write(LF,*)
         Write(LF,*) ' >>> Exit Fock <<< '
         Write(LF,*)
      End If
      END

      SUBROUTINE FOCK_update(F,FI,FP,D,P,Q,FINT,CMO)
!This subroutine is supposed to add the dft portions of the mcpdft fock
!matrix to the Fock matrix pieces that have already been built for the
!CASSCF portion.

!
!     RASSCF program version IBM-3090: SX section
!
!     Calculation of the MCSCF fock matrix F(eq.(7) in I.J.Q.C.S14,175)
!     FP is the matrix FI+FA (FP is FA at entrance)
!     F is stored as a symmetry blocked square matrix, by columns.
!     Note that F contains all elements, also the zero elements
!     occurring when the first index is secondary.
!     F is used to construct the Brillouin elements and the first row
!     of the super-CI Hamiltonian, while FP is used as the effective
!     one-electron operator in the construction of the super-CI
!     interaction matrix.
!
!     No longer computed the Brillouin elements!
!
!          ********** IBM-3090 MOLCASs Release: 90 02 22 **********
!
      use mspdft, only: dogradmspd, iFxyMS, iIntS
      use mcpdft_output, only: debug, lf, iPrLoc

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FI(*),FP(*),D(*),P(*),Q(*),FINT(*),F(*),CMO(*)
      integer ISTSQ(8),ISTAV(8),iTF

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
      Character*16 ROUTINE
      Parameter (ROUTINE='FOCK    ')
#include "WrkSpc.fh"

C
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF

      Call GetMem('fockt','ALLO','REAL',iTF,NTOT4)
      Call dcopy_(ntot4,[0d0],0,Work(iTF),1)


      ISTSQ(1)=0
      ISTAV(1)=0
      DO iSym=2,nSym
         ISTSQ(iSym) = ISTSQ(iSym-1) + nBas(iSym-1)**2
         ISTAV(iSym) = ISTAV(iSym-1) + nBas(iSym-1)*nAsh(iSym-1)
      End Do
C *****************************************

!      Call GetMem('ONTOPT','ALLO','Real',iTEOTP,NFINT)
!      Call GetMem('ONTOPO','ALLO','Real',iOEOTP,NTOT1)
      !Read in the one- and two- electron potentials.
!      Call Get_dArray('ONTOPT',work(iTEOTP),NFINT)
!      Call Get_dArray('ONTOPO',work(iOEOTP),NTOT1)


!I think the best way forward is to construct FI, FA (MO basis) using
!the potentials (v_pqrs and V_pq) instead of the integrals.  I think we
!want to use the full 2-body density matrix and ExFac = 1.  If we have
!FI, FA, and Q, then we can use the prescription from fock.f to
!construct the Focc term (the part of the fock matrix that we want).



!******************************************************************
!
! Build the FA terms using the potentials v.
!
!******************************************************************

!      ExFac_tmp = 1.0d0
!      Call Upd_FA_m(Work(iTEOTP),FP,D,ExFac_tmp)
!Check - does this regenerate FA if the regular integrals are passed?
!FP should contain the Fock matrix contribution that we want.


!******************************************************************
!
! Build the FI terms using the potentials v and V.
!
!******************************************************************


!      iOff1 = 0
!      Do ISYM=1,NSYM
!        Do iOrb=1,norb(iSym)
!          do jOrb=1,iOrb
!should we follow the guide of ftwo.f?
!for starters, I don't seem to have all the necessary 2-body potentials,
!right?

!          end do
!        end do
!      end do


c     add FI to FA to obtain FP
      CALL DAXPY_(NTOT3,1.0D0,FI,1,FP,1)
C     LOOP OVER ALL SYMMETRY BLOCKS

      ISTFCK=0
      ISTFP=0
      ISTD=0
      IX1=0
      ISTZ=0
C
* A long loop over symmetry
      DO ISYM=1,NSYM
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NO=NORB(ISYM)
       NO2=(NO**2+NO)/2
       N1=0
       N2=0
       IF(NO.EQ.0) GO TO 90
       CALL FZERO(Work(iTF-1+ISTFCK+1),NO**2)

!    First index in F is inactive

       IF(NIO.NE.0) THEN
        DO NP=1,NO
         DO NI=1,NIO
          N1=MAX(NP,NI)
          N2=MIN(NP,NI)
          Work(iTF-1+ISTFCK+NO*(NP-1)+NI)=2*FP(ISTFP+(N1**2-N1)/2+N2)
         END DO
        END DO
       ENDIF
c
c      first index in F active
c
       IF(NAO.NE.0) THEN

        ISTP=ISTORP(ISYM)+1
        JSTF=ISTORD(ISYM)+1
        NUVX=(ISTORP(ISYM+1)-ISTORP(ISYM))/NAO
c
c          first compute the Q-matrix (equation (19))
c
c          Q(m,v) = sum_wxy  (m|wxy) * P(wxy,v)
c
c          P is packed in xy and pre-multiplied by 2
c                            and reordered

        CALL DGEMM_('N','N',
     &              NO,NAO,NUVX,
     &              1.0d0,FINT(JSTF),NO,
     &              P(ISTP),NUVX,
     &              0.0d0,Q,NO)


!Now Q should contain the additional 2-electron part of the fock matrix
!for mcpdft, for the active region, at least.

!We should also have contributions from terms like FI and FA, too.
!FA takes care of the 1-RDM/2e- integral terms?
!FI takes care of the one-body hamiltonian and the occ/occ and occ/act
!contributions.

        E2eP=0d0
        DO NT=1,NAO
         NTT=(NT-1)*NO+NIO+NT
         E2eP=E2eP+0.5D0*Q(NTT)
        END DO
c
c       Fock matrix
c
        NTM=0
        DO NT=1,NAO
         DO NM=1,NO
          NTM=NTM+1
          QNTM=Q(NTM)
          DO NV=1,NAO
           NVI=NV+NIO
           NTV=ITRI(MAX(NT,NV))+MIN(NT,NV)+ISTD
           NVM=ITRI(MAX(NVI,NM))+MIN(NVI,NM)+ISTFP
           QNTM=QNTM+D(NTV)*FI(NVM)
          END DO
          Work(iTF-1+ISTFCK+NO*(NM-1)+NT+NIO)=QNTM
         END DO
        END DO
       ENDIF

CGLM        call recprt('Q-mat',' ',Q(1),NO,NAO)

c
c       active-active interaction term in the RASSCF energy
c
c
* End of long loop over symmetry
90     CONTINUE
       ISTFCK=ISTFCK+NO**2
       ISTFP=ISTFP+NO2
       ISTD=ISTD+(NAO**2+NAO)/2
       IX1=IX1+NBAS(ISYM)
       ISTZ=ISTZ+(NAO**2-NAO)/2
      END DO


!Now, add all components to the original Fock matrix
!      Call DAXPY_()
!      Call DAXPY_()
c
C
c     Calculate Fock matrix for occupied orbitals.
C

      If ( iPrLev.ge.DEBUG ) then
      write(6,*) 'old fock terms:'
      do i=1,Ntot4
        write(6,*) F(i)
      end do
      write(6,*) 'new fock terms to add:'
      do i=1,Ntot4
        write(6,*) Work(itF-1+i)
      end do
      end if
      Call DAXPY_(NTOT4,1.0d0,Work(iTF),1,F,1)
!I am going to add the Fock matrix temporarily to the Runfile.  I don't
!want to construct it again in MCLR in the case of gradients.
      If (iPrLev >= DEBUG ) then
        Write(LF,'(A)')'MC-PDFT Generalized Fock-matrix in MO-basis'
        ipFMCSCF=1
        Do iSym=1,nSym
           nOr=nOrb(iSym)
           Call RecPrt(' ',' ',F(ipFMCSCF),nOr,nOr)
           ipFMCSCF=ipFMCSCF+nOr*nOr
        End Do
      End If

!For MCLR
      IF(DoGradMSPD) THEN
       CALL DCopy_(nTot4,F,1,WORK(iFxyMS+(iIntS-1)*nTot4),1)
      ELSE
       Call put_dArray('Fock_PDFT',F,ntot4)
      END IF

      CALL FOCKOC_m(Q,F,CMO)
C
      Call GetMem('fockt','Free','REAL',iTF,NTOT4)
C
      RETURN
      END
