      module mod_solv

      implicit real*8 (a-h,o-z)

#include "rasdim.fh"
#include "rasscf.fh"

      real*8, allocatable :: savmod(:)

      contains

      subroutine mod_solvx(ifinal,ener_mod,CMO)

      use general_data, only : iSpin, nActEl, nSym, nTot1, nTot2,
     &    nBas, nIsh, nAsh, nFro, JOBIPH
      use Constants, only: Zero, One

      implicit real*8 (a-h,o-z)

#include "rasdim.fh"
#include "rasscf.fh"
#include "output_ras.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      real*8, intent(inout) :: ener_mod(*)
      real*8, intent(in) :: CMO(*)

      logical :: First,Dff,lRF,Do_DFT
      real*8 :: rdum(1)
!
!     This should be called from cictl.f after SCF is converged
!     Modify the total energy so that electrostatic solvation energies
!     are expressed better
!
      if (ifinal==1) then
        CALL mma_allocate(savmod,nRoots,Label='savmod')
        savmod = 0.0d+00
        if (iPCMroot==0) then
          !! nothing is done
          !! Although it is possible to differentiate the energy,
          !! this formulation is somewhat strange:
          !! \Delta F = Tr(D^I \cdot (V^N + V^e))
          !!          - Tr(D^{A,SA} \cdot V^e)/2
          !!          + Tr(D^{A,S} \cdot (V^N + V^e))
          !! If D^{A,SA} = D^{A,S}, \Delta F coincides with iPCMRoot <= -1
          !! def_solv = 3 in MCLR
        else if (iPCMroot==-1 .or. iPCMroot==-2 .or. iPCMroot==-3 ) then
          Call GetMem('DtmpI','Allo','Real',iTmp3,nTot2)
          Call GetMem('htmp','Allo','Real',iTmp5,nTot2)
          Call GetMem('gtmp','Allo','Real',iTmp6,nTot2)

          Call Get_dArray('D1ao',Work(iTmp3),nTot1)
          Call dCopy_(nTot2,[Zero],0,Work(iTmp5),1)
          Call dCopy_(nTot2,[Zero],0,Work(iTmp6),1)
*
          First=.True.
          Dff=.False.
          Do_DFT=.True.
          lRF=.True.

          !! iTmp5: V^N
          !! iTmp6: V^e
          !! iTmp3: D^SA (= D^I + D^{A,SA})
          iCharge=Int(Tot_Charge)
C         write (*,*)" ntot2 = ", ntot2
C         write (6,*) "AO density in mod_solv"
C         call sqprt(work(itmp3),nbas(1))
          Call DrvXV(Work(iTmp5),Work(iTmp6),Work(iTmp3),
     &               PotNuc,nTot1,First,Dff,NonEq,lRF,
     &               KSDFT,ExFac,iCharge,iSpin,rdum,rdum,
     &               1,DFTFOCK,Do_DFT)
C         write (*,*) "nuclear contributions"
C         call sqprt(work(itmp5),nbas(1))
C         write (*,*) "electron contributions"
C         call sqprt(work(itmp6),nbas(1))

          !! The integral contribution for iPCMRoot = -2 is V^N + V^e,
          !! so add V^e = V^N + V^e
          if (iPCMRoot==-2)
     *      call daxpy_(nTot1,1.0d+00,Work(iTmp5),1,Work(iTmp6),1)
          if (iPCMRoot==-3)
     *      call dcopy_(nTot1,Work(iTmp5),1,Work(iTmp6),1)

          !! AO -> MO
          LMOP = 1
          ISTFA = 1
          ISTFP = 1
          ntAsh = 0
          do ISYM=1,NSYM
            NB = NBAS(ISYM)
            NA = NASH(ISYM)
            ntAsh = ntAsh + NA
            LMOP1 = LMOP+NB*(NISH(ISYM)+NFRO(ISYM))
            if (NA /= 0) then
              call SQUARE(Work(iTmp6),Work(iTmp5),1,NB,NB)
              call DGEMM_('N','N',NB,NA,NB,
     *                    One,Work(iTmp5),NB,CMO(LMOP1),NB,
     *                    Zero,Work(iTmp3),NB)
              call DGEMM_('T','N',NA,NA,NB,
     *                       One,Work(iTmp3),NB,CMO(LMOP1),NB,
     *                       Zero,Work(iTmp6+ISTFA-1),NA)
              ISTFA = ISTFA+ITRI(NA+1)
            end if
            LMOP = LMOP+NB**2
            ISTFP = ISTFP+ITRI(NB+1)
          end do
C         write (*,*) "MO transofmred"
C         call sqprt(work(itmp6),nash(1))

C          write (*,*) "nacpar = ", nacpar
          Call GetMem('D^SA','Allo','Real',ipDSA,NACPAR)
          Call GetMem('D^SS','Allo','Real',ipDSS,NACPAR)

          call dcopy_(NACPAR,[Zero],0,Work(ipDSA),1)
          jDisk = IADR15(3)
          Do i=1,NROOTS
            Call DDafile(JOBIPH,2,Work(ipDSS),NACPAR,jDisk)
            Call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
            Call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
            Call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
            call daxpy_(NACPAR,WEIGHT(i),Work(ipDSS),1,Work(ipDSA),1)
C         write (*,*) "SS density:", i
C         call sqprt(work(ipDSS),nAsh(1))
          End Do
C         write (*,*) "SA density"
C         call sqprt(work(ipDSA),nAsh(1))

          jDisk = IADR15(3)
          do iRoots = 1, nRoots
            !! iPCMRoot = -1
            !!   \Delta F = D^SS*V(D^SA)/2
            !!   correction: Tr(D^{A,SA} - D^{A,S}) \cdot V^e)/2
            !!   def_solv = 4 in MCLR
            !! iPCMRoot = -2
            !!   \Delta F = D^SA*V(D^SA)/2 (this is not state-specific)
            !!   correction: Tr(D^{A,SA} - D^{A,S}) \cdot (V^N + V^e)
            !! iPCMRoot = -3
            !!   \Delta F = Tr((D^{A,SS} - D^{A,SA}/2) \cdot V^e)
            !!            + Tr(D^{A,SA} \cdot V^N)
            !!   correction: Tr(D^{A,SA} - D^{A,S}) \cdot (V^N + V^e)
            Call DDafile(JOBIPH,2,Work(ipDSS),NACPAR,jDisk)
            Call DDafile(JOBIPH,0,rdum,NACPAR,jDisk)
            Call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
            Call DDafile(JOBIPH,0,rdum,NACPR2,jDisk)
C         write (*,*) "SS density again:", iRoots
C         call sqprt(work(ipDSS),nAsh(1))

            call dscal_(NACPAR,-1.0d+00,Work(ipDSS),1)
            call daxpy_(NACPAR,+1.0d+00,Work(ipDSA),1,Work(ipDSS),1)
            if (iPCMRoot==-1)
     *        call dscal_(NACPAR,+0.5d+00,Work(ipDSS),1)

            call square(Work(ipDSS),Work(iTmp5),1,ntAsh,ntAsh) !?

C          write (*,*) "iroot = ", iroots
C          write (*,*) "original = ", ener_mod(iroots)
C          write (*,*) "correction = ",
C    *      ddot_(ntAsh**2,Work(itmp5),1,Work(iTmp6),1)
C           ener_mod(iRoots) = ener_mod(iRoots)
C    *        + ddot_(ntAsh**2,Work(itmp5),1,Work(iTmp6),1)
C          write (*,*) "energy for iroot = ", ener_mod(iroots)

C          write (*,*) "iroot = ", iroots
C          write (*,*) "original = ", ener(iroots,iter)
C          write (*,*) "correction = ",
C    *      ddot_(ntAsh**2,Work(itmp5),1,Work(iTmp6),1)
C           ener(iRoots,iter) = ener(iRoots,iter)
C    *        + ddot_(ntAsh**2,Work(itmp5),1,Work(iTmp6),1)
C          write (*,*) "energy for iroot = ", ener(iroots,iter)

           savmod(iroots) = ddot_(ntAsh**2,Work(itmp5),1,Work(iTmp6),1)
          end do
          Call GetMem('DtmpI','Free','Real',iTmp3,nTot2)
          Call GetMem('gtmp','Free','Real',iTmp6,nTot2)
          Call GetMem('htmp','Free','Real',iTmp5,nTot2)
          Call GetMem('D^SA','Free','Real',ipDSA,NACPAR)
          Call GetMem('D^SS','Free','Real',ipDSS,NACPAR)
        else
          write (LF,'(X,"Unknown iPCMRoot in mod_solv")')
          write (LF,'(X,"Nothing is done rather than abort")')
        end if
      else if (ifinal==2) then
        do iroots = 1, nroots
          ener(iroots,iter) = ener(iroots,iter) + savmod(iroots)
        end do
        CALL mma_deallocate(savmod)
      end if

      end subroutine mod_solvx

      end module mod_solv
