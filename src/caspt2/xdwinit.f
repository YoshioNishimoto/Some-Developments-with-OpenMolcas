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
* Copyright (C) 2012, Per Ake Malmqvist                                *
*               2019, Stefano Battaglia                                *
************************************************************************
      subroutine xdwinit(Heff,H0,U0)

      use definitions, only: wp, iwp, u6
      use caspt2_output, only: iPrGlb, usual, verbose, debug, insane
      use caspt2_gradient, only: do_grad

      implicit none

#include "rasdim.fh"
#include "caspt2.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"

      Real(kind=wp),intent(inout) :: Heff(Nstate,Nstate)
      Real(kind=wp),intent(inout) :: H0(Nstate,Nstate)
      Real(kind=wp),intent(inout) :: U0(Nstate,Nstate)

      Real(kind=wp) :: wgt,FIJ

      Integer(kind=iwp) :: LCIRef,LCIXMS,LDAVE,iState,iDisk,I,J


* Allocate memory for CI array state averaged 1-RDM
      call getmem('LCI','ALLO','REAL',LCI,NCONF)
      call getmem('LDAVE','ALLO','REAL',LDAVE,NDREF)
      call dcopy_(NDREF,[0.0_wp],0,WORK(LDAVE),1)

* Set the weight for the density averaging
      wgt = 1.0_wp/real(Nstate,kind=wp)

* Loop over all states to compute the state-average density matrix
      do Istate=1,Nstate

        if (ISCF.NE.0) then
* Special case for a single Slater determinant
          WORK(LCI)=1.0_wp
        else
* Get the CI array
          call loadCI(WORK(LCI), Istate)
        end if

* Compute 1-particle active density matrix GAMMA1
        call POLY1(WORK(LCI))

* Restructure GAMMA1 as DREF array
        call GETDREF(WORK(LDREF))

* Average the density
        call DAXPY_(NDREF,wgt,WORK(LDREF),1,WORK(LDAVE),1)

      end do

      if (IPRGLB.GE.INSANE) then
        write(u6,*)' State-average 1-RDM'
        do I=1,NASHT
          write(u6,'(1x,14f10.6)')(WORK(LDAVE+(I*(I-1))/2+J-1),J=1,I)
        end do
        write(u6,*)
      end if

* Copy the state-average 1-RDM into LDREF and release memory
* for both the CI array and LDAVE
      call dcopy_(NDREF,WORK(LDAVE),1,WORK(LDREF),1)
      call getmem('LCI','FREE','REAL',LCI,NCONF)
      call getmem('LDAVE','FREE','REAL',LDAVE,NDREF)

* Load CASSCF MO coefficients
      call getmem('LCMO','ALLO','REAL',LCMO,NCMO)
      iDisk=IAD1M(1)
      call ddafile(LUONEM,2,WORK(LCMO),NCMO,iDisk)

* Build the state-average Fock matrix in MO basis
      if (IfChol) then
* INTCTL2 uses TraCho2 and FMatCho to get matrices in MO basis.
        call INTCTL2(.false.)
      else
* INTCTL1 uses TRAONE and FOCK_RPT2 to get the matrices in MO basis.
        call INTCTL1(WORK(LCMO))
      end if

* Loop again over all states to compute H0 in the model space
* Loop over ket functions
      do J=1,Nstate
* Loop over bra functions
        do I=1,Nstate
* Compute matrix element <I|F|J> and store it into H0
          FIJ = 0.0_wp
          call FOPAB(WORK(LFIFA),I,J,FIJ)
          H0(I,J) = FIJ
        end do
      end do
* End of loop over states

      if (IPRGLB.ge.USUAL) then
        write(u6,*)
        write(u6,*)' H0 in the original model space basis:'
        call prettyprint(H0,Nstate,Nstate)
      end if

* Diagonalize H0 in the model space
      call eigen(H0,U0,Nstate)

* Transform the Fock matrix in the new basis
      call transmat(H0,U0,Nstate)
        if (IPRGLB.ge.USUAL) then
          write(u6,*)' H0 eigenvectors:'
          call prettyprint(U0,Nstate,Nstate)
        end if
        if (IPRGLB.ge.DEBUG) then
          write(u6,*)' H0 in the rotated model space basis:'
          call prettyprint(H0,Nstate,Nstate)
        end if

* As well as Heff
      call transmat(Heff,U0,Nstate)
      if (IPRGLB.ge.VERBOSE) then
        write(u6,*)' Heff[1] in the rotated model space basis:'
        call prettyprint(Heff,Nstate,Nstate)
      end if

* Mix the CI arrays according to the H0 eigenvectors. Assume we can
* put all the original ones in memory, but put the resulting vectors
* one by one in a buffer.
        if (IPRGLB.ge.VERBOSE) then
          write(u6,'(A)')' The CASSCF states are now rotated'//
     &                  ' according to the H0 eigenvectors'
          write(u6,*)
        end if

      call getmem('CIREF','ALLO','REAL',LCIref,Nstate*Nconf)
* Load the CI arrays into memory
      do I=1,Nstate
        call loadCI(WORK(LCIref+Nconf*(I-1)),I)
      end do

      call getmem('CIXMS','ALLO','REAL',LCIXMS,Nconf)
      do J=1,Nstate
* Transform the states
        call dgemm_('N','N',Nconf,1,Nstate,
     &              1.0_wp,WORK(LCIREF),Nconf,U0(:,J),Nstate,
     &              0.0_wp,WORK(LCIXMS),Nconf)

* Write the rotated CI coefficients back into LUCIEX and REPLACE the
* original unrotated CASSCF states. Note that the original states
* are still available in the JobIph file
        call writeCI(WORK(LCIXMS),J)

        if (IPRGLB.ge.VERBOSE) then
          write(u6,'(1x,a,i3)')
     &    ' The CI coefficients of rotated model state nr. ',MSTATE(J)
          call PRWF_CP2(STSYM,NCONF,WORK(LCIXMS),CITHR)
        end if
      end do
C
      If (do_grad) call dcopy_(NCMO,WORK(LCMO),1,WORK(LCMOPT2),1)

* Release all memory
      call getmem('CIREF','FREE','REAL',LCIREF,Nstate*NCONF)
      call getmem('CIXMS','FREE','REAL',LCIXMS,NCONF)
      call getmem('LCMO','FREE','REAL',LCMO,NCMO)


      return
      end

