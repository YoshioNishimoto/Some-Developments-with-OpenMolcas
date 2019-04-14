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
* Copyright (C) 2019, Giovanni Li Manni                                *
************************************************************************
      Subroutine orbsort(iUHF)
************************************************************************
*                                                                      *
* Purpose: It sorts orbitals of INPORB2 as in INPORB1 based on overlap.*
*          Sorted orbitals are printed in file.SortOrb                 *
*          If mixing of orbitals in INPORB2 has occurred and made      *
*          sorting not obvious, a WARNING is given and sorting won't be*
*          performed.                                                  *
*                                                                      *
*          In the future sorting will be possible by manual user input.*
*                                                                      *
*          G. Li Manni, Max-Planck Institute Stuttgart, April 2019.    *
*                                                                      *
* Example of input:                                                    *
*                                                                      *
*   &GATEWAY                                                           *
*    coord                                                             *
*     file.xyz                                                         *
*    basis                                                             *
*     ano-rcc-mb                                                       *
*   &SEWARD                                                            *
*    oneonly                                                           *
*                                                                      *
*   >>COPY $CurrDir/scf_reference.ScfOrb  INPORB1                      *
*   >>COPY $CurrDir/scf_ToBeSorted.ScfOrb INPORB2                      *
*   &EXPBAS                                                            *
*    NoEx                                                              *
*    Sort                                                              *
*   END                                                                *
*                                                                      *
************************************************************************

      Implicit Real*8 (a-h,o-z)
      integer iUHF, iErr, iDummy, iRc
      integer IndT(7,8)
      integer IndType(56)
      Logical AddFragments, Found, Debug
      real Dummy
      character vTitle*80
#include "itmax.fh"
#include "info.fh"
#include "info_expbas.fh"
#include "WrkSpc.fh"

      Debug=.false.
*                                                                      *
************************************************************************
*                                                                      *
      write(6,*) 'Sorting option enabled.'
      write(6,*) 'Attempting sorting of INPORB2 as INPORB1...'

      if(iUHF.eq.1) then
        write(6,*) 'SORT keyword not implemented for UHF wf!'
        Call Abend()
      end if
*                                                                      *
************************************************************************
*                                                                      *
*     Read the characteristics of all different basis sets,            *
*     provide each atom with a nuclear charge and establish            *
*     a link between an atom and its basis set                         *
*                                                                      *
*     This call will also fill info.fh and the dynamic storage in      *
*     Work(ipInf)                                                      *
*
      AddFragments=.true.

      If (AddFragments) Then
        Call Inter1_FAIEMP(AtomLabel,iBas_Lab,Coor,Znuc,nAtom,ipInf)
      Else
        Call Inter1       (AtomLabel,iBas_Lab,Coor,Znuc,nAtom,ipInf)
      End If
      Call Qpg_iArray('nOrb',Found,nData)
      If (Found) Then
         Call Get_iArray('nOrb',nOrb,nData)
      Else
         Call iCopy(nIrrep,nBas,1,nOrb,1)
      End If
*                                                                      *
************************************************************************
*     Compute memory requirements and allocate memory                  *
*
      nB=0
      Do iSym=0,nIrrep-1
       nB=nB+nBas(iSym)
      End Do
      Call GetMem('CMO1','ALLO','REAL',ipCMO1,nB**2)
      Call GetMem('CMO2','ALLO','REAL',ipCMO2,nB**2)
      Call GetMem('OCC2','ALLO','REAL',ipOCC2,nB)
      Call GetMem('EOR2','ALLO','REAL',ipEOrb2,nB)
      Call GetMem('IND2','ALLO','INTE',ipIndt2,nB)
      Call GetMem('AOSM','ALLO','REAL',ipAOSmat,nB*(nB+1)/2)

      If(Debug) then
        write(6,*) 'DEBUG: nIrrep, nBasTot, nBas(iSym) =',nIrrep,nB,nBas
      end if
*                                                                      *
************************************************************************
*                   Read Reference CMOs from INPORB1 file              *
*                                                                      *
      Call RdVec('INPORB1',1,'C',nIrrep,nBas,nBas,
     &           Work(ipCMO1),Dummy,Dummy,iDummy,vTitle,1,iErr)
      if(iErr.ne.0) then
        write(6,*) 'Sorry: Something went wrong when reading INPORB1'
        iRc=1
        return
      end if
      write(6,*) 'INPORB1 Title:', vTitle
*                                                                      *
************************************************************************
*                   Read CMOs to be sorted from INPORB2 file           *
*                                                                      *
* In output Work(ipCMO2) is destroyed and replaced by the sorted  CMOs.*
*                                                                      *
      Call RdVec('INPORB2',1,'COEI',nIrrep,nBas,nBas,
     &           Work(ipCMO2),Work(ipOCC2),Work(ipEOrb2),
     &           iWork(ipIndt2),vTitle,1,iErr)
      if(iErr.ne.0) then
        write(6,*) 'Sorry: Something went wrong when reading INPORB2'
        iRc=1
        return
      end if
      write(6,*) 'INPORB2 Title:', vTitle
*                                                                      *
************************** Read AO Overlap Mat *************************
*   Stored and read in Lower Triangular form                           *
*                                                                      *
      iErr=-1
      iOpt=6
      iComp=1
      iSyLbl=1
*      Label='Mltpl  0'
      Call RdOne(iErr,iOpt,'Mltpl  0',iComp,Work(ipAOSmat),iSyLbl)
      if(iErr.ne.0) then
        write(6,*) 'Sorry: Something wrong in reading AO overlap matrix'
        iRc=1
        return
      end if

      If(Debug) then
         Write(6,*)
         Write(6,*) ' Overlap in AO basis in EXPBAS'
         Write(6,*) ' ---------------------'
         Write(6,*)
         iOff=0
         Do iSym = 0,nIrrep-1
          iBas = nBas(iSym)
          Call TriPrt(' ','(5G17.11)',Work(ipAOSmat+iOff),iBas)
          iOff = iOff + (iBas*iBas+iBas)/2
         End Do
      end if
*                                                                      *
***************************** Get Swapping-Matrix **********************
*  We will compute C(sorted) as:                                       *
*         C(sorted) =    C(tosort) * W                                 *
*  where the swapping matrix, W, is obtained as:                       *
*                     W = C(tosort)^T * S * C(ref)                     *
*                                                                      *
      ipOrb1  = ipCMO1
      ipOrb2  = ipCMO2
      ipOvlp  = ipAOSmat
      zero = 0.0d0
      Do iSym = 0,nIrrep-1
        nB = nBas(iSym)
        if(nB.gt.0) then
          Call GetMem('SCRA','ALLO','Real',lScra,nB*nB)
          Call GetMem('SCR1','ALLO','Real',lScr1,nB*nB)
          Call GetMem('ISCR','ALLO','INTE',ipIntScr,nB*nB)
*  Build the Overlap matrix in its square form.                        *
          Call Square(Work(ipOvlp),Work(lScra),1,nB,nB)
* Here we compute W^T, stored in Work(lScra) after two calls to DGEMM_ *
          Call dCopy_(nB*nB,zero,0,Work(lScr1),1)
          Call DGEMM_('T','N',
     &                nB,nB,nB,
     &                1.0d0,Work(lScra),nB,
     &                Work(ipOrb1),nB,
     &                0.0d0,Work(lScr1),nB)
          Call dCopy_(nB*nB,zero,0,Work(lScra),1)
          Call DGEMM_('T','N',
     &                nB,nB,nB,
     &                1.0d0,Work(ipOrb2),nB,
     &                Work(lScr1),nB,
     &                0.0d0,Work(lScra),nB)
          If(Debug) then
            Write(6,*) ' ** Print Swap-matrix Symmetry Block (real) **'
            ioff=0
            do i = 1, nB
              write(6, '(68F5.2)') (Work(lScra+ioff+j-1),j=1,nB)
             ioff = ioff + nB
            end do
          End if
* W is rounded to nearest integer, as to serve for pure swapping.      *
* It will be stored back as real in Work(lScra) as real                *
          do i = 0, nB*nB - 1
              iWork(ipIntScr+i) = nint(Work(lScra+i))
              Work(lScra+i)     = real(iWork(ipIntScr+i))
          end do
* Check for rowsums=0 or/and row-sums>1. If positive throw a WARNING   *
* as it is simptomatic of large orbital mixing in C(tosort)            *
* with respect to C(ref) CMOs.                                         *
          ioff = 0
          do i = 0, nB-1
            iSum = 0
            do j = 0, nB-1
              iSum = iSum + iWork(ipIntScr+iOff+j)
            end do
            if(isum.eq.0.or.isum.gt.1) then
              write(6,*) '*********************************************'
              write(6,*) '********** WARNING WARNING WARNING **********'
              write(6,*) '*********************************************'
              write(6,*)
     &  'Orbitals differ too much!! Relative re-ordering not obvious!'
              write(6,*) 'Manual sorting is required! Sorry!'
              write(6,*)
              Write(6,*) ' ** Swap-matrix Symmetry Block (int) **'
              ioff2=0
              do ii = 1, nB
                write(6, '(68I2)') ((iWork(ipIntScr+ioff2+j-1)),j=1,nB)
                ioff2 = ioff2 + nB
              end do
              Call Abend()
            end if
            ioff= ioff + nB
          end do
          If(Debug) then
            Write(6,*) ' ** Print Swap-matrix Symmetry Block (int) **'
            ioff=0
            do i = 1, nB
              write(6, '(68I2)') ((iWork(ipIntScr+ioff+j-1)),j=1,nB)
             ioff = ioff + nB
            end do
            Write(6,*) ' ** Print Swap-matrix Symmetry Block (real) **'
            ioff=0
            do i = 1, nB
              write(6, '(68F4.1)') (Work(lScra+ioff+j-1),j=1,nB)
             ioff = ioff + nB
            end do
          end if
          Call GetMem('ISCR','FREE','INTE',ipIntScr,nB*nB)
* Here we compute C(sorted) = C(tosort)*W... Stored in Work(lScr1)   *
          Call dCopy_(nB*nB,zero,0,Work(lScr1),1)
          Call DGEMM_('N','N',
     &                nB,nB,nB,
     &                1.0d0,Work(ipOrb2),nB,
     &                Work(lScra),nB,
     &                0.0d0,Work(lScr1),nB)
          If(Debug) then
             Write(6,*) ' Sorted CMOs :'
             call wrtmat(Work(lScr1),nB,nB,nB,nB)
          end if
* The symmetry-block results (Work(lScr1)) are copied to
* (and in so doing overwriting) the Work(ipOrb2) array.
          Call dCopy_(nB*nB,Work(lScr1),1,Work(ipOrb2),1)
          Call GetMem('SCRA','FREE','Real',lScra,nB*nB)
          Call GetMem('SCR1','FREE','Real',lScr1,nB*nB)
          ipOvlp = ipOvlp + nB*(nB+1)/2
          ipOrb1 = ipOrb1 + nB**2
          ipOrb2 = ipOrb2 + nB**2
        end if
      End do
*                                                                      *
********************** Make typeindex information **********************
*   Pretty lame typeindex loop. The type index should be read from     *
*   INPORB2 and sorted according  to the swapping matrix.              *
*   This feature is not supported by the WrVec yet.                    *
      do i = 1, 7
        do j = 1,8
          IndT(i,j) = 0
        end do
      end do

      j=ipIndt2-1
      Do iSym=1,nIrrep
        Do k=1,nBas(iSym-1)
          kIndT = iWork(j+k)
          IndT(kIndT,iSym)=IndT(kIndT,iSym)+1
*          write(6,*) 'IndT =', IndT(kIndT,iSym)
        End Do
        j=j+nBas(iSym)
      End Do

      icount = 0
      do i = 1,7
        do j = 1,8
          icount = icount +1
          IndType(icount) = IndT(j,i)
*          write(6,*) 'IndType =', IndType(icount)
        end do
      end do
*                                                                      *
****************** Write Sorted CMOs into SortOrb file *****************
*                                                                      *
* For now OccNumb, OrbEnergy, spaceIndices are copied from the INPORB2.*
* TO BE DONE: They should be swapped as for the orbitals!              *
       VTitle = 'Sorted Orbitals SORTORB'
       Call WrVec('SORTORB',1,'COEI',nIrrep,nBas,nBas,Work(ipCMO2),
     &            Work(iPOcc2),Work(ipEorb2),IndType,VTitle)
*                                                                      *
***************************** Clean-up memory **************************
*                                                                      *
      Call GetMem('CMO1','FREE','REAL',ipCMO1,nB**2)
      Call GetMem('CMO2','FREE','REAL',ipCMO2,nB**2)
      Call GetMem('OCC2','FREE','REAL',ipOCC2,nB)
      Call GetMem('EOR2','FREE','REAL',ipEOrb2,nB)
      Call GetMem('IND2','FREE','INTE',ipIndt2,nB)
      Call GetMem('AOSM','FREE','REAL',ipAOSmat,nB*(nB+1)/2)
      Call ClsSew

      Return
      End
