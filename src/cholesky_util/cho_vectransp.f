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
      SUBROUTINE Cho_VecTransp(Vec,Jin,Jfi,iSym,iRed,iPass)
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: MyRank, nProcs
      use ChoSwp, only: InfVec_G, IndRed
      use ChoArr, only: iL2G
#endif
      Implicit Real*8 (a-h,o-z)
      Real*8   Vec(*)
      Integer  Jin, Jfi, iSym, iRed, iPass

      Character*13 SecNam
      Parameter (SecNam = 'Cho_VecTransp')

#if defined (_MOLCAS_MPP_)
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choglob.fh"
#include "stdalloc.fh"
#include "mafdecls.fh"

#ifndef _GA_
#include "WrkSpc.fh"
#endif

      Logical LocDbg
#if defined (_DEBUGPRINT_)
      Parameter (LocDbg = .true.)
#else
      Parameter (LocDbg = .false.)
#endif


      External  ga_create_irreg, ga_destroy
      Logical   ga_create_irreg, ga_destroy, ok
      Integer   g_a
CVVP:2014 DGA is here
#ifndef _GA_
      Logical   ga_create_local
      Integer   ga_local_woff,nelm,iGAL
      External  ga_local_woff,ga_create_local
#endif
      Integer, Allocatable:: Map(:), iAdrLG(:,:), iVecR(:),
     &                       nRSL(:), MapRS2RS(:)
      Real*8, Allocatable:: VecR(:,:)
***************************************************************

      If (.not.Cho_Real_Par) Then
         If (LocDbg) Then
            Write(6,'(A,A,A)') 'Illegal call to ',SecNam,':'
            Write(6,*)
     &      'Should only be called in parallel, but Cho_Real_Par = ',
     &      Cho_Real_Par
         End If
         Call Cho_Quit('Illegal call to '//SecNam,103)
      End If

      If (iRed .eq. 2) Then
         jRed=3
      Else If (iRed .eq. 3) Then
         jRed=2
      Else
         Call Cho_Quit('iRed must be 2 or 3 in '//SecNam,104)
      End If
      nRS_l = nnBstR(iSym,iRed)   ! local  red set dimension
      nRS_g = nnBstR_G(iSym,iRed) ! global red set dimension
      nV = Jfi - Jin + 1
      nVR = 0

      call mma_allocate(iVecR,nV,Label='iVecR')
      Call cho_p_distrib_vec(Jin,Jfi,iVecR,nVR)
      Call mma_allocate(VecR,nRS_g,nVR+1,Label='VecR')

      Call mma_allocate(nRSL,nProcs,Label='nRSL')
      nRSL(:)=0
      nRSL(1+MyRank) = nRS_l  ! MyRank starts from 0
      Call Cho_GAIGOP(nRSL,nProcs,'+')

      MxRSL=nRSL(1)
      Do i=2,nProcs
         MxRSL=max(MxRSL,nRSL(i))
      End Do
      Call mma_allocate(iAdrLG,MxRSL,nProcs,Label='iAdrLG')

      Call mma_allocate(Map,nProcs,Label='Map')
      nProcs_eff = 0
      iStart = 1
      myStart = 0
      Do i=1,nProcs
         If (nRSL(i) .gt. 0) Then
            nProcs_eff = nProcs_eff + 1
            Map(nProcs_eff) = iStart
            If ((i-1) .eq. myRank) myStart = iStart
            iStart = iStart + nRSL(i)
         End If
      End Do

      If (LocDbg) Then
         Write(LuPri,*)
         Write(LuPri,*) SecNam,': debug info.'
         Write(LuPri,*) '#nodes: ',nProcs,'  myRank: ',myRank
         Write(LuPri,*) '#contributing nodes: ',nProcs_eff
         Write(LuPri,*) 'Symmetry block: ',iSym
         Write(LuPri,*) 'On this node:'
         Write(LuPri,*) 'Vector dimension : ',nRS_l
         Write(LuPri,*) 'Number of vectors: ',nV,' (',Jin,'-',Jfi,')'
         Write(LuPri,*) 'Global vector dimension : ',nRS_g
         Write(LuPri,*) 'Number of global vectors: ',nVR
         Write(Lupri,*) 'MAP:'
         Write(LuPri,*) (map(i),i=1,nProcs_eff)
      End If
CVVP:2014 Local rather than Global
#ifdef _GA_
      ok = ga_create_irreg(mt_dbl,nRS_g,nV,'Ga_Vec',Map,
     &                     nProcs_eff,1,1,g_a)
#else
      ok = ga_create_local(mt_dbl,nRS_g,nV,'Ga_Vec',g_a)
#endif
      If (.not. ok) Call Cho_Quit(SecNam//': ga_create_irreg error',101)

      If (nRS_l .gt. 0) Then
         myEnd = myStart + nRS_l - 1
#ifdef _GA_
         Call ga_put(g_a,myStart,myEnd,1,nV,Vec,nRS_l)
#else
CVVP:2014 the minimal latency and scalable putC call
         Call ga_putc(g_a,myStart,myEnd,1,nV,Vec,nRS_l)
#endif
      End If
#ifndef _GA_
      nelm=nRS_g*nV
      iGAL=ga_local_woff(g_a)
      Call Cho_GAdGOP(Work(iGAL),nelm,'+')
#else
      Call Cho_GASync()
#endif
      Jin0 = Jin - 1
      Do i=1,nVR
         jv=iVecR(i) - Jin0
#ifdef _GA_
         Call ga_get(g_a,1,nRS_g,jv,jv,VecR(:,i),nRS_g)
#else
CVVP:2014 the minimal latency and scalable getC call
         Call ga_getc(g_a,1,nRS_g,jv,jv,VecR(:,i),nRS_g)
#endif
      End Do

      ok = ga_destroy(g_a)
      If (.not. ok) Call Cho_Quit(SecNam//': ga_destroy error',101)

C --- write the reordered vec on disk

      Call Cho_P_IndxSwp()
      irc=-1
      Call Cho_X_RSCopy(irc,1,jRed)
      If (irc .ne. 0) Then
         Call Cho_Quit(SecNam//
     &                 ': Non-zero return code from Cho_X_RSCopy',
     &                 104)
      End If

      Call mma_allocate(MapRS2RS,nnBstR(iSym,1),Label='MapRS2RS')
      Call Cho_RS2RS(mapRS2RS,SIZE(mapRS2RS),jRed,iRed,iPass,iSym)
      Call Cho_P_IndxSwp()

      iAdrLG(:,:)=0
      Do i = 1,nRS_l
         i1 = IndRed(iiBstR(iSym,iRed)+i,iRed) ! addr in local rs1
         j1 = iL2G(i1) ! addr in global rs1
         j = mapRS2RS(j1-iiBstR_G(iSym,1)) ! addr in glob. rs
         iAdrLG(i,myRank+1) = j
      End Do
      Call Cho_GAIGOP(iAdrLG,SIZE(iAdrLG),'+')

      Call mma_deallocate(MapRS2RS)

      If (LocDbg) Then
         iCount=0
         Do iNode=1,nProcs
            Do iRSL=1,nRSL(iNode)
               iCount=iCount+1
            End Do
         End Do
         If (iCount .ne. nRS_g) Then
            Call Cho_Quit('nRSL error in '//SecNam,104)
         End If
      End If

      Do j=1,nVR
         Call dCopy_(nRS_g,VecR(:,j),1,VecR(:,nVr+1),1)
         iCount=0
         Do iNode=1,nProcs
            Do iRSL=1,nRSL(iNode)
               VecR(iAdrLG(iRSL,iNode),j)=VecR(1+iCount,nVR+1)
               iCount=iCount+1
            End Do
         End Do
      End Do

      If (CHO_ADRVEC .ne. 1) THEN ! only WA files!!
         Call Cho_Quit('CHO_ADRVEC error in '//SecNam,102)
      Else ! write to disk and update InfVec_G(*,3,iSym)
         iVec1=myNumCho(iSym)+1
         lTot=nRS_g*nVR
         If (lTot .gt. 0) Then
            iOpt=1
            iAdr=InfVec_G(iVec1,3,iSym)
            Call dDAfile(LuCho_G(iSym),iOpt,VecR,lTot,iAdr)
         End If
         Do iVec = 1,nVR
            jVec = iVec1 + iVec - 1
            If (jVec .lt. MaxVec) Then
               InfVec_G(jVec+1,3,iSym)= InfVec_G(jVec,3,iSym) + nRS_g
            End If
         End Do
         LastV = myNumCho(iSym) + nVR
         If (LastV .gt. MaxVec) Then
            Call Cho_Quit('Max. number of vectors exceeded in '//SecNam,
     &                    104)
         End If
      End If
      myNumCho(iSym) = myNumCho(iSym) + nVR

C --- deallocations

      Call mma_deallocate(Map)
      Call mma_deallocate(iAdrLG)
      Call mma_deallocate(nRSL)
      Call mma_deallocate(VecR)
      Call mma_deallocate(iVecR)
#else
      Call Cho_Quit(SecNam//
     &              ' should never be called in serial installation',
     &              103)
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Vec)
         Call Unused_integer(Jin)
         Call Unused_integer(Jfi)
         Call Unused_integer(iSym)
         Call Unused_integer(iRed)
         Call Unused_integer(iPass)
      End If
#endif

      End
