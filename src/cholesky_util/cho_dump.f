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
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine Cho_Dump(irc,Lunit)
C
C     T.B. Pedersen, March 2005.
C
C     Purpose: print all entries in include files
C              choorb.fh
C              cholesky.fh
C              choptr.fh
C              chosubscr.fh
C
C     On input, Lunit is the logical unit to print to...
C
      Implicit None
      Integer irc, Lunit
#include "choorb.fh"
#include "choprint.fh"
#include "chosubscr.fh"
#include "cholesky.fh"

      Character*8 SecNam
      Parameter (SecNam = 'Cho_Dump')

      Integer i, j

      irc = 0

      Write(Lunit,*)
      Write(Lunit,*)
      Write(Lunit,*) '>>> Output from ',SecNam,':'
      Write(Lunit,*)
      Call Cho_Flush(Lunit)

C     choorb.fh.
C     -----------

      Write(Lunit,*) '*** Contents of choorb.fh:'
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'iBas : ',(iBas(i),i=1,8)
      Write(Lunit,*) 'nBas : ',(nBas(i),i=1,8)
      Write(Lunit,*) 'XnBas: ',(XnBas(i),i=1,8)
      Write(Lunit,*) 'nBasT: ',nBasT
      Write(Lunit,*)
      Call Cho_Flush(Lunit)

C     cholesky.fh.
C     -------------

      Write(Lunit,*) '*** Contents of cholesky.fh:'
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'ThrDef        : ',ThrDef
      Write(Lunit,*) 'ThrCom        : ',ThrCom
      Write(Lunit,*) 'ThrDiag       : ',ThrDiag
      Write(Lunit,*) 'Tol_DiaChk    : ',Tol_DiaChk
      Write(Lunit,*) 'ThrNeg        : ',ThrNeg
      Write(Lunit,*) 'WarNeg        : ',WarNeg
      Write(Lunit,*) 'TooNeg        : ',TooNeg
      Write(Lunit,*) 'nSym          : ',nSym
      Write(Lunit,*) 'lBuf          : ',lBuf
      Write(Lunit,*) 'MinQual       : ',MinQual
      Write(Lunit,*) 'MaxQual       : ',MaxQual
      Write(Lunit,*) 'IFCSew        : ',IFCSew
      Write(Lunit,*) 'Mode_Screen   : ',Mode_Screen
      Write(Lunit,*) 'Cho_DecAlg    : ',Cho_DecAlg
      Write(Lunit,*) 'Cho_DecAlg_Def: ',Cho_DecAlg_Def
      Write(Lunit,*) 'iAlQua        : ',iAlQua
      Write(Lunit,*) 'MxShPr        : ',MxShPr
      Write(Lunit,*) 'ModRst        : ',ModRst
      Write(Lunit,*) 'Run_Mode      : ',Run_Mode
      Write(Lunit,*) 'ScDiag        : ',ScDiag
      Write(Lunit,*) 'ChkOnly       : ',ChkOnly
      Write(Lunit,*) 'Cho_IntChk    : ',Cho_IntChk
      Write(Lunit,*) 'Cho_MinChk    : ',Cho_MinChk
      Write(Lunit,*) 'Cho_UseAbs    : ',Cho_UseAbs
      Write(Lunit,*) 'Cho_TrcNeg    : ',Cho_TrcNeg
      Write(Lunit,*) 'Cho_ReOrd     : ',Cho_ReOrd
      Write(Lunit,*) 'Cho_DiaChk    : ',Cho_DiaChk
      Write(Lunit,*) 'Cho_TstScreen : ',Cho_TstScreen
      Write(Lunit,*) 'Cho_1Center   : ',Cho_1Center
      Write(Lunit,*) 'Cho_No2Center : ',Cho_No2Center
      Write(Lunit,*) 'Cho_PreScreen : ',Cho_PreScreen
      Write(Lunit,*) 'Cho_SimP      : ',Cho_SimP
      Write(Lunit,*) 'Cho_Fake_Par  : ',Cho_Fake_Par
      Write(Lunit,*) 'RstDia        : ',RstDia
      Write(Lunit,*) 'RstCho        : ',RstCho
      Write(Lunit,*) 'Did_DecDrv    : ',Did_DecDrv
      Write(Lunit,*) 'HaltIt        : ',HaltIt
      Write(Lunit,*) 'Trace_Idle    : ',Trace_Idle
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'LuCho: ',(LuCho(i),i=1,8)
      Write(Lunit,*) 'LuSel: ',(LuSel(i),i=1,8)
      Write(Lunit,*) 'LuTmp: ',(LuTmp(i),i=1,8)
      Write(Lunit,*) 'LuPri: ',LuPri
      Write(Lunit,*) 'LuScr: ',LuScr
      Write(Lunit,*) 'LuRed: ',LuRed
      Write(Lunit,*) 'LuRst: ',LuRst
      Write(Lunit,*) 'LuMap: ',LuMap
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'nShell   : ',nShell
      Write(Lunit,*) 'nnShl_Tot: ',nnShl_Tot
      Write(Lunit,*) 'nnShl    : ',nnShl
      Write(Lunit,*) 'MxORSh   : ',MxORSh
      Write(Lunit,*) 'Mx2Sh    : ',Mx2Sh
      Do j = 1,3
         Write(Lunit,*) 'iiBstR : ',(iiBstR(i,j),i=1,8)
      End Do
      Do j = 1,3
         Write(Lunit,*) 'nnBstR : ',(nnBstR(i,j),i=1,8)
      End Do
      Write(Lunit,*) 'nnBstRT: ',(nnBstRT(i),i=1,3)
      Write(Lunit,*) 'mmBstRT: ',mmBstRT
      Write(Lunit,*) 'nQual  : ',(nQual(i),i=1,8)
      Write(Lunit,*) 'iOffQ  : ',(iOffQ(i),i=1,8)
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'DiaMax: ',(DiaMax(i),i=1,8)
      Write(Lunit,*) 'DiaMin: ',(DiaMin(i),i=1,8)
      Write(Lunit,*) 'Damp  : ',(Damp(i),i=1,2)
      Write(Lunit,*) 'Span  : ',Span
      Write(Lunit,*) 'XlDiag: ',XlDiag
      Write(Lunit,*) 'DiaMnZ: ',DiaMnZ
      Write(Lunit,*) 'Thr_PreScreen: ',Thr_PreScreen
      Write(Lunit,*) 'iABMnZ: ',iABMnZ
      Write(Lunit,*) 'nnZTot: ',nnZTot
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'NumCho   : ',(NumCho(i),i=1,8)
      Write(Lunit,*) 'NumChT   : ',NumChT
      Write(Lunit,*) 'MaxVec   : ',MaxVec
      Write(Lunit,*) 'MaxRed   : ',MaxRed
      Write(Lunit,*) 'BlockSize: ',BlockSize
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'ShA     : ',ShA
      Write(Lunit,*) 'ShB     : ',ShB
      Write(Lunit,*) 'ShAB    : ',ShAB
      Write(Lunit,*) 'ShC     : ',ShC
      Write(Lunit,*) 'ShD     : ',ShD
      Write(Lunit,*) 'ShCD    : ',ShCD
      Write(Lunit,*) 'nColAB  : ',nColAB
      Write(Lunit,*) 'iOff_Col: ',(iOff_Col(i),i=1,8)
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'XThrCom    : ',XThrCom
      Write(Lunit,*) 'XThrDiag   : ',XThrDiag
      Write(Lunit,*) 'XDamp      : ',(XDamp(i),i=1,2)
      Write(Lunit,*) 'XSpan      : ',XSpan
      Write(Lunit,*) 'XThrNeg    : ',XThrNeg
      Write(Lunit,*) 'XWarNeg    : ',XWarNeg
      Write(Lunit,*) 'XTooNeg    : ',XTooNeg
      Write(Lunit,*) 'XnSym      : ',XnSym
      Write(Lunit,*) 'XnShell    : ',XnShell
      Write(Lunit,*) 'XnnShl     : ',XnnShl
      Write(Lunit,*) 'XnPass     : ',XnPass
      Write(Lunit,*) 'XCho_AdrVec: ',XCho_AdrVec
      Write(Lunit,*) 'XScDiag    : ',XScDiag
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Do j = 1,nChkQ+1
         Write(Lunit,*) 'iChkQ   : ',(iChkQ(i,j),i=1,4)
      End Do
      Write(Lunit,*) 'nCol_Chk: ',nCol_Chk
      Do j = 1,nSection
         Write(Lunit,*) 'TimSec  : ',(TimSec(i,j),i=1,4)
      End Do
      Do j = 1,nInteg
         Write(Lunit,*) 'tInteg  : ',(tInteg(i,j),i=1,2)
      End Do
      Do j = 1,nDecom
         Write(Lunit,*) 'tDecom  : ',(tDecom(i,j),i=1,2)
      End Do
      Write(Lunit,*) 'tDecDrv : ',(tDecDrv(i),i=1,2)
      Write(Lunit,*)
      Write(Lunit,*) 'nVecRS1: ',(nVecRS1(i),i=1,8)
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'Cho_AdrVec: ',Cho_AdrVec
      Write(Lunit,*) 'Cho_IOVec : ',Cho_IOVec
      Write(Lunit,*) 'nSys_Call : ',nSys_Call
      Write(Lunit,*) 'nDGM_Call : ',nDGM_Call
      Write(Lunit,*) 'N1_VecRd  : ',N1_VecRd
      Write(Lunit,*) 'N2_VecRd  : ',N2_VecRd
      Write(Lunit,*) 'N_Subtr   : ',N_Subtr
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'N1_Qual: ',N1_Qual
      Write(Lunit,*) 'N2_Qual: ',N2_Qual
      Write(Lunit,*)
      Write(Lunit,*) 'Frac_ChVBuf: ',Frac_ChVBuf
      Call Cho_Flush(Lunit)

C     choptr.fh.
C     -----------

      Call Cho_PrintPointers(irc,Lunit)
      If (irc .ne. 0) Return

C     chosubscr.fh.
C     --------------

      Write(Lunit,*) '*** Contents of chosubscr.fh:'
      Write(Lunit,*)
      Write(Lunit,*) 'Cho_SScreen: ',Cho_SScreen
      Write(Lunit,*) 'SSTau      : ',SSTau
      Write(Lunit,*) 'SubScrStat : ',(SubScrStat(i),i=1,2)
      Write(Lunit,*) 'DSubScr    : ',ip_DSubScr,l_DSubScr
      Write(Lunit,*) 'DSPNm      : ',ip_DSPNm,l_DSPNm

      End
      SubRoutine Cho_PrintPointers(irc,Lunit)
C
C     Purpose: print all entries in choptr.fh to Lunit.
C
      Implicit None
      Integer irc, Lunit
#include "choptr.fh"
      Integer nAlloc

      nAlloc = 0  ! allocation counter

      Write(Lunit,*) '*** Contents of choptr.fh:'
      Write(Lunit,*) '    (pointer,dimension)'
      Write(Lunit,*)
      Call Cho_Flush(Lunit)
      Write(Lunit,*) 'InfRed  : ',ip_InfRed,l_InfRed
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'InfVec  : ',ip_InfVec,l_InfVec
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'IndRed  : ',ip_IndRed,l_IndRed
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'IndRSh  : ',ip_IndRsh,l_IndRsh
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iScr    : ',ip_iScr,l_iScr
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iiBstRSh: ',ip_iiBstRSh,l_iiBstRSh
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'nnBstRSh: ',ip_nnBstRSh,l_nnBstRSh
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'IntMap  : ',ip_IntMap,l_IntMap
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'nDimRS  : ',ip_nDimRS,l_nDimRS
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iRS2F   : ',ip_iRS2F,l_iRS2F
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iSOShl  : ',ip_iSOShl,l_iSOShl
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iShlSO  : ',ip_iShlSO,l_iShlSO
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iQuab   : ',ip_iQuab,l_iQuab
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iBasSh  : ',ip_iBasSh,l_iBasSh
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'nBasSh  : ',ip_nBasSh,l_nBasSh
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'nBstSh  : ',ip_nBstSh,l_nBstSh
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iAtomShl: ',ip_iAtomShl,l_iAtomShl
      nAlloc    = nAlloc + 1
      Write(Lunit,*) 'iSP2F   : ',ip_iSP2F,l_iSP2F
      nAlloc    = nAlloc + 1
      Write(Lunit,*)
      Call Cho_Flush(Lunit)

      irc = CHO_NALLOC - nAlloc

      End
