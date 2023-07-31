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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine Cho_X_SetInc(irc)
C
C     T.B. Pedersen, July 2004.
C
C     Purpose: define all entries in include files
C              choprint.fh
C              choorb.fh
C              cholesky.fh
C              chosew.fh
C              chovecbuf.f90
C              chosubscr.fh
C              chpari.fh
C              cho_para_info.fh
C              and some in the Module choarr.f90
C
      use ChoArr, only: nDim_Batch, nQual_L, n_MySP
      use ChoBkm, only:  nRow_BkmVec, nCol_BkmVec,
     &                   nRow_BkmThr, nCol_BkmThr
      use ChoVecBuf, only: ip_CHVBUF_SYM, l_CHVBUF_SYM,
     &                     ip_CHVBFI_SYM, l_CHVBFI_SYM,
     &                     nVec_in_Buf
      use ChoSubScr, only: Cho_SScreen, SSTau, SubScrStat, SSNorm
      Implicit None
      Integer irc
#include "choorb.fh"
#include "choprint.fh"
#include "cholesky.fh"
#include "chpari.fh"
#include "cho_para_info.fh"

      Integer iLarge
      Parameter (iLarge = 99999999)

      Real*8 Large, Small
      Parameter (Large = 1.0D15, Small = 1.0D-15)

C     Set return code.
C     ----------------

      irc = 0

C     choprint.fh.
C     -------------

      iPrint = -iLarge

C     choorb.fh.
C     -----------

      Call iZero(iBas,8)
      Call iZero(nBas,8)
      Call iZero(XnBas,8)
      nBasT = 0

C     cholesky.fh.
C     -------------

      ThrCom  = Large
      ThrDiag = Large
      Tol_DiaChk = -Large
      ThrNeg  = Large
      WarNeg  = Large
      TooNeg  = Large
      nSym    = -iLarge
      lBuf    = -iLarge
      MinQual = iLarge
      MaxQUal = iLarge
      IFCSew  = -iLarge
      Mode_Screen = -iLarge
      Cho_DecAlg  = -iLarge
      Cho_DecAlg_Def = -iLarge
      iAlQua  = iLarge
      MxShPr  = iLarge
      ModRst  = iLarge
      Run_Mode = iLarge
      ScDiag  = .false.
      ChkOnly = .false.
      Cho_IntChk = .false.
      Cho_MinChk = .false.
      Cho_UseAbs = .false.
      Cho_TrcNeg = .false.
      Cho_ReOrd  = .false.
      Cho_DiaChk = .false.
      Cho_TstScreen = .false.
      Cho_1Center = .false.
      Cho_No2Center = .false.
      Cho_PreScreen = .false.
      Cho_SimP = .false.
      Cho_Fake_Par = .false.
      RstDia  = .false.
      RstCho  = .false.
      Did_DecDrv = .false.
      HaltIt  = .false.
      Trace_Idle = .false.

      Call iZero(LuCho,8)
      Call iZero(LuSel,8)
      Call iZero(LuTmp,8)
      LuPri = 0
      LuScr = 0
      LuRed = 0
      LuRst = 0
      LuMap = 0

      nShell = 0
      nnShl_Tot  = 0
      nnShl = 0
      MxORSh = 0
      Mx2Sh  = 0
      Call iZero(iiBstR,8*3)
      Call iZero(nnBstR,8*3)
      Call iZero(nnBstRT,3)
      mmBstRT = 0
      nQual_L(:)=0
      Call iZero(iOffQ,8)

      Call FZero(DiaMax,8)
      Call FZero(DiaMaxT,8)
      Call FZero(DiaMin,8)
      Call FZero(Damp,2)
      Span   = Large
      XlDiag = Large
      DiaMnZ = Large
      Thr_PreScreen = -Large
      iABMnZ = -iLarge
      nnZTot = 0

      Call iZero(NumCho,8)
      NumChT = 0
      MaxVec = 0
      MaxRed = 0
      BlockSize = -iLarge

      ShA = -iLarge
      ShB = -iLarge
      ShC = -iLarge
      ShD = -iLarge
      ShAB = -iLarge
      ShCD = -iLarge
      nColAB = -iLarge
      Call iZero(iOff_Col,8)

      XThrCom  = Large
      XThrDiag = Large
      Call FZero(XDamp,2)
      XSpan    = Large
      XThrNeg  = Large
      XWarNeg  = Large
      XTooNeg  = Large
      XnSym    = 0
      XnShell  = 0
      XnnShl   = 0
      XnPass   = 0
      XScDiag  = .false.
      XCho_AdrVec = -iLarge

      Call iZero(iChkQ,4*(nChkQ+1))
      nCol_Chk = -iLarge
      Call FZero(TimSec,4*nSection)
      Call FZero(tInteg,2*nInteg)
      Call FZero(tDecom,2*nDecom)
      Call FZero(tMisc,2*nMisc)
      Call FZero(tDecDrv,2)

      Call iZero(nVecRS1,8)

      Cho_AdrVec= -iLarge
      Cho_IOVec = -iLarge
      nSys_Call = 0
      nDGM_Call = 0
      N1_VecRd  = 0
      N2_VecRd  = 0
      N_Subtr   = 0

      N1_Qual = -iLarge
      N2_Qual = iLarge

      Frac_ChVBuf = 0.0d0

      nDim_Batch(:)=0

      nQual_L(:)=0

      n_MySP=0

      Cho_SimRI = .false.
      Thr_SimRI = -Large

C     chovecbuf.f90.
C     --------------

      Call iZero(ip_ChVBuf_Sym,8)
      Call iZero(l_ChVBuf_Sym,8)
      Call iZero(ip_ChVBfI_Sym,8)
      Call iZero(l_ChVBfI_Sym,8)
      Call iZero(nVec_in_Buf,8)

C     chosubscr.fh.
C     --------------

      Cho_SScreen = .false.
      SSTau       = 0.0d0
      SubScrStat(1) = 0.0d0
      SubScrStat(2) = 0.0d0
      SSNorm      = 'tbp'

C     chpari.fh.
C     -----------

      Call iZero(NumCho_Bak,8)

C     cho_para_info.fh.
C     ------------------

      Cho_Real_Par = .false.

C     chobkm.f90
C     -----------

      nRow_BkmVec=0
      nCol_BkmVec=0
      nRow_BkmThr=0
      nCol_BkmThr=0

      End
