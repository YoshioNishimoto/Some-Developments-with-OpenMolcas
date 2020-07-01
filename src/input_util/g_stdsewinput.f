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
      Subroutine StdSewInput(Info,nInfo,LuRd,ifnr,mdc,iShll,BasisTypes,
     &                       STDINP,lSTDINP,iErr,DInf,nDInf)
************************************************************************
* This is a simplified copy of the BASI section of RdCtl_Seward that   *
* reads the string vector STDINP with the standard seward input        *
* generated by ZMatrixConverter.                                       *
************************************************************************
      use Basis_Info
      Implicit Real*8 (a-h,o-z)
*
#include "itmax.fh"
#include "info.fh"
#include "SysDef.fh"
#include "rctfld.fh"
#include "real.fh"
#include "print.fh"
#include "gateway.fh"
#include "stdalloc.fh"
c
c     IRELAE = 0  .... DKH
c            = 1  .... DK1
c            = 2  .... DK2
c            = 3  .... DK3
c            = 4  .... DK3full
c            = 11 .... RESC
c            = 21 .... ZORA
c            = 22 .... ZORA(FP)
c            = 23 .... IORA
CAWMR
c     NB: The IRELAE flag has been extended to account for
c         arbitrary-order DKH with different parametrizations!
c         IMPORTANT: new arbitrary-order DKH routines are only
c                    called for IRELAE values LARGER than 1000.
CAWMR
c
#include "relae.fh"
      Integer, Parameter:: nBuff=10000
      Real*8, Allocatable:: Buffer(:)
      Real*8 DInf(nDInf)
      Common /AMFn/ iAMFn
      Common /delete/ kDel(0:MxAng,MxDc)
*
      Character Key*180, KWord*180,            BSLbl*80, Fname*256,
     &          DefNm*13, Ref(2)*80, dbas*4
      Integer StayAlone, nDel(MxAng), BasisTypes(4)
*
      Character*180 Line, STDINP(mxAtom*2) ! CGGn
      Character*256 Basis_lib ! CGGd , INT2CHAR, CHAR4
      Character*5  Symbol                 ! CGGn
*
CGGd      Data WellRad/-1.22D0,-3.20D0,-6.20D0/
      Data StayAlone/0/
*
#include "angstr.fh"
      Data DefNm/'basis_library'/ ! CGGd,
*                                                                      *
************************************************************************
*                                                                      *
#include "getbs_interface.fh"
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Buffer,nBuff,Label='Buffer')
      iErr=0
      iRout=3
*                                                                      *
************************************************************************
*                                                                      *
      LuWr=6
*
      imix=0
      itype=0
*
      ipExp(1) = Info
      BasisTypes(1)=0
      BasisTypes(2)=0
      BasisTypes(3)=0
      BasisTypes(4)=0
*                                                                      *
****** BASI ************************************************************
*                                                                      *
      iSTDINP = 2

10    nCnttp = nCnttp + 1
      If (nCnttp.gt.Mxdbsc) Then
         Write (LuWr,*) ' Increase Mxdbsc'
         iErr=1
         Return
      End If
*
*     Read the basis set label
*
      Key = STDINP(iSTDINP)
      BSLbl = Key(1:80)
*     Call UpCase(BSLbl)
      LenBSL=Len(BSLbl)
      Last=iCLast(BSLbl,LenBSL)
      Indx=Index(BSLbl,'/')
      If (Indx.eq.0) Then
       call WhichMolcas(Basis_lib)
       if(Basis_lib(1:1).ne.' ') then
         StayAlone=1
         ib=index(Basis_lib,' ')-1
         if(ib.lt.1)
     *    Call SysAbendMsg('rdCtl','Too long PATH to MOLCAS',' ')
         Fname=Basis_lib(1:ib)//'/basis_library'
       else
         Fname=DefNm
       endif
       Indx = Last+1
       Bsl(nCnttp)=BSLbl
      Else
         Fname= BSLbl(Indx+2:Last)
         If (Fname.eq.' ') Then
            Call WarningMessage(2,
     &                     ' No basis set library specified for'
     &                   //';BSLbl='//BSLbl//';Fname='//Fname)
            Call Quit_OnUserError()
         End If
 1919    If (Fname(1:1).eq.' ') Then
            Fname(1:79)=Fname(2:80)
            Fname(80:80) = ' '
            Go To 1919
         End If
         Bsl(nCnttp)=BSLbl(1:Indx-1)
      End If
*
      n=INDEX(Bsl(nCnttp),' ')
      Bsl(nCnttp)(n:n+5)='.....'
*
      If (Show.and.nPrint(2).ge.6) Then
         Write (LuWr,*)
         Write (LuWr,*)
         Write(LuWr,'(1X,A,I5,A,A)')
     &           'Basis Set ',nCnttp,' Label: ', BSLbl(1:Indx-1)
         Write(LuWr,'(1X,A,A)') 'Basis set is read from library:',
     *         Fname(1:index(Fname,' '))
      End if
*
      jShll = iShll
      SODK(nCnttp)=.False.
      AuxCnttp(nCnttp)=.False.
      Bsl_Old(nCnttp)=Bsl(nCnttp)
      mdciCnttp(nCnttp)=mdc
      Call GetBS(Fname,Bsl(nCnttp),Indx-1,lAng,ipExp,
     &           ipCff,ipCff_Cntrct,ipCff_Prim,ipFockOp,
     &           nExp,nBasis,nBasis_Cntrct,MxShll,iShll,
     &           MxAng,Charge(nCnttp),
     &           iAtmNr(nCnttp),BLine,Ref, PAM2(nCnttp),
     &           ipPAM2xp(nCnttp),ipPAM2cf(nCnttp),nPAM2(nCnttp),
     &           FockOp(nCnttp),
     &           ECP(nCnttp),NoPairL(nCnttp),SODK(nCnttp),
     &           ipM1xp(nCnttp),ipM1cf(nCnttp),nM1(nCnttp),
     &           ipM2xp(nCnttp),ipM2cf(nCnttp),nM2(nCnttp),ipBk,
     &           CrRep(nCnttp),nProj,nAIMP,ipAkl,ip_Occ,iOptn,
     &           UnNorm,nDel,
     &            nVal,   nPrj,   nSRO,   nSOC,  nPP,
     &           ipVal_, ipPrj_, ipSRO_, ipSOC_,ipPP_,
     &           LuRd,BasisTypes,AuxCnttp(nCnttp),
     &           nFragType(nCnttp),nFragCoor(nCnttp),nFragEner(nCnttp),
     &           nFragDens(nCnttp),ipFragType(nCnttp),ipFragCoor(nCnttp)
     &           ,ipFragEner(nCnttp),ipFragCoef(nCnttp),IsMM(nCnttp),
     &           STDINP,iSTDINP,.True.,.true.,' ',
     &           DInf,nDInf,nCnttp)
*
      Do_FckInt = Do_FckInt .and. FockOp(nCnttp)
      If (itype.eq.0) Then
         If (BasisTypes(3).eq.1 .or. BasisTypes(3).eq.2)
     &       iType=BasisTypes(3)
      Else
         If (BasisTypes(3).eq.1 .or. BasisTypes(3).eq.2) Then
            If (BasisTypes(3).ne.iType) Then
               imix=1
               BasisTypes(3)=-1
            End If
            iType=BasisTypes(3)
         End If
      End If
      If (itype.eq.1) ifnr=1
      If (itype.eq.2) ifnr=0
*
      If (nSOC.gt.-1) Then
         Do l = 1, MxAng
            kDel(l,nCnttp)=nDel(l)
         End Do
      End If
      If (Show.and.nPrint(2).ge.6 .and.
     &   Ref(1).ne.BLine .and. Ref(2).ne.Bline) Then
         Write (LuWr,'(1x,a)')  'Basis Set Reference(s):'
         If (Ref(1).ne.BLine) Write (LuWr,'(5x,a)') Ref(1)
         If (Ref(2).ne.BLine) Write (LuWr,'(5x,a)') Ref(2)
         Write (LuWr,*)
         Write (LuWr,*)
      End If
      lPAM2 = lPAM2 .or. PAM2(nCnttp)
      ECP(nCnttp)=(nPP+nPrj+nSRO+nSOC+nM1(nCnttp)+nM2(nCnttp)).ne.0
      lPP=lPP .or. nPP.ne.0
      lECP = lECP .or. ECP(nCnttp)
      lNoPair = lNoPair .or. NoPairL(nCnttp)
*
      iAngMx=Max(iAngMx,lAng)
*     No transformation needed for s and p shells
      Transf(jShll+1)=.False.
      Prjct(jShll+1)=.False.
      Transf(jShll+2)=.False.
      Prjct(jShll+2)=.False.
      pChrg(nCnttp)=.False.
      Fixed(nCnttp)=.False.
      nOpt(nCnttp) = iOptn
      ipVal(nCnttp) = ipVal_
      ipPrj(nCnttp) = ipPrj_
      ipSRO(nCnttp) = ipSRO_
      ipSOC(nCnttp) = ipSOC_
      ipPP(nCnttp)  = ipPP_
      nVal_Shells(nCnttp) = nVal
      nPrj_Shells(nCnttp) = nPrj
      nSRO_Shells(nCnttp) = nSRO
      nSOC_Shells(nCnttp) = nSOC
      nPP_Shells(nCnttp)  = nPP
      nTot_Shells(nCnttp) = nVal+nPrj+nSRO+nSOC+nPP
      nCnt = 0
      lAux = lAux .or. AuxCnttp(nCnttp)
      If (AuxCnttp(nCnttp)) Then
         Do iSh = jShll+1, iShll
            AuxShell(iSh)=.True.
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the effective radius of this center
*
      iAng = 0
      Thrshld_R=1.0D-08
      Do iSh = ipVal_, ipVal_+nVal-1
         RMax_R=Zero
         Do iPrim = 0, nExp(iSh)-1
            ValExp = DInf(ipExp(iSh)+iPrim)
            RMax_R = Max(RMax_R,
     &                   Eval_RMax(ValExp,iAng,Thrshld_R))
         End Do
         RMax_Shll(iSh)=RMax_R
C        Write (LuWr,*) 'RMax_R=',RMax_R
         iAng = iAng + 1
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Here we will have to fix that the 6-31G family of basis sets
*     should by default be used with 6 d-functions rather than 5.
*
      KWord=BSLbl(1:Indx-1)
      Call UpCase(KWord)
      If (INDEX(KWord,'6-31G').ne.0) Then
         Do iSh = jShll+3, iShll
            Prjct(iSh)=.False.
            Transf(iSh)=.False.
         End Do
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
100   iSTDINP = iSTDINP + 1
      Line = STDINP(iSTDINP)
      KWord = Line
      Call UpCase(KWord)
      If (KWord(1:4).eq.'END ') Then
         If (nCnt.eq.0) Then
            Call WarningMessage(2,' Input error, no center specified!')
            Call Quit_OnUserError()
         End If
         dbsc(nCnttp)%nCntr = nCnt
!        call allocate(dbsc(nCnttp)%Coor(1:3,nCnt))
         call mma_allocate(dbsc(nCnttp)%Coor,3,nCnt,Label='dbsc:C')
         Call DCopy_(3*nCnt,Buffer,1,dbsc(nCnttp)%Coor(1,1),1)
         mdc = mdc + nCnt
*        Compute the number of elements stored in the dynamic memory
*        so far.
         nInfo = ipExp(iShll+1) - Info
* the next line seems to convince IBM XLF 6.1 to forgo its otherwise
* crass behaviour. Who can tell why? Peter Knowles, 7/99
         ninfo_stupid = nInfo
         Go To 900
      End If
*
*     Read Coordinates
*
      nCnt = nCnt + 1
      If (mdc+nCnt.gt.Mxdc) Then
         Call WarningMessage(2,' RdCtl: Increase Mxdc')
         Write (LuWr,*) '        Mxdc=',Mxdc
         Call Quit_OnUserError()
      End If
      iend=Index(KWord,' ')
      If (iEnd.gt.LENIN+1) Then
         Write (6,*) 'Warning: the label ', KWord(1:iEnd),
     &               ' will be truncated to ',LENIN,' characters!'
      End If
      LblCnt(mdc+nCnt) = KWord(1:Min(LENIN,iend-1))
      dbas=Trim(LblCnt(mdc+nCnt)(1:LENIN))
      Call Upcase(dbas)
      If (dbas.eq.'DBAS') Then
         Call WarningMessage(2,' RdCtl: ZMAT does not work with DBAS')
         Call Quit_OnUserError()
      End If
      If (mdc+nCnt.gt.1)
     &   Call ChkLbl(LblCnt(mdc+nCnt),LblCnt,mdc+nCnt-1)
      iOff=1+(nCnt-1)*3
C      print *,line
      Read (Line,'(A5)') Symbol
      Read (Line(6:),*) (Buffer(iOff+i),i=0,2) ! CGGn
      If (Index(KWord,'ANGSTROM').ne.0) Then
         Do i = 0, 2
            Buffer(iOff+i) = Buffer(iOff+i)/angstr
         End Do
      End If
*
      GoTo 100

900   iSTDINP = iSTDINP + 2
      If (iSTDINP.LT.lSTDINP) Go to 10

      If (iAngMx.lt.0) Then
         Write (6,*) ' There is an error somewhere in the input!'
         Write (6,*) 'iAngMx.lt.'
         iErr=1
         Return
      End If
      Call mma_deallocate(Buffer)
*
      Return
      End
