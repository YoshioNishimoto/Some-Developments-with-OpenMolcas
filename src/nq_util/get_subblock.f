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
* Copyright (C) 1999, Roland Lindh                                     *
************************************************************************
      Subroutine Get_Subblock(Kernel,Func,ixyz,
     &                        Maps2p,list_s,list_exp,list_bas,
     &                        nShell,nSym, list_p,R2_trial,nNQ,
     &                        FckInt,nFckDim,nFckInt,nD,
     &                        mGrid,nP2_ontop,Do_Mo,
     &                        P2Unzip,D1mo,D1Unzip,nD1mo,P2_ontop,
     &                        Do_Grad,Grad,nGrad,List_G,IndGrd,iTab,
     &                        Temp,mGrad,mAO,mdRho_dR,
     &                        LTEG_DB,PDFTPot1,PDFTFocI,PDFTFocA)
************************************************************************
*                                                                      *
* Object: to generate the list of the shell and exponent that have an  *
*         influence on a subblock                                      *
*                                                                      *
*     Called from: Drvnq_                                              *
*                                                                      *
*     Author: Roland Lindh,                                            *
*             Dept of Chemical Physics,                                *
*             University of Lund, Sweden                               *
*             August 1999                                              *
************************************************************************
      use iSD_data
      use Basis_Info
      use Center_Info
      use nq_Grid, only: Grid, Weights, TabAO, Grid_AO, Dens_AO,
     &                   TabAO_Pack, Ind_Grd, dRho_dR, TabAO_Short,
     &                   kAO, iBfn_Index
      use nq_MO, only: DoIt
      Implicit Real*8 (A-H,O-Z)
      External Kernel
#include "itmax.fh"
#include "nq_info.fh"
#include "Molcas.fh"
#include "nsd.fh"
#include "setup.fh"
#include "real.fh"
#include "grid_on_disk.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "debug.fh"
#include "ksdft.fh"
      Integer Maps2p(nShell,0:nSym-1), list_s(2,*), List_G(3,*),
     &        list_exp(nSym*nShell), list_bas(2,nSym*nShell),
     &        list_p(nNQ), iTab(4,mGrad),IndGrd(mGrad)
      Real*8 R2_trial(nNQ), FckInt(nFckInt,nFckDim),
     &       Grad(nGrad), Temp(mGrad), D1mo(nD1mo),
     &       P2_ontop(nP2_ontop,mGrid), Roots(3,3),
     &       xyz0(3,2),PDFTPot1(npot1),PDFTFocI(nPot1),PDFTFocA(nPot1)
      Logical InBox(MxAtom), Do_Grad, More_to_come
      Logical Do_Mo
      Integer LTEG_DB
      Real*8,DIMENSION(NASHT4)::P2Unzip
      Real*8,DIMENSION(NASHT**2)::D1Unzip
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
#include "nq_structure.fh"
      declare_ip_coor
      declare_ip_r_max
      declare_ip_angular
      iGridInfo(i,iNQ)=iWork(ip_GridInfo+(iNQ-1)*2+i-1)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (Debug) Then
         Write(6,*) 'Enter Get_Subblock'
         Write(6,*) 'ip_nR_Eff GET_SBK',ip_nR_Eff
      End If
#endif
*
*-----Resolve triplet index
*
      iyz = 1 + (ixyz-1)/nx
      ix  = ixyz - (iyz-1)*nx
      iz  = 1 + (iyz-1)/ny
      iy  = iyz - (iz-1)*ny
*
*-----Get the extreme coordinates of the box.
*
      x_min_= x_min+DBLE(ix-2)*Block_Size
      x_max_= x_min_+Block_Size
      y_min_= y_min+DBLE(iy-2)*Block_Size
      y_max_= y_min_+Block_Size
      z_min_= z_min+DBLE(iz-2)*Block_Size
      z_max_= z_min_+Block_Size
      If (ix.eq.1 ) x_min_=-1.0D99
      If (ix.eq.nx) x_max_= 1.0D99
      If (iy.eq.1 ) y_min_=-1.0D99
      If (iy.eq.ny) y_max_= 1.0D99
      If (iz.eq.1 ) z_min_=-1.0D99
      If (iz.eq.nz) z_max_= 1.0D99
*                                                                      *
#ifdef _DEBUGPRINT_
      If (Debug) Then
         Write(6,*)
         Write(6,*) 'Block_Size=',Block_Size
         Write(6,*) 'ix,iy,iz=',ix,iy,iz
         Write(6,*) 'x_min_,x_max_',x_min_,x_max_
         Write(6,*) 'y_min_,y_max_',y_min_,y_max_
         Write(6,*) 'z_min_,z_max_',z_min_,z_max_
         Write(6,*) 'nNQ=',nNQ
      End If
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----Generate list over atoms which contribute to the subblock.
*                                                                      *
************************************************************************
*                                                                      *
      ilist_p=0
      Do 10 iNQ=1,nNQ
         InBox(iNQ)=.False.
*--------Get the coordinates of the partitionning
         x_NQ =Work(ip_Coor(iNQ)  )
         y_NQ =Work(ip_Coor(iNQ)+1)
         z_NQ =Work(ip_Coor(iNQ)+2)
*
*        1) center is in the box
*
         If   ((x_NQ.ge.x_min_).and.(x_NQ.le.x_max_)
     &    .and.(y_NQ.ge.y_min_).and.(y_NQ.le.y_max_)
     &    .and.(z_NQ.ge.z_min_).and.(z_NQ.le.z_max_))
     &     InBox(iNQ)=.True.
         If (InBox(iNQ)) Then
            ilist_p=ilist_p+1
            list_p(ilist_p)=iNQ
            GoTo 10
         EndIf
*
*        2) atomic grid of this center extends inside the box.
*
         RMax = Work(ip_R_Max(iNQ))
         t1=(x_NQ-x_min_)/(x_max_-x_min_)
         If (t1.lt.Zero) t1=Zero
         If (t1.gt.One ) t1=One
         t2=(y_NQ-y_min_)/(y_max_-y_min_)
         If (t2.lt.Zero) t2=Zero
         If (t2.gt.One ) t2=One
         t3=(z_NQ-z_min_)/(z_max_-z_min_)
         If (t3.lt.Zero) t3=Zero
         If (t3.gt.One ) t3=One
         R2_Trial(iNQ)= (x_NQ-(x_max_-x_min_)*t1-x_min_)**2
     &                + (y_NQ-(y_max_-y_min_)*t2-y_min_)**2
     &                + (z_NQ-(z_max_-z_min_)*t3-z_min_)**2
         If (R2_Trial(iNQ).le.RMax**2) Then
            ilist_p=ilist_p+1
            list_p(ilist_p)=iNQ
         EndIf
 10   Continue
      nlist_p=ilist_p
      If (nlist_p.eq.0) return
#ifdef _DEBUGPRINT_
      If (debug) Write (6,*) 'Get_Subblock: List_p:',List_p
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----Generate list over shells which contribute to the subblock.
*                                                                      *
************************************************************************
*                                                                      *
      ilist_s=0
      Do iShell=1,nShell
#ifdef _DEBUGPRINT_
         If (debug) Write (6,*) 'iShell,nShell=',iShell,nShell
#endif
         NrExp =iSD( 5,iShell)
         iAng  =iSD( 1,iShell)
         iShll =iSD( 0,iShell)
         NrBas =iSD( 3,iShell)
         mdci  =iSD(10,iShell)
         nDegi=nSym/dc(mdci)%nStab
*
         Do jSym = 0, nDegi-1
            iSym=dc(mdci)%iCoSet(jSym,0)
#ifdef _DEBUGPRINT_
            If (debug) Write (6,*) 'iSym,nDegi-1=',iSym,nDegi-1
#endif
*
            iNQ=Maps2p(iShell,NrOpr(iSym))
            RMax_NQ = Work(ip_R_Max(iNQ))
#ifdef _DEBUGPRINT_
            If (debug) Then
               Write (6,*) 'iNQ=',iNQ
               Write (6,*) 'RMax_NQ=',RMax_NQ
               Write (6,*) 'InBox(iNQ)=',InBox(iNQ)
            End If
#endif
*
*           1) the center of this shell is inside the box
*
            If (InBox(iNQ)) Then
               ilist_s=ilist_s+1
               list_s(1,ilist_s)=iShell
               list_s(2,ilist_s)=iSym
               list_exp(ilist_s)=NrExp
               list_bas(1,ilist_s)=NrBas
               GoTo 20
            End If
#ifdef _DEBUGPRINT_
            If (debug) Write (6,*) 'Passed here!'
#endif
*
*           2) the Gaussian has a grid point which extends inside the
*              box. The Gaussians are ordered from the most diffuse to
*              the most contracted.
*
            nExpTmp=0
            Do iExp=1,NrExp
*------------- Get the value of the exponent
               ValExp=Shells(iShll)%Exp(iExp)
*------------- If the exponent has an influence then increase the
*              number of actives exponents for this shell, else
*              there is no other active exponent (they ar ordered)
               RMax=Min(Eval_RMax(ValExp,iAng,Threshold),RMax_NQ)
#ifdef _DEBUGPRINT_
               If (Debug) Then
                  Write (6,*) 'iShell,iNQ=',iShell,iNQ
                  Write (6,*) 'ValExp=',ValExp
                  Write (6,*) 'RMax_NQ=',RMax_NQ
                  Write (6,*) 'RMax_Exp=',
     &                         Eval_RMax(ValExp,iAng,Threshold)
                  Write (6,*) 'RMax=',RMax
                  Write (6,*) 'R2_Trial(iNQ),RMax**2=',
     &                         R2_Trial(iNQ),RMax**2
               End If
#endif
               If (R2_Trial(iNQ).gt.RMax**2) Go To 99
               nExpTmp=nExpTmp+1
            End Do              !iExp
 99         Continue
            If (nExpTmp.ne.0) Then
               ilist_s=ilist_s+1
               list_s(1,ilist_s)=iShell
               list_s(2,ilist_s)=iSym
               list_exp(ilist_s)=nExpTmp
*
*              Examine if contracted basis functions can be ignored.
*              This will be the case for segmented basis sets.
*
               list_bas(1,ilist_s)=nBas_Eff(NrExp,NrBas,
     &                                      Shells(iShll)%Exp,
     &                                      Shells(iShll)%pCff,
     &                                      list_exp(ilist_s))
            End If
 20         Continue
         End Do ! iSym
      End Do    ! iShell
      nlist_s=ilist_s
#ifdef _DEBUGPRINT_
      If (Debug) Write (6,*) 'nList_s,nList_p=',nList_s,nList_p
#endif
      If (nList_s*nList_p.eq.0) return
*                                                                      *
************************************************************************
*                                                                      *
*     Generate index arrays to address the density matrix, which is in
*     the full shell and the reduced shell over which in the basis
*     functions will be evaluated.
*
      nIndex=0
      Do ilist_s=1, nlist_s
         iShell=list_s(1,ilist_s)
         NrBas_Eff=list_bas(1,ilist_s)
         iCmp  = iSD( 2,iShell)
         nIndex=nIndex + NrBas_Eff*iCmp
      End Do
*
      Call GetMem('Index','Allo','Inte',ipIndex,nIndex)
*
      iIndex=1
      nAOs=0
      nAOs_Eff=0
      Do ilist_s=1, nlist_s
         iShell=list_s(1,ilist_s)
         NrBas =iSD( 3,iShell)
         NrBas_Eff=list_bas(1,ilist_s)
         iCmp  = iSD( 2,iShell)
         nAOs=nAOs+NrBas*iCmp
         nAOs_Eff=nAOs_Eff+NrBas_Eff*iCmp
         list_bas(2,ilist_s)=iIndex
         Call Do_Index(iWork(ipIndex-1+iIndex),NrBas,NrBas_Eff,iCmp)
         iIndex=iIndex + NrBas_Eff*iCmp
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (Debug) Then
         write(6,*) 'Contribution to the subblock :'
         write(6,*) 'NQ :',(list_p(ilist_p)  ,ilist_p=1,nlist_p)
         write(6,*) 'Sh :',(list_s(1,ilist_s),ilist_s=1,nlist_s)
         write(6,*) '   :',(list_s(2,ilist_s),ilist_s=1,nlist_s)
         write(6,*) 'Exp:',(list_exp(ilist_s),ilist_s=1,nlist_s)
      End If
#endif
*
      nBfn=0
      Do iList_s = 1, nList_s
         iSkal    =list_s(1,ilist_s)
         NrBas_Eff=list_bas(1,ilist_s)
         iCmp  = iSD( 2,iSkal)
         nBfn=nBfn+NrBas_Eff*iCmp
      End Do
*
      Call Allocate_iWork(ipTabAO,2*(nlist_s+1))
      Call mma_Allocate(Dens_AO,nBfn,nBfn,nD,Label='Dens_AO')
      Call mma_Allocate(Ind_Grd,3,nBfn,Label='Ind_Grd')
      Call mma_Allocate(iBfn_Index,6,nBfn,Label='iBfn_Index')
      Ind_Grd(:,:)=0
      iBfn_Index(:,:)=0
*
      If ((Functional_Type.eq.CASDFT_Type).or.Do_MO) Then
         nMOs=SIZE(DoIt)
         nTabMO=mAO*nMOs*mGrid
         Call Allocate_Work(ipTabMO,nTabMO)
         nTabSO=mAO*nMOs*mGrid
         Call Allocate_Work(ipTabSO,nTabSO)
      Else
         ipTabMO=ip_Dummy
         ipTabSO=ip_Dummy
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Generate indexation of which shells contributes to which centers
*     and center index for each gradient contribution which is computed.
*
      nGrad_Eff=0
      If (Do_Grad) Then
         Call ICopy(3*nShell*nSym,[0],0,List_G,1)
         Do ilist_s = 1, nlist_s
            iShell=list_s(1,ilist_s)
            iSym  =list_s(2,ilist_s)
            mdci  =iSD(10,iShell)
            iNQ = Maps2p(iShell,NrOpr(iSym))
            Do iCar=0,2
               If ((iSD(16+iCar,iShell).ne.0 .or.
     &              iSD(12,iShell).eq.1) .and.
     &             List_G(1+iCar,ilist_s).eq.0) Then
                  nGrad_Eff=nGrad_Eff+1
*
*                 For pseudo centers note that there will not be a
*                 gradient computed for this center.
*
                  iPseudo=iSD(12,iShell)
                  If (iPseudo.eq.0) Then
                     IndGrd(nGrad_Eff)=iSD(16+iCar,iShell)
                  Else
                     IndGrd(nGrad_Eff)=-1
                  End If
                  List_G(1+iCar,ilist_s)=nGrad_Eff
                  iTab(1,nGrad_Eff)=iCar+1
                  iTab(3,nGrad_Eff)=iNQ
                  kNQ=Maps2p(iShell,0)
                  Xref=Work(ip_Coor(kNQ)+iCar)
                  X   =Work(ip_Coor(iNQ)+iCar)
                  If (X.eq.Xref) Then
                     iTab(4,nGrad_Eff)=dc(mdci)%nStab
                  Else
                     iTab(4,nGrad_Eff)=-dc(mdci)%nStab
                  End If
*
*---------------- Find all other shells which contribute to the same
*                 gradient.
*
                  Do jlist_s = ilist_s+1, nlist_s
                     jShell=list_s(1,jlist_s)
                     If ( (iSD(16+iCar,iShell).eq.
     &                     iSD(16+iCar,jShell))
     &                    .and.
     &                    (iSym.eq.list_s(2,jlist_s))
     &                  ) Then
                        List_G(1+iCar,jlist_s)=nGrad_Eff
                     End If
                  End Do
*
*
*------- Include derivatives which will be used for
*        the translational invariance equation but which do not
*        contribute directly to a symmetry adapted gradient
*
               Else If (iSD(16+iCar,iShell).eq.0 .and.
     &                 List_G(1+iCar,ilist_s).eq.0) Then
                       nGrad_Eff=nGrad_Eff+1
                       IndGrd(nGrad_Eff)=-1
                       List_G(1+iCar,ilist_s)=nGrad_Eff
                       iTab(1,nGrad_Eff)=iCar+1
                       iTab(3,nGrad_Eff)=iNQ
                       iTab(4,nGrad_Eff)=dc(mdci)%nStab
*
*--------------------- Find all other shells which contribute to the same
*                      gradient.
*
                       Do jlist_s = ilist_s+1, nlist_s
                          jShell=list_s(1,jlist_s)
                          jSym  =list_s(2,jlist_s)
                          jNQ = Maps2p(jShell,NrOpr(jSym))
                          If (iNQ.eq.jNQ) Then
                             List_G(1+iCar,jlist_s)=nGrad_Eff
                          End If
                       End Do
               End If
            End Do
         End Do
*
         If (Grid_Type.eq.Moving_Grid) Then
            ndW_dR=nGrad_Eff*mGrid
            Call GetMem('dW_dR','Allo','Real',ip_dW_dR,ndW_dR)
            ndW_Temp=3*nlist_p
            Call GetMem('dW_Temp','Allo','Real',ip_dW_Temp,ndW_Temp)
            ndPB=3*nlist_p**2
            Call GetMem('dPB','Allo','Real',ip_dPB,ndPB)
         Else
            ip_dW_dR  =ip_Dummy
            ip_dW_Temp=ip_Dummy
            ip_dPB    =ip_Dummy
         End If
      Else
         ip_dW_dR  =ip_Dummy
         ip_dW_Temp=ip_Dummy
         ip_dPB    =ip_Dummy
      End If
      If (Do_Grad.and.nGrad_Eff.eq.0) Go To 998
      If (Grid_Status.eq.Use_Old) Go To 997
         Call ICopy(3*nBatch_Max,[0],0,iBatchInfo,1)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*-----For each partition active in the subblock create the grid
*     and perform the integration on it.
*                                                                      *
************************************************************************
*                                                                      *
      number_of_grid_points=0
      nBatch = 0
      Do ilist_p=1,nlist_p
         iNQ=list_p(ilist_p)
#ifdef _DEBUGPRINT_
         If (Debug) Write (6,*) 'ilist_p=',ilist_p
         If (Debug) Write (6,*) 'Get_SubBlock: iNQ=',iNQ
#endif
*
*------- Select which gradient contributions that should be computed.
*        For basis functions which have the center common with the grid
*        do not compute any contribution.
*
         If (Do_Grad) Then
            Call ICopy(nGrad_Eff,[On],0,iTab(2,1),4)
            If (Grid_Type.eq.Moving_Grid) Then
               Do iGrad = 1, nGrad_Eff
                  jNQ = iTab(3,iGrad)
                  If (iNQ.eq.jNQ) iTab(2,iGrad)=Off
               End Do
            End If
#ifdef _DEBUGPRINT_
            If (Debug) Then
             Write (6,*)
             Write (6,'(A,24I3)') '       i =',(       i ,i=1,nGrad_Eff)
             Write (6,'(A,24I3)') 'iTab(1,i)=',(iTab(1,i),i=1,nGrad_Eff)
             Write (6,'(A,24I3)') 'iTab(2,i)=',(iTab(2,i),i=1,nGrad_Eff)
             Write (6,'(A,24I3)') 'iTab(3,i)=',(iTab(3,i),i=1,nGrad_Eff)
             Write (6,'(A,24I3)') 'iTab(4,i)=',(iTab(4,i),i=1,nGrad_Eff)
             Write (6,*) 'IndGrd=',IndGrd
             Write (6,*)
            End If
#endif
*
         End If
*
*--------Get the coordinates of the partition
         x_NQ =Work(ip_Coor(iNQ)  )
         y_NQ =Work(ip_Coor(iNQ)+1)
         z_NQ =Work(ip_Coor(iNQ)+2)
*--------Get the maximum radius on which we have to integrate for the
*        partition
         RMax=Work(ip_R_Max(iNQ))
*
         Call Box_On_Sphere(x_Min_-x_NQ,x_Max_-x_NQ, y_Min_-y_NQ,
     &                      y_Max_-y_NQ, z_Min_-z_NQ,z_Max_-z_NQ,
     &                      xyz0(1,1), xyz0(1,2),
     &                      xyz0(2,1), xyz0(2,2),
     &                      xyz0(3,1), xyz0(3,2))
*                                                                      *
************************************************************************
*                                                                      *
*------- Establish R_Box_Max and R_Box_Min, the longest and the shortest
*        distance from the origin of the atomic grid to a point in the
*        box
*
         R_Box_Max=Zero
         R_Box_Min=RMax
*
         x_box_min = x_min_ - x_NQ
         x_box_max = x_max_ - x_NQ
         y_box_min = y_min_ - y_NQ
         y_box_max = y_max_ - y_NQ
         z_box_min = z_min_ - z_NQ
         z_box_max = z_max_ - z_NQ
*
         Roots(1,1)=x_box_min
         Roots(2,1)=x_box_max
         If (x_box_max*x_box_min.lt.Zero) Then
            nx_Roots=3
            Roots(3,1)=Zero
         Else
            nx_Roots=2
         End If
*
         Roots(1,2)=y_box_min
         Roots(2,2)=y_box_max
         If (y_box_max*y_box_min.lt.Zero) Then
            ny_Roots=3
            Roots(3,2)=Zero
         Else
            ny_Roots=2
         End If
*
         Roots(1,3)=z_box_min
         Roots(2,3)=z_box_max
         If (z_box_max*z_box_min.lt.Zero) Then
            nz_Roots=3
            Roots(3,3)=Zero
         Else
            nz_Roots=2
         End If
*
*        Check all stationary points
*
         Do ix = 1, nx_Roots
            x = Roots(ix,1)
            Do iy = 1, ny_Roots
               y = Roots(iy,2)
               Do iz = 1, nz_Roots
                  z = Roots(iz,3)
*
                  r=Sqrt(x**2+y**2+z**2)
*
                  R_Box_Max=Max(R_Box_Max,r)
                  R_Box_Min=Min(R_Box_Min,r)
*
              End Do
           End Do
        End Do
*
        If (Abs(R_Box_Min).lt.1.0D-12) R_Box_Min=Zero
        R_Box_Max=R_Box_Max+1.0D-15
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
         If (Debug) write(6,*) 'Get_Subblock ----> Subblock'
#endif
c
c        Note that in gradient calculations we process the grid points
c        for each atomic grid seperately in order to used the
c        translational invariance on the atomic contributions to the
c        gradient.
c
         nTotGP_Save = nTotGP
         ip_iA=ip_of_iWork_d(Work(ip_Angular(iNQ)))
         ip_A=iWork(ip_iA)
         nR_Eff=iWork(ip_nR_eff-1+iNQ)
         Call Subblock(iNQ,x_NQ,y_NQ,z_NQ,InBox(iNQ),
     &                 x_min_,x_max_, y_min_,y_max_, z_min_,z_max_,
     &                 list_p,nlist_p,Grid,Weights,mGrid,.True.,
     &                 number_of_grid_points,R_Box_Min,R_Box_Max,
     &                 iList_p,xyz0,iWork(ip_A),nR_Eff)
         nTotGP = nTotGP_Save
*
#ifdef _DEBUGPRINT_
         If (Debug) write(6,*) 'Subblock ----> Get_Subblock'
#endif
      End Do
      iWork(ip_GridInfo+(ixyz-1)*2)=iDisk_Grid
      iWork(ip_GridInfo+(ixyz-1)*2+1)=nBatch
      Call iDaFile(Lu_Grid,1,iBatchInfo,3*nBatch,iDisk_Grid)
*                                                                      *
************************************************************************
*                                                                      *
*     Process grid points on file
*
 997  Continue
*
      iDisk_Grid =iGridInfo(1,ixyz)
      nBatch=iGridInfo(2,ixyz)
      Call iDaFile(Lu_Grid,2,iBatchInfo,3*nBatch,iDisk_Grid)
*
      iBatch = 0
      nogp=0
 888     Continue
         iBatch = iBatch + 1
         If (iBatch.gt.nBatch) Go To 996
 887     Continue
         jDisk_Grid=           iBatchInfo(1,iBatch)
         number_of_grid_points=iBatchInfo(2,iBatch)
*
         iNQ=                  iBatchInfo(3,iBatch)
#ifdef _DEBUGPRINT_
         If (Debug) Then
            Write (6,*)
            Write (6,*)  'iNQ=',iNQ
            Write (6,*)
         End If
#endif
         ilist_p=-1
         Do klist_p = 1, nlist_p
            If (List_p(klist_p).eq.iNQ) ilist_p=klist_p
         End Do
*
         If (nogp+number_of_grid_points.le.mGrid) Then
            Call dDaFile(Lu_Grid,2,Grid(1,nogp+1),
     &                   3*number_of_grid_points,jDisk_Grid)
            Call dDaFile(Lu_Grid,2,Weights(nogp+1),
     &                   number_of_grid_points,jDisk_Grid)
            nogp = nogp+number_of_grid_points
*
*           If this is not a gradient evaluation read next buffer if the
*           current one is not the last one.
*
            More_to_Come=.False.
            If (.Not.Do_Grad.and.iBatch.ne.nBatch) Go To 888
         Else
            More_to_Come=.True.
         End If
*
*        Here if it is a gradient evaluation or we have a buffer to
*        process.
*
         If (Do_Grad) Then
            Call ICopy(nGrad_Eff,[On],0,iTab(2,1),4)
            If (Grid_Type.eq.Moving_Grid) Then
               Do iGrad = 1, nGrad_Eff
                  jNQ = iTab(3,iGrad)
                  If (iNQ.eq.jNQ) iTab(2,iGrad)=Off
               End Do
*
*------------- Generate derivative with respect to the weights
*              if needed.
*
               Call dWdR(Grid,ilist_p,Weights,list_p,nlist_p,
     &                   Work(ip_dW_dR),nGrad_Eff,iTab,Work(ip_dW_Temp),
     &                   Work(ip_dPB),number_of_grid_points)
            End If
         End If
*
         Call mma_Allocate(Grid_AO,kAO,nogp,nBfn,nD,Label='Grid_AO')
         Call mma_Allocate(TabAO,mAO,nogp,nBfn,Label='TabAO')
         If (Do_Grad) Call mma_Allocate(TabAO_Short,kAO,nogp,nBfn,
     &                                  Label='TabAO_Short')
         TabAO_Pack(1:mAO*nogp*nBfn) => TabAO(:,:,:)
         If (Do_Grad) Then
           Call mma_allocate(dRho_dR,mdRho_dR,nogp,
     &                                  nGrad_eff,Label='dRho_dR')
         Else
           Call mma_allocate(dRho_dR,1,1,1,Label='dRho_dR')
         End If

         Call Do_Batch(Kernel,Func,nogp,
     &                 list_s,nlist_s,List_Exp,List_Bas,
     &                 iWork(ipIndex),nIndex,
     &                 FckInt,nFckDim,nFckInt,
     &                 iWork(ipTabAO),mAO,nSym,nD,nP2_ontop,
     &                 Do_Mo,
     &                 Work(ipTabMO),Work(ipTabSO),nMOs,
     &                 P2unzip,D1mo,D1Unzip,nd1mo,P2_ontop,
     &                 Do_Grad,Grad,nGrad,
     &                 mdRho_dR,nGrad_Eff,
     &                 list_g,IndGrd,iTab,Temp,
     &                 Work(ip_dW_dR),iNQ,
     &                 LTEG_DB,PDFTPot1,PDFTFocI,PDFTFocA)
*
         If (Allocated(dRho_dR)) Call mma_deallocate(dRho_dR)
         If (Allocated(TabAO_Short)) Call mma_deallocate(TabAO_Short)
         TabAO_Pack => Null()
         Call mma_deallocate(TabAO)
         Call mma_deallocate(Grid_AO)

         nTotGP=nTotGP+nogp
*        update the "LuGridFile":
         do i=1,nogp
         write(LuGridFile,'(3ES24.14,1x,ES24.14)')
     &               (Grid(l,i),l=1,3), Weights(i)
         enddo
         nogp=0
         If (More_To_Come) Go To 887
         Go To 888
 996     Continue
*
*                                                                      *
************************************************************************
*                                                                      *
 998  Continue
*
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deAllocate(iBfn_Index)
      Call mma_deAllocate(Ind_Grd)
      Call GetMem('Index','Free','Real',ipIndex,nIndex)
      Call Free_iWork(ipTabAO)
      Call mma_deallocate(Dens_AO)
      If (ipTabMO.ne.ip_Dummy) Call Free_Work(ipTabMO)
      If (ipTabSO.ne.ip_Dummy) Call Free_Work(ipTabSO)
      If (Do_Grad.and.Grid_Type.eq.Moving_Grid) Then
         Call GetMem('dPB','Free','Real',ip_dPB,ndPB)
         Call GetMem('dW_Temp','Free','Real',ip_dW_Temp,ndW_Temp)
         Call GetMem('dW_dR','Free','Real',ip_dW_dR,ndW_dR)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      End
