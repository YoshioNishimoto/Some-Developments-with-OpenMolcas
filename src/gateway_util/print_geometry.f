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
* Copyright (C) 2006, Roland Lindh                                     *
************************************************************************
      SubRoutine Print_Geometry(iOpt,DInf,nDInf)
************************************************************************
*                                                                      *
*     Object: to print the molecular coordinates, bonds, angles and    *
*             torsional angles.                                        *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
*             September '06                                            *
************************************************************************
      use Period
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "print.fh"
      Character help_c*1
      Character FMT*16
      Real*8 DInf(nDInf)
      Real*8, Dimension (:,:), Allocatable :: Centr
#include "angstr.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      iPrint = nPrint(iRout)
      If (iPrint.eq.0) Return
      Call qEnter('Print_Geometry')
      LuWr=6
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Centr,3,mCentr)
*                                                                      *
************************************************************************
*                                                                      *

      Write (LuWr,*)
      Call CollapseOutput(1,'   Molecular structure info:')
      Write (LuWr,'(3X,A)') '   -------------------------'
      Write (LuWr,*)
      if(iOpt.eq.0) FMT='(19X,A)'
      if(iOpt.eq.1) FMT='(11X,A)'
*
      Write (LuWr,FMT)
     &       ' ************************************************ '
      if(iOpt.eq.0) then
      Write (LuWr,FMT)
     &       ' **** Cartesian Coordinates / Bohr, Angstrom **** '
      else
      Write (LuWr,FMT)
     &       ' **** Cartesian Coordinates / Angstrom       **** '
            endif
      Write (LuWr,FMT)
     &       ' ************************************************ '
      Write (LuWr,*)
      if(iOpt.eq.0) then
      Write (LuWr,'(A,A,A)')  '     Center  Label ',
     &             '               x              y              z',
     &  '                     x              y              z'
      else
      Write (LuWr,'(A,A,A)')  '     Center  Label ',
     &             '               x              y              z'
      endif
*                                                                      *
************************************************************************
*                                                                      *
      ndc = 0
      nc = 1
      Do jCnttp = 1, nCnttp
         mCnt = nCntr(jCnttp)
         If (AuxCnttp(jCnttp).or.FragCnttp(jCnttp))Then
            ndc = ndc + mCnt
            Go To 32
         End If
         jxyz = ipCntr(jCnttp)
         Do jCnt = 1, mCnt
            ndc = ndc + 1
            x1 = Dinf(jxyz)
            y1 = Dinf(jxyz+1)
            z1 = Dinf(jxyz+2)
            Do i = 0, nIrrep/nStab(ndc) - 1
               Facx=DBLE(iPhase(1,iCoset(i,0,ndc)))
               Facy=DBLE(iPhase(2,iCoset(i,0,ndc)))
               Facz=DBLE(iPhase(3,iCoset(i,0,ndc)))
               If (Show) Then
                  help_c = ' '
                  If(Cell_l) Then
                     Do j=1,lthCell
                         If(AdCell(j).EQ.ndc) help_c = '*'
                     End Do
                  End If
                  if(iOpt.eq.0) then
                     Write (LuWr,
     &                   '(6X,I3,A1,5X,A,3F15.6,7X,3F15.6)')
     &                    nc, help_c, LblCnt(ndc),
     &                    x1*Facx, y1*Facy, z1*Facz,
     &                    x1*Facx*angstr,
     &                    y1*Facy*angstr,
     &                    z1*Facz*angstr
                  else
                     Write (LuWr,
     &                   '(6X,I3,A1,5X,A,3F15.6)')
     &                    nc, help_c, LblCnt(ndc),
     &                    x1*Facx*angstr,
     &                    y1*Facy*angstr,
     &                    z1*Facz*angstr
                  endif
               End If
               Centr(1,nc) = x1*Facx
               Centr(2,nc) = y1*Facy
               Centr(3,nc) = z1*Facz
               nchr=iAtmNr(jCnttp)
               nchr=iAtmNr(jCnttp)
               if (nc.gt.8*mxdc) Then
                  Call WarningMessage(2,'lblxxx too small')
                  Call Abend()
               End If
               lblxxx(nc)=lblcnt(ndc)(1:LENIN)
               nc = nc + 1
            End Do
            jxyz = jxyz + 3
         End Do
32       Continue
      End Do
      nc=nc-1
*                                                                      *
************************************************************************
*                                                                      *
*     Compute distances
*
      If (mCentr.le.2) Go To 55
      Call Dstncs(lblxxx,Centr,nc,
     &            angstr,Max_Center,6)
      If (.Not.Expert) Call DstChk(Centr,lblxxx,nc)
*
*     Compute valence bond angels
*
      If (iPrint.lt.5.or.mCentr.lt.3.or.iOpt.eq.1) Go To 55
      Call Angles(lblxxx,Centr,nc,rtrnc,Max_Center)
*
*     Compute dihedral angles
*
      If (iPrint.lt.5.or.mCentr.lt.4) Go To 55
      Call Dihedr(lblxxx,Centr,nc,rtrnc,Max_Center)
 55   Continue
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(Centr)
*                                                                      *
************************************************************************
*                                                                      *
      Call CollapseOutput(0,'   Molecular structure info:')
      Write (LuWr,*)
      Call qExit('Print_Geometry')
      Return
      End
