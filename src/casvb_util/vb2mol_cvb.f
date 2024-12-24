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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine vb2mol_cvb(vecvb,vecmol,isyml)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension vecvb(ndet),vecmol(*)

      iwr=0
      call icomb_cvb(norb,nalf,nsa)
      call icomb_cvb(norb,nbet,nsb)
      i1 = mstacki_cvb(nsa)
      i2 = mstacki_cvb(nsb)
      i3 = mstacki_cvb(mxirrep)
      i4 = mstacki_cvb(mxirrep)
      call mol2vb2_cvb(vecvb,vecmol,isyml,0d0,iwr,
     >  iwork(i1),iwork(i2),iwork(i3),iwork(i4),nsa,nsb)
      ibasemx=max(ibasemx,mstackr_cvb(0))
      call mfreei_cvb(i1)
      return
      end
      subroutine mol2vb_cvb(vecvb,vecmol,isyml)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension vecvb(ndet),vecmol(*)

      iwr=1
      call icomb_cvb(norb,nalf,nsa)
      call icomb_cvb(norb,nbet,nsb)
      i1 = mstacki_cvb(nsa)
      i2 = mstacki_cvb(nsb)
      i3 = mstacki_cvb(mxirrep)
      i4 = mstacki_cvb(mxirrep)
      call mol2vb2_cvb(vecvb,vecmol,isyml,0d0,iwr,
     >  iwork(i1),iwork(i2),iwork(i3),iwork(i4),nsa,nsb)
      ibasemx=max(ibasemx,mstackr_cvb(0))
      call mfreei_cvb(i1)
      return
      end

      subroutine mol2vbma_cvb(vecvb,vecmol,isyml,fac)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension vecvb(ndet),vecmol(*)
      iwr=2
      call icomb_cvb(norb,nalf,nsa)
      call icomb_cvb(norb,nbet,nsb)
      i1 = mstacki_cvb(nsa)
      i2 = mstacki_cvb(nsb)
      i3 = mstacki_cvb(mxirrep)
      i4 = mstacki_cvb(mxirrep)
      call mol2vb2_cvb(vecvb,vecmol,isyml,fac,iwr,
     >  iwork(i1),iwork(i2),iwork(i3),iwork(i4),nsa,nsb)
      ibasemx=max(ibasemx,mstackr_cvb(0))
      call mfreei_cvb(i1)
      return
      end
