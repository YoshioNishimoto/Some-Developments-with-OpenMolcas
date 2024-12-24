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
c  *************************************************************
c  ** Routines for imposing constraints on VB wfn. parameters **
c  *************************************************************
c  *********************
c  ** Set-up routines **
c  *********************
      subroutine construc_cvb(tconstr,ipermzeta)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"

#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension tconstr(nvb,nvb),ipermzeta(norb,nzeta)

      i1 = mstackr_cvb(norb*norb)
      i2 = mstackr_cvb(norb*norb)
      i3 = mstackr_cvb(norb*norb)
      call setipermzeta_cvb(ipermzeta,
     >  work(lv(1)),work(ls(1)),iwork(ls(13)),
     >  work(i1),work(i2),work(i3))
      call mfreer_cvb(i1)
      if(iconstruc.eq.2)call construc2_cvb(tconstr)
      return
      end
