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
      subroutine axexbsol_cvb(ap,rhsp,itdav,maxdav,nfrdim,
     >  solp,solp_res,eig,eig_res)
      implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
      dimension ap(maxdav,maxdav),rhsp(maxdav)
      dimension solp(maxdav),solp_res(maxdav)

      i1 = mstackr_cvb(itdav)
      i2 = mstackr_cvb(itdav*itdav)
      i3 = mstackr_cvb(itdav)
      i4 = mstackr_cvb(itdav)
      i5 = mstackr_cvb(itdav)
      call axexbsol2_cvb(ap,rhsp,itdav,maxdav,nfrdim,
     >  solp,solp_res,eig,eig_res,
     >  work(i1),work(i2),work(i3),work(i4),work(i5))
      call mfreer_cvb(i1)
      return
      end
