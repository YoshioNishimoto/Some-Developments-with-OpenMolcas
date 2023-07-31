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
      subroutine span1_cvb(c,nvec,s,n,metr)
      implicit real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "span_comcvb.fh"
      dimension c(n,nvec),s(*)

      nvremain=nvec
100   nvmove=min(nvremain,nvecmx-nvtot)
      if(nvmove.eq.0.and.nvremain.gt.0)then
        write(6,*)' Fatal error in SPAN_CVB!',nvmove,nvremain
        call abend_cvb()
      endif
      call fmove_cvb(c(1,1+nvec-nvremain),work(nvtot*n+iaddr),n*nvmove)
      nvtot=nvtot+nvmove
      if(nvtot.eq.nvecmx)call span_cvb(work(iaddr),nvtot,nvtot,s,n,metr)
      nvremain=nvremain-nvmove
      if(nvremain.gt.0)goto 100
      return
      end
