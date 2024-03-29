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
      subroutine ncset_cvb(ic)
      implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "loopcntr_cvb.fh"
      external istkprobe_cvb
      logical istkprobe_cvb

      if(istkprobe_cvb(istackrep))then
        call istkpop_cvb(istackrep,nc_zeroed)
        call istkpop_cvb(istackrep,nconvinone)
        if(ic.eq.0.or.ic.eq.1)then
          nconvinone=nconvinone+1
        elseif(ic.gt.1)then
          nconvinone=0
          nc_zeroed=1
        else
          nconvinone=-1
          nc_zeroed=1
        endif
        call istkpush_cvb(istackrep,nconvinone)
        call istkpush_cvb(istackrep,nc_zeroed)
      endif
      return
      end
