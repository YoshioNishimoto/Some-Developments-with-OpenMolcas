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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
c  *********************************************************************
c  *                                                                   *
c  *  DEV2B  := calculate two-electron Hessian                         *
c  *                                                                   *
c  *********************************************************************
      subroutine dev2b_cvb(v1,v2,cfrom,hessorb,hesst,
     >  oaa2,aa1,gx,grad2)
c  Calculates V1 EijEkl CFROM and V2 EijEkl CFROM
      implicit real*8 (a-h,o-z)
#include "ext_cvb.fh"
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "malloc_cvb.fh"
      dimension hessorb(nprorb,nprorb),hesst(norb*norb,norb*norb)
      dimension gx(norb,norb),grad2(nprorb),v1(*),v2(*),cfrom(*)

      iv1=nint(v1(1))
      iv2=nint(v2(1))
      icfrom=nint(cfrom(1))
      n_2el=n_2el+2
      if(iform_ci(icfrom).ne.0)then
        write(6,*)' Unsupported format in DEV2B :',iform_ci(icfrom)
        call abend_cvb()
      endif

      call dev2b_2_cvb(w(iaddr_ci(iv1)),w(iaddr_ci(iv2)),
     >  w(iaddr_ci(icfrom)),
     >  hessorb,hesst,oaa2,aa1,nprorb,
     >  iw(ll(1)),iw(ll(2)),iw(ll(3)),iw(ll(4)),iw(ll(5)),iw(ll(6)),
     >  w(ll(9)),w(ll(10)),
     >  iw(ll(11)),iw(ll(12)),iw(ll(13)),iw(ll(14)),npvb,
     >  gx,grad2,
     >  nda,ndb,n1a,n1b,nam1,nbm1,norb,projcas,sc,absym(3))
      return
      end
