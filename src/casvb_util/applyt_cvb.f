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
c  ************************************
c  ** Routines involving CI and ORBS **
c  ************************************
c  *********************************************************************
c  *                                                                   *
c  *  APPLYT := transform CI vector IVEC according to orbital          *
c  *         := transformation in IGJORB (from GAUSSJ)                 *
c  *                                                                   *
c  *  :> IVEC:      starting CI vector                                 *
c  *  :> IGJORB:    simple orbital updates                             *
c  *  :> N_APPLYT   stats                                              *
c  *  :> I1ALF:     string handling                                    *
c  *  :> I1BET:           do.                                          *
c  *  :> IATO:            do.                                          *
c  *  :> IBTO:            do.                                          *
c  *  :> PHATO:           do.                                          *
c  *  :> PHBTO:           do.                                          *
c  *                                                                   *
c  *  :< IVEC:      transformed CI vector                              *
c  *  :< N_APPLYT   stats                                              *
c  *                                                                   *
c  *********************************************************************
      subroutine iapplyt_cvb(cvec,igjorb)
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "WrkSpc.fh"
      dimension igjorb(*),cvec(*)

      call iapplyt_cvb_internal(igjorb)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine iapplyt_cvb_internal(igjorb)
      use iso_c_binding
      integer, target :: igjorb(*)
      real*8, pointer :: gjorb(:)
      ivec=nint(cvec(1))
      n_applyt=n_applyt+1
      ioff=idbl_cvb(norb*norb)
      if(iform_ci(ivec).eq.0)then
        call permci_cvb(cvec,igjorb(1+ioff))
        call c_f_pointer(c_loc(igjorb(1)),gjorb,[1])
        call applyt2_cvb(work(iaddr_ci(ivec)),
     >    gjorb,igjorb(1+norb+ioff),
     >    iwork(ll(1)),iwork(ll(2)),iwork(ll(5)),iwork(ll(6)),
     >    work(ll(9)),work(ll(10)))
        nullify(gjorb)
      else
        write(6,*)' Unsupported format in APPLYT :',iform_ci(ivec)
        call abend_cvb()
      endif

      call setcnt2_cvb(ivec,0)
      return
      end subroutine iapplyt_cvb_internal
*
      end
*
      subroutine applyt_cvb(cvec,gjorb)
      implicit real*8 (a-h,o-z)
      dimension gjorb(*),cvec(*)
*
      call applyt_cvb_internal(gjorb)
*
*     This is to allow type punning without an explicit interface
      contains
      subroutine applyt_cvb_internal(gjorb)
      use iso_c_binding
      real*8, target :: gjorb(*)
      integer, pointer :: igjorb(:)
      call c_f_pointer(c_loc(gjorb(1)),igjorb,[1])
      call iapplyt_cvb(cvec,igjorb)
      nullify(igjorb)
      end subroutine applyt_cvb_internal
*
      end
