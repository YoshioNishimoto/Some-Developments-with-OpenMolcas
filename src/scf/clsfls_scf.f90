!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************
      Subroutine ClsFls_SCF()
!***********************************************************************
!                                                                      *
!     purpose: Close files after SCF calculations                      *
!                                                                      *
!***********************************************************************
!
#ifdef _HDF5_
      Use mh5, Only: mh5_close_file
#endif
      use InfSCF
      use Files
      use SCFWfn
      Implicit Real*8 (a-h,o-z)
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!
!---  close two-electron integral file --------------------------------*
      If (.Not.DSCF .And. .Not.DoCholesky) Then
         iRc=-1
         Call ClsOrd(iRc)
         If (iRc.ne.0) Then
            Write (6,*) 'ClsFls: Error closing ORDINT'
            Call Abend()
         End If
      End If
!
!---  close DNSMAT, dVxcdR, TWOHAM and GRADIENT -----------------------*
      Call DaClos(LuDSt)
      Call DaClos(LuOSt)
      Call DaClos(LuTSt)
      Call DaClos(LuGrd)
!
!---  close 2nd order updatefiles       -------------------------------*
      Call DaClos(LuDGd)
      Call DaClos(Lux)
      Call DaClos(LuDel)
      Call DaClos(Luy)
!
#ifdef _HDF5_
      call mh5_close_file(wfn_fileid)
#endif

      End Subroutine ClsFls_SCF