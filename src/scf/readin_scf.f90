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
!               1995, Martin Schuetz                                   *
!***********************************************************************
      SubRoutine ReadIn_SCF(SIntTh)
!***********************************************************************
!                                                                      *
!     purpose: Read input to SCF: one-electron integrals, informations *
!              about basis set and molecule and options to SCF.        *
!                                                                      *
!***********************************************************************
      use InfSCF, only: DSCF, nDisc, EThr, KSDFT, nCore, TimFld
      Implicit None
!
      Real*8 SIntTh

      Real*8 CPU1, CPU2, Tim1, Tim2, Tim3
      Logical PkMode
!
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
!                                                                      *
!***********************************************************************
!                                                                      *
!     read one electron integral file header
!
      Call R1IBas()
!                                                                      *
!***********************************************************************
!                                                                      *
!     read input
!
      Call RdInp_SCF()
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Allocate memory for SCF procedure
!
      Call MemAlo()
!                                                                      *
!***********************************************************************
!                                                                      *
!     read one-electron integrals
!
      Call R1IntA()
!                                                                      *
!***********************************************************************
!                                                                      *
!     Initialize seward
!
      Call IniSew_scf(DSCF,EThr,SIntTh,KSDFT)
!                                                                      *
!***********************************************************************
!                                                                      *
!     setup for direct or conventional integral calculations
!
      If (DSCF) Then
         Call Set_Basis_Mode('Valence')
         Call Setup_iSD()
         Call AlloK2()
         Call Free_iSD()
!
!------- Initiate integral packing for semi-direct implementation
!
         If (nDisc.ne.0) Then
            PkMode=.True.
            Call Ini_PkR8(PkMode)
         End If
!
!------- Allocate buffers for semi-direct SCF
!
         Call IniBuf(nDisc,nCore)
!
      Else
!
!------ Read the header of the two-electron integral file
!
        Call Rd2Int_SCF()
!
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 1) = TimFld( 1) + (Cpu2 - Cpu1)
      End subroutine ReadIn_SCF



      Subroutine Ini_PkR8(PkMode)
      use Gateway_Info, only: PkAcc
      Logical PkMode
!
      Call inipkr8(PkAcc,PkMode)
!
      End Subroutine Ini_PkR8
