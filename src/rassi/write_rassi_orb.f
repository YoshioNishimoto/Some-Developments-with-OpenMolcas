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
* Copyright (C) 2020, Bruno Tenorio                                    *
************************************************************************
<<<<<<< HEAD
!     Subroutine to correctly bunch together CMO matrices orbitals
!     and pass them to the molden_dysorb interface for .molden export
=======

!     Subroutine to correctly bunch together CMO matrices orbitals
!     and pass them to the molden_dysorb interface for .molden export

>>>>>>> 0893210f7 (New Cholesky orthonormalization of CMO2)
      SUBROUTINE WRITE_RASSI_ORB(NBASFN,NOSHN,NCMOA,CMOA,NUM,ist,jst)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "Molcas.fh"
#include "cntrl.fh"
#include "WrkSpc.fh"
#include "rassi.fh"
#include "prgm.fh"
#include "symmul.fh"
#include "Files.fh"

      INTEGER   NBASFN,NCMOA,NOSHN
<<<<<<< HEAD
      INTEGER   DYSCIND,NBSFA
      INTEGER coeff1 , coeff2
      DIMENSION DYSAMPS(NOSHN,NOSHN)
=======
      INTEGER   NBSFA
      INTEGER coeff1 , coeff2
>>>>>>> 0893210f7 (New Cholesky orthonormalization of CMO2)
      DIMENSION DYSEN(NOSHN)
      DIMENSION AMPS(NBASFN)
      DIMENSION CMOA(NCMOA)
      real*8, Allocatable:: CMOA1(:)
      INTEGER IBIO,IZZ,SYMOFF,BIOOFF,IBIOFF

      Character*30 Filename
      Character*3 NUM1 , NUM2
<<<<<<< HEAD
      Character*30 FNW
=======
>>>>>>> 0893210f7 (New Cholesky orthonormalization of CMO2)
      INTEGER NUM,ist,jst
      WRITE(NUM1,'(I3.3)') ist
      WRITE(NUM2,'(I3.3)') jst

      NBSFA=NBASFN*NOSHN
      ALLOCATE(CMOA1(NBSFA))
      DYSEN=0.0D0 ! Orbital energies
      AMPS=0.0D0 ! Transition amplitudes (shown as occupations)
      CMOA1=0.0D0
<<<<<<< HEAD
=======

>>>>>>> 0893210f7 (New Cholesky orthonormalization of CMO2)
      SYMOFF=0
      IBIOFF=0
      IZZOFF=0
      IBSA=0
      IBSD=0
      DO ISY1=1,NSYM
        NO1=NOSH(ISY1)
        NA1=NASH(ISY1)
        NB1=NBASF(ISY1)
<<<<<<< HEAD
=======

>>>>>>> 0893210f7 (New Cholesky orthonormalization of CMO2)
        IF(NA1.GT.0) THEN
        DO IBIO=1,NO1
         DO IZZ=1,NB1
          BIOOFF=(IBIO-1)*NB1
          coeff1=IZZ+IZZOFF+IBSA
          coeff2=SYMOFF+BIOOFF+IZZ
          CMOA1(coeff1)=CMOA(coeff2)
         END DO
         IBSA=IBSA+NBASFN
        END DO
<<<<<<< HEAD
        END IF
=======

        END IF

>>>>>>> 0893210f7 (New Cholesky orthonormalization of CMO2)
        IZZOFF=NB1+IZZOFF
        IBIOFF=NO1+IBIOFF
        SYMOFF=(NO1*NB1)+SYMOFF
      END DO

! If at least one orbital was found, export it/them
      IF(NUM.Eq.2) THEN
       !WRITE(FNM,'(A25)') 'rassi.cmo2.molden.'//NUM1//'_'//NUM2
       !FNM='rassi.cmo2.molden.'//NUM1//'_'//NUM2
      Write(filename,'(A18,A3,A1,A3)')'rassi.cmo2.molden.',NUM1,'_',NUM2
      Call Molden_DysOrb(filename,DYSEN,AMPS,CMOA1,NOSHN,NBASFN)
      END IF
      IF(NUM.Eq.1) THEN
       !WRITE(FNM,'(A25)') 'rassi.cmo1.molden.'//NUM1//'_'//NUM2
       !FNM='rassi.cmo1.molden.'//NUM1//'_'//NUM2
      Write(filename,'(A18,A3,A1,A3)')'rassi.cmo1.molden.',NUM1,'_',NUM2
      Call Molden_DysOrb(filename,DYSEN,AMPS,CMOA1,NOSHN,NBASFN)
      END IF

      RETURN
      END
