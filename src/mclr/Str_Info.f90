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
! Copyright (C) 2020, Roland Lindh                                     *
!***********************************************************************
!        OCSTR        :        Offsets for occupation of strings
!        STREO        :        reordering array
!        STSM         :        Symmetry of each string
!        STCL         :        Class of each string
!        NSTSO        :        Number of strings per symmetry and occupation
!        ISTSO        :        Offset of strings per symmetry and occupation
!        EL1          :        Number of electrons in RAS1 per sub type
!        EL3          :        Number of electrons in RAS3 per sub type
!        ACTP         :        is sub-type active
!        Z            :         Lexical adressing of arrays
!        EL123        :        -"-    But array
!        STSTMI       :        Explicit offsets and lengths
!        STSTMN       :                  "
!        STSTM        :        ?
!        NDMAP        :        Down mappings of strings containing the same number of electrons
!        NUMAP        :          Up mappings of strings containing the same number of electrons


! Not used
!        COBSM        :        Symmetry of conjugated orbitals
!        NIFSJ        :
!        IFSJ         :
!        IFSJO        :
!        STSTX        :        Symmetry of excitation connecting strings of given symmetry

Module Str_Info

Implicit None
Private
Public::String_Info, Str, NSTTYP_MAX, ITYP_DUMMY, NSTTYP, NELEC, MNRS1, MXRS1, MNRS3, MXRS3, IZORR, ISTTP, iuniqmp, iuniqtp, &
        IAZTP,IBZTP,IARTP,IBRTP, NZSTTP,NRSTTP, IATPM1,IATPM2,IBTPM1,IBTPM2, ISTAC, NOCTYP, NSTFTP, INUMAP, INDMAP, MXNSTR,  &
        DFTP, CFTP, DTOC, Storage, CNSM


Type String_Info
     Sequence
     Integer, Pointer:: OCSTR(:)=>Null()
     Integer, Allocatable:: OCSTR_hidden(:)
     Integer, Pointer:: STREO(:)=>Null()
     Integer, Allocatable:: STREO_hidden(:)
     Integer, Pointer:: STSM(:)=>Null()
     Integer, Allocatable:: STSM_hidden(:)
     Integer, Pointer:: STCL(:)=>Null()
     Integer, Allocatable:: STCL_hidden(:)
     Integer, Pointer:: NSTSO(:)=>Null()
     Integer, Allocatable:: NSTSO_hidden(:)
     Integer, Pointer:: ISTSO(:)=>Null()
     Integer, Allocatable:: ISTSO_hidden(:)
     Integer, Pointer:: EL1(:)=>Null()
     Integer, Allocatable:: EL1_hidden(:)
     Integer, Pointer:: EL3(:)=>Null()
     Integer, Allocatable:: EL3_hidden(:)
     Integer, Pointer:: ACTP(:)=>Null()
     Integer, Allocatable:: ACTP_hidden(:)
     Integer, Pointer:: Z(:)=>Null()
     Integer, Allocatable:: Z_hidden(:)
     Integer, Pointer:: EL123(:)=>Null()
     Integer, Allocatable:: EL123_hidden(:)
     Integer, Allocatable:: STSTMI(:)
     Integer, Allocatable:: STSTMN(:)
     Integer, Pointer:: STSTM(:,:)=>Null()
     Integer, Allocatable:: STSTM_hidden(:,:)
     Integer, Allocatable:: NUMAP(:)
     Integer, Allocatable:: NDMAP(:)
End Type String_Info

Type (String_Info), Allocatable, Target:: Str(:)
!     Integer, Allocatable:: COBSM(:)
!     Integer, Allocatable:: NIFSJ(:)
!     Integer, Allocatable:: IFSJ(:)
!     Integer, Allocatable:: IFSJO(:)
!     Integer, Allocatable:: STSTX(:)

!             INITITIALIZED IN STRTYP
!        NSTTYP       :        Number of string types
!        MNRS1        :        Min ras1
!        MXRS1        :        Max ras1
!        MNRS3        :        Min ras3
!        MXRS3        :        Max ras3
!        NELEC        :        Number of electrons
!        IZORR        :        Zero orde y/n
!        IAZTP        :        Pointer to alpha types
!        IBZTP        :        Pointer to beta types
!
!not in use (just zero order space)
!        IARTP        :        Give type nr to certain exc
!        IBRTP        :        Give type nr to certain exc
!
!        NZSTTP       :        Not in use
!        NRSTTP       :        Not in use
!
!        ISTTP        :        Space (0=zero order)
!        iuniqmp      :        Unique types (not necessary here just 0order space)
Integer, Parameter:: NSTTYP_MAX=6+1   ! "+1" is the dummy layer
Integer:: ITYP_DUMMY=0
Integer i
Integer::     NSTTYP,                                         &
              NELEC(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],         &
              MNRS1(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],         &
              MXRS1(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],         &
              MNRS3(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],         &
              MXRS3(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],         &
              IZORR(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],         &
              ISTTP(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],         &
              iuniqmp(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],       &
              iuniqtp(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],       &
              IAZTP,IBZTP,IARTP(3,10),IBRTP(3,10),            &
              NZSTTP,NRSTTP,                                  &
              IATPM1,IATPM2,IBTPM1,IBTPM2

!        ISTAC                 : Stringtype maping; a(or a+) i -> istac(j,1(2))
!        NOCTYP                : Number of occupation classes for given type
!        NSTFTP                : Number of strings of this type
!        INUMAP                : Mapping of string type to next more general type
!        INDMAP                : Mapping of string type to next more restricted type
!        MXNSTR                : Largest number of strings of given sym and type
!
Integer:: ISTAC(NSTTYP_MAX,2),                          &
          NOCTYP(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],      &
          NSTFTP(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],      &
          INUMAP(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],      &
          INDMAP(NSTTYP_MAX)=[(0,i=1,NSTTYP_MAX)],      &
          MXNSTR

#include "detdim.fh"
!     DFTP          :        OPEN SHELL DETERMINANTS OF PROTO TYPE
!     CFTP          :        BRANCHING DIAGRAMS FOR PROTO TYPES
!     DTOC          :        CSF-DET TRANSFORMATION FOR PROTO TYPES
!     CNSM(:)%ICONF :        NCNSM  CONFIGURATION EXPANSIONS
!     CNSM(I)%ICTS  :        adress of determinant I in STRING ordering for
!                            determinant I in CSF ordering
!                            reference symmetry IREFSM.
Integer, Allocatable:: DFTP(:)
Integer, Allocatable:: CFTP(:)
Real*8,  Allocatable:: DTOC(:)
Type Storage
  Integer, Allocatable:: ICONF(:)
  Integer, Allocatable:: ICTS(:)
End Type Storage
Type (Storage) :: CNSM(MXCNSM)

End Module Str_Info
