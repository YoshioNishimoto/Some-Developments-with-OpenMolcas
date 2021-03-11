************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE IJKL(INTSYM,INDX,C,S,FIJKL)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),INDX(*),C(*),S(*),
     *          FIJKL(*)
*
      JSYM(L)=JSUNP(INTSYM,L)
*------
* POW: Unnecessary but warning stopping initialization
      fini=1.0d30
*------
      ICHK=0
      NIJ=IROW(LN+1)
      NIJKL=NIJ*(NIJ+1)/2
      CALL FZERO(FIJKL,NIJKL)
      IADR=LASTAD(1)
201   CALL dDAFILE(Lu_70,2,VALSRT,NSRTMX,IADR)
      CALL iDAFILE(Lu_70,2,INDSRT,NSRTMX+2,IADR)
      LENGTH=INDSRT(NSRTMX+1)
      IADR=INDSRT(NSRTMX+2)
*       IF(LENGTH.GT.0) CALL SCATTER(LENGTH,FIJKL,INDSRT,VALSRT)
      do i=1,length
        FIJKL(INDSRT(i))=VALSRT(i)
      end do
      IF(IADR.NE.-1) GO TO 201
      IADD10=IAD10(5)
100   CALL dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
      CALL iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
      LEN=ICOP1(nCOP+1)
      IF(LEN.EQ.0)GO TO 100
      IF(LEN.LT.0)GO TO 200
      DO 10 IN=1,LEN
      IND=ICOP1(IN)
      IF(ICHK.NE.0) THEN
        ICHK=0
        INDI=IND
*        IP=MOD(INDI,2**8)
*        JP=MOD(INDI/2**8,2**8)
*        KP=MOD(INDI/2**16,2**8)
*        LP=MOD(INDI/2**24,2**8)
      IP=IBITS(INDI, 0,8)
      JP=IBITS(INDI, 8,8)
      KP=IBITS(INDI,16,8)
      LP=IBITS(INDI,24,8)
        NIJ=IROW(IP)+JP
        NKL=IROW(KP)+LP
        IND=NIJ*(NIJ-1)/2+NKL
        FINI=FIJKL(IND)
        GOTO 10
      END IF
      IF(IND.EQ.0) THEN
        ICHK=1
        GOTO 10
      END IF
*      IVL=MOD(IND,2**6)
*      IC2=MOD(IND/2**6,2**13)
*      IC1=MOD(IND/2**19,2**13)
      IVL=IBITS(IND, 0, 6)
      IC2=IBITS(IND, 6,13)
      IC1=IBITS(IND,19,13)
      COPI=COP(IN)*FINI
      IF(IVL.EQ.0) THEN
        S(IC1)=S(IC1)+COPI*C(IC2)
        S(IC2)=S(IC2)+COPI*C(IC1)
        GO TO 10
      END IF
      INDA=IRC(IVL)+IC1
      INDB=IRC(IVL)+IC2
      NA=INDX(INDA)
      NB=INDX(INDB)
      NS1=JSYM(INDA)
      NS1L=MUL(NS1,LSYM)
      INUM=NVIR(NS1L)
      IF(IVL.GE.2)INUM=NVPAIR(NS1L)
      CALL DAXPY_(INUM,COPI,C(NB+1),1,S(NA+1),1)
      CALL DAXPY_(INUM,COPI,C(NA+1),1,S(NB+1),1)
10    CONTINUE
      GO TO 100
200   CONTINUE
      RETURN
      END
