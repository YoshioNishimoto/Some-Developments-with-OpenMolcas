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
      SUBROUTINE SMMAT2(PROP,PRMAT,NSS,PRLBL,IPRCMP,IJ)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) PRLBL
      DIMENSION PRMAT(NSS,NSS)
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SMMAT')
      Integer IJ(4)
#include "SysDef.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
      DIMENSION PROP(NSTATE,NSTATE,NPROP)
      REAL*8, EXTERNAL :: DCLEBS
*
      IPRNUM=-1
C IFSPIN takes values the values 0,1,2
C 0 = spin free property
C 1 = spin operator (S)
C 2 = spin dependent property, triplet operator
      IFSPIN=0

      DO IPROP=1,NPROP
         IF (PRLBL.EQ.PNAME(IPROP)) THEN
            IF (PRLBL(1:5).eq.'TMOS0') THEN
               IFSPIN=2
               IPRNUM=IPROP
               EXIT
            ELSE
               IFSPIN=0
               IF (IPRCMP.EQ.ICOMP(IPROP)) THEN
                  IPRNUM=IPROP
                  EXIT
               END IF
            END IF
         ELSE IF (PRLBL(1:4).EQ.'SPIN') THEN
            IFSPIN=1
            IPRNUM=0
            EXIT
         END IF
      END DO
      IF (IPRNUM.EQ.-1) THEN
         Write (6,*) 'SMMAT_CHECK, Abend IPRNUM.EQ.-1'
         Write (6,*) 'SMMAT_CHECK, PRLBL=','>',PRLBL,'<'
         Call Abend()
      ENDIF
C Mapping from spin states to spin-free state and to spin:
      ISS=0
      DO ISTATE=1,IJ(3)-1
        JOB1=iWork(lJBNUM+ISTATE-1)
        MPLET1=MLTPLT(JOB1)
        ISS=ISS+MPLET1
      End Do
      DO ISTATE=IJ(3),IJ(4)
       JOB1=iWork(lJBNUM+ISTATE-1)
       MPLET1=MLTPLT(JOB1)
       S1=0.5D0*DBLE(MPLET1-1)

       DO MSPROJ1=-MPLET1+1,MPLET1-1,2
        SM1=0.5D0*DBLE(MSPROJ1)
        ISS=ISS+1

        JSS=0
        DO JSTATE=1,IJ(1)-1
         JOB2=iWork(lJBNUM+JSTATE-1)
         MPLET2=MLTPLT(JOB2)
         JSS=JSS+MPLET2
        End Do
        DO JSTATE=IJ(1),IJ(2)
         JOB2=iWork(lJBNUM+JSTATE-1)
         MPLET2=MLTPLT(JOB2)
         S2=0.5D0*DBLE(MPLET2-1)

         DO MSPROJ2=-MPLET2+1,MPLET2-1,2
          SM2=0.5D0*DBLE(MSPROJ2)
          JSS=JSS+1

          IF (IFSPIN.EQ.0 .AND. IPRNUM.NE.0) THEN
                  IF (MPLET1.EQ.MPLET2 .AND. MSPROJ1.EQ.MSPROJ2) THEN
                          PRMAT(ISS,JSS)=PROP(ISTATE,JSTATE,IPRNUM)
                  END IF
          ELSE IF (IFSPIN.EQ.1 .AND. IPRNUM.EQ.0) THEN
                  SXMER=0.0D0
                  SYMEI=0.0D0
                  SZMER=0.0D0
                  SMINUS=0.0D0
                  SPLUS=0.0D0
                  IF((ISTATE.EQ.JSTATE).AND.(MPLET1.eq.MPLET2)) THEN
                    IF (MSPROJ1.eq.MSPROJ2-2) THEN
                      SMINUS=SQRT((S1+SM2)*(S1-SM1))
                      SXMER= 0.5D0*SMINUS
                      SYMEI= 0.5D0*SMINUS
                    ELSE IF (MSPROJ1.eq.MSPROJ2) THEN
                      SZMER=SM1
                    ELSE IF (MSPROJ1.eq.MSPROJ2+2) THEN
                      SPLUS=SQRT((S1+SM1)*(S1-SM2))
                      SXMER= 0.5D0*SPLUS
                      SYMEI=-0.5D0*SPLUS
                    END IF
                    IF (IPRCMP.EQ.1) THEN
                            PRMAT(ISS,JSS)=SXMER
                    ELSE IF (IPRCMP.EQ.2) THEN
                            PRMAT(ISS,JSS)=SYMEI
                    ELSE IF (IPRCMP.EQ.3) THEN
                            PRMAT(ISS,JSS)=SZMER
                    END IF
                  END IF
          ELSE IF (IFSPIN.EQ.2) THEN
C 1-electron triplet operator so only Delta S =0,+-1 and Delta MS =0,+-1
C Notice S1,S2 and SM1,SM2 need not to be integers
C Hence MPLET1,MPLET2 and MSPROJ1,MSPROJ2 are used
C Notice SMINUS and SPLUS is interchanged compared to above
C Notice that the Y part is imaginary
C
C see section 3 (Spin_orbit coupling in RASSI) in
C P A Malmqvist, et. al CPL, 357 (2002) 230-240
C for details
C
C Note that we work on the x, y, and z components at this time.
C
C On page 234 we have the notation V^{AB}(x), that is the
C potential has Cartesian components. Here, however, this is
C partitioned in a slightly different way since we have that
C V^{AB}(x)=(k x e_l)_x V^{AB}. We will only handle the
C V^{AB} part.
C
C Hence, we will compute the contributions to T(i), i=x,y,z
C here and form the inner product
C (k x e_l)_i V^{AB} . T(i)
C outside the code.
C
C
C Note that this code stricktly follows the code of soeig.f where
C the term L.S is added to the Hamiltonian.
C
                  FACT=1.0D0/SQRT(DBLE(MPLET1))
                  IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
                  CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
                  CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
                  CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
                  CGX= SQRT(0.5D0)*(CGM-CGP)
                  CGY=-SQRT(0.5D0)*(CGM+CGP)
*
                  EXPKR=PROP(ISTATE,JSTATE,IPRNUM)
*
                  IF (IPRCMP.EQ.1) THEN
                    EXPKR=EXPKR*CGX
                  ELSE IF (IPRCMP.EQ.2) THEN
                    EXPKR=EXPKR*CGY
                  ELSE IF (IPRCMP.EQ.3) THEN
                    EXPKR=EXPKR*CG0
                  END IF
                  PRMAT(ISS,JSS)= EXPKR
          END IF

         END DO
        END DO
       END DO
      END DO

      RETURN
      END
