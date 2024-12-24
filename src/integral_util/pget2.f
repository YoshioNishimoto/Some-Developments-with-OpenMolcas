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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine PGet2(iCmp,iBas,jBas,kBas,lBas,
     &                  Shijij, iAO, iAOst, nijkl,PSO,nPSO,
     &                  DSO,DSSO,nDSO,ExFac,CoulFac,PMax)
!***********************************************************************
!                                                                      *
!  Object: to assemble the 2nd order density matrix of a SCF wave      *
!          function from the 1st order density matrix.                 *
!                                                                      *
!          The indices has been scrambled before calling this routine. *
!          Hence we must take special care in order to regain the can- *
!          onical order.                                               *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             January '92.                                             *
!***********************************************************************
      use SOAO_Info, only: iAOtSO
      use Basis_Info, only: nBas
      use Symmetry_Info, only: nIrrep
      use Constants, only: Zero, Quart
#ifdef _DEBUGPRINT_
      use pso_stuff, only: iD0Lbl, D0
#endif
      Implicit None
      Integer iBas, jBas, kBas, lBas, nijkl, nPSO, nDSO
      Real*8 PSO(nijkl,nPSO), DSO(nDSO), DSSO(nDSO)
      Integer iCmp(4), iAO(4), iAOst(4)
      Logical Shijij
      Real*8 ExFac, CoulFac, PMax

!     Local Array
      Integer iSym(0:7), jSym(0:7), kSym(0:7), lSym(0:7)
      Integer, external:: iPntSO
      Integer lOper, MemSO2, i1, i2, i3, i4, j, niSym, njSym, nkSym,
     &        nlSym, iS, jS, kS, lS, j1, j2, j3, j123, j4, iSO, jSO,
     &        kSO, lSO, iSOi, jSOj, kSOk, lSOl, IndI, IndJ, IndK, IndL,
     &        j12, mijkl, iAOi, jAOj, kAOk, lAOl, ipntIJ, ipntKL,
     &        ipntIK, ipntIL, ipntJK, ipntJL, IndIJ, IndKL, IndIK,
     &        IndIL, IndJK, IndJL
      Real*8 t14, Temp
#ifdef _DEBUGPRINT_
      Integer iComp
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
      iComp = 1
      Call PrMtrx(' In PGet2:DSO ',[iD0Lbl],iComp,1,D0)
#endif
      lOper=1
      t14 = Quart * ExFac
      PMax=Zero
!
!-----Quadruple loop over elements of the basis functions angular
!     description.
!     Observe that we will walk through the memory in AOInt in a
!     sequential way.
!
      MemSO2 = 0
      Do 100 i1 = 1, iCmp(1)
         niSym = 0
         Do 101 j = 0, nIrrep-1
            If (iAOtSO(iAO(1)+i1,j)>0) Then
               iSym(niSym) = j
               niSym = niSym + 1
            End if
101      Continue
         Do 200 i2 = 1, iCmp(2)
            njSym = 0
            Do 201 j = 0, nIrrep-1
               If (iAOtSO(iAO(2)+i2,j)>0) Then
                  jSym(njSym) = j
                  njSym = njSym + 1
               End If
201         Continue
            Do 300 i3 = 1, iCmp(3)
               nkSym = 0
               Do 301 j = 0, nIrrep-1
                  If (iAOtSO(iAO(3)+i3,j)>0) Then
                     kSym(nkSym) = j
                     nkSym = nkSym + 1
                  End If
301            Continue
               Do 400 i4 = 1, iCmp(4)
                  nlSym = 0
                  Do 401 j = 0, nIrrep-1
                     If (iAOtSO(iAO(4)+i4,j)>0) Then
                        lSym(nlSym) = j
                        nlSym = nlSym + 1
                     End If
401               Continue
!
!------Loop over irreps which are spanned by the basis function.
!
       Do 110 is = 0, niSym-1
          j1 = iSym(is)
!
          Do 210 js = 0, njSym-1
             j2 = jSym(js)
             j12 = iEor(j1,j2)
!
             Do 310 ks = 0, nkSym-1
                j3 = kSym(ks)
                j123 = iEor(j12,j3)
                Do 410 ls = 0, nlSym-1
                   j4 = lSym(ls)
                   If (j123.ne.j4) Go To 410
!
                MemSO2 = MemSO2 + 1
!
!               Unfold the way the eight indices have been reordered.
                iSO = iAOtSO(iAO(1)+i1,j1)+iAOst(1)
                jSO = iAOtSO(iAO(2)+i2,j2)+iAOst(2)
                kSO = iAOtSO(iAO(3)+i3,j3)+iAOst(3)
                lSO = iAOtSO(iAO(4)+i4,j4)+iAOst(4)
!
                If (j1.ne.j2 .and. j1.ne.j3 .and. j1.ne.j4) Then
!------------------all irreps are different and the 2nd order density
!                  matrix will be identical to zero for a SCF type wave
!                  function.
                   call dcopy_(nijkl,[Zero],0,PSO(1,MemSO2),1)
                   Go To 310
                End If
!
                mijkl = 0
                Do 120 lAOl = 0, lBas-1
                   lSOl = lSO + lAOl
                   Do 220 kAOk = 0, kBas-1
                      kSOk = kSO + kAOk
                      Do 320 jAOj = 0, jBas-1
                         jSOj = jSO + jAOj
                         Do 420 iAOi = 0, iBas-1
                            iSOi = iSO + iAOi
                            mijkl = mijkl + 1
!
!---------------------------Contribution D(ij)*D(kl) to P(ijkl)
                            If (j1.eq.j2) Then
!------------------------------j3.eq.j4 also
                               Indi=Max(iSOi,jSOj)
                               Indj=iSOi+jSOj-Indi
                               Indk=Max(kSOk,lSOl)
                               Indl=kSOk+lSOl-Indk
                               iPntij=iPntSO(j1,j2,lOper,nbas)
                               iPntkl=iPntSO(j3,j4,lOper,nbas)
                               Indij=iPntij+(Indi-1)*Indi/2+Indj
                               Indkl=iPntkl+(Indk-1)*Indk/2+Indl
                               temp=DSO(Indij)*DSO(Indkl)*Coulfac
                            Else
                               temp = Zero
                            End If
!
!---------------------------Contribution -1/4*D(ik)*D(jl) to P(ijkl)
                            If (j1.eq.j3) Then
!------------------------------j2.eq.j4 also
                               Indi=Max(iSOi,kSOk)
                               Indk=iSOi+kSOk-Indi
                               Indj=Max(jSOj,lSOl)
                               Indl=jSOj+lSOl-Indj
                               iPntik=iPntSO(j1,j3,lOper,nbas)
                               iPntjl=iPntSO(j2,j4,lOper,nbas)
                               Indik=iPntik+(Indi-1)*Indi/2+Indk
                               Indjl=iPntjl+(Indj-1)*Indj/2+Indl
                               temp=temp-t14*(
     &                              DSO(Indik)*DSO(Indjl)
     &                             +DSSO(Indik)*DSSO(Indjl)
     &                                       )
                            End If
!
!---------------------------Contribution -1/4*D(il)*D(jk) to P(ijkl)
                            If (j1.eq.j4) Then
!------------------------------j2.eq.j3 also
                               Indi=Max(iSOi,lSOl)
                               Indl=iSOi+lSOl-Indi
                               Indj=Max(jSOj,kSOk)
                               Indk=jSOj+kSOk-Indj
                               iPntil=iPntSO(j1,j4,lOper,nbas)
                               iPntjk=iPntSO(j2,j3,lOper,nbas)
                               Indil=iPntil+(Indi-1)*Indi/2+Indl
                               Indjk=iPntjk+(Indj-1)*Indj/2+Indk
                               temp=temp-t14*(
     &                              DSO(Indil)*DSO(Indjk)
     &                             +DSSO(Indil)*DSSO(Indjk)
     &                                       )
                            End If
!
                            PMax=Max(PMax,Abs(Temp))
                            PSO(mijkl,MemSO2) =  temp
!
 420                     Continue
 320                  Continue
 220               Continue
 120            Continue
!
 410            Continue
 310         Continue
 210      Continue
 110   Continue
!
!
 400           Continue
 300        Continue
 200     Continue
 100  Continue
      If (nPSO.ne.MemSO2) Then
         Call WarningMessage(2,' PGet2: nPSO.ne.MemSO2')
         Write (6,*) nPSO, MemSO2
         Call Abend()
      End If
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' In PGet2:PSO ',' ',PSO,nijkl,nPSO)
#endif
      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_logical(Shijij)
      End SubRoutine PGet2
