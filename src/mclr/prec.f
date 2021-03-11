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
* Copyright (C) 1996,1997, Anders Bernhardsson                         *
************************************************************************
      SubRoutine Prec(rpre,idsym)
************************************************************************
*
*  idsym, symmetry of orbital hessian of interest
*  CMtx preconditioner
*
*
* The orbital hessian is dominated of elements that couples
*
* kappa  -> kappa      where i is occupied and p,q is general.
*      ip        iq
*
* we therefore approximate the hessian with thoose diagonal
* terms in the preconditioner
*
*  Anders Bernhardsson 96
*
*     active; active,general is needed for rasscf calculation
*     and is not coded yet (ugly bastard) (970109, AB )
************************************************************************
      Implicit Real*8(a-h,o-z)
#include "Pointers.fh"
#include "Input.fh"
#include "WrkSpc.fh"
#include "machine.fh"
      Real*8 rpre(*)
*
      Call Prec_internal(rpre)
*
*     This is to allow type punning without an explicit interface
      Contains
      Subroutine Prec_internal(rpre)
      Use Iso_C_Binding
      Real*8, Target :: rpre(*)
      Integer, Pointer :: ipre(:)
      nmm=0
      nmmm=0
      Do iS=1,nSym
         nMM=Max(nMM,nAsh(is)+nIsh(iS))
         nMMM=Max(nmmM,nBas(is))
      End Do
      n2=nMMM**2
      nmmm=((nmmm-1)/nRec+1)*nRec
      nmm=nmm*nMMM
      nmm=nmm**2

*
      Call GetMem('JInt','Allo','Real',ipJ,n2)
      Call GetMem('KInt','Allo','Real',ipK,n2)
      Call GetMem('Scr ','Allo','Real',ipS,n2)
      If (ActRot.or.nRs1(iDSym).ne.0.or.nRs3(iDSym).ne.0) Then
C       write (*,*) "ntash = ", ntash
        Call GetMem('AINT','ALLO','REAL',ipActInt,ntAsh**4)
        Call Precaaa_Pre(Work(ipActInt),Work(ipJ),Work(ipS))
      End If
*
      ip=1
      iAdr=0
      iAdr2=0
      Do iS=1,nSym
         jS=iEOr(is-1,iDSym-1)+1
         nD=nOrb(js)-nIsh(jS)
         ni=nBas(js)**2
         sign=1.0d0
         Call GetMem('2TEMP2','ALLO','REAL',ipTemp2,ni)
         Call GetMem('3TEMP3','ALLO','REAL',ipTemp3,ni)
         Call GetMem('1Temp1','MAX','Real',ipT1,nTemp)
         nTemp=Min(nmm,nTemp/2)
         Call GetMem('1Temp1','ALLO','Real',ipTemp1,2*nTemp)
         ipScr=ipTemp1+nTemp
         If (nd.eq.0) Goto 100
         Do iB=1,nIsh(iS)
            call dcopy_(nD**2,[0.0d0],0,Work(ipTemp3),1)
            ibb=nOrb(is)*(ib-1)+ib-2
*
**          Cholesky code
*
            If (newCho) Then
              Call preci_cho(ib,is,jS,nD,Work(ipTemp3),
     &                       nOrb(is),nOrb(js),
     &                       Work(ipFIMO+ipCM(is)+ibb),
     &                       Work(ipFAMO+ipCM(is)+ibb),
     &                       Work(ipF0sqMO+ipCM(is)+ibb),
     &                       Work(ipFIMO+ipCM(js)-1),
     &                       Work(ipFAMO+ipCM(js)-1),
     &                       Work(ipF0sqMO+ipCM(js)-1),sign,
     &                       Work(ipJ),Work(ipK),Work(ipS),n2,
     &                       iAdr) ! OK

            Else
            If (iMethod.eq.2) Then
*                                                                      *
************************************************************************
*                                                                      *
*              G
*               iaib
*
               If (nash(js).gt.0)
     &            Call Preciaa(ib,is,js,nd,Work(ipTemp3),
     &                         nOrb(is),nOrb(js),
     &                         Work(ipFIMO+ipCM(is)+ibb),
     &                         Work(ipFAMO+ipCM(is)+ibb),
     &                         Work(ipF0sqMO+ipCM(is)+ibb),
     &                         Work(ipFIMO+ipCM(js)-1),
     &                         Work(ipFAMO+ipCM(js)-1),
     &                         Work(ipF0sqMO+ipCM(js)-1),sign,
     &                         Work(ipJ),Work(ipK),Work(ipS),n2) ! OK
*                                                                      *
************************************************************************
*                                                                      *
*              G
*               ipia
*
               If ((nOrb(js)-nish(js)-nash(js))*nash(js).gt.0)
     &            Call Preciba(ib,is,js,nd,Work(ipTemp3),nOrb(js),
     &                         Work(ipFIMO+ipCM(js)-1),
     &                         Work(ipFAMO+ipCM(js)-1),
     &                         Work(ipF0sqMO+ipCM(js)-1),sign,
     &                         Work(ipJ),Work(ipK),Work(ipS),n2) ! OK
*
            End If
*                                                                      *
************************************************************************
*                                                                      *
*           G
*            ipiq
*
            If ((nOrb(js)-nish(js)-nash(js)) .gt.0)
     &         Call Precibb(ib,is,js,nd,Work(ipTemp3),
     &                      nbas(js),norb(js),
     &                      Work(ipTemp1),Work(ipScr),Work(ipTemp2),
     &                      Work(ipFiMo+ipCM(is)+ibb),
     &                      Work(ipFAMO+ipcm(is)+ibb),  ! OK
     &                      Work(ipFiMo+ipCM(js)-1),
     &                      Work(ipFAMO+ipcm(js)-1),sign)  ! OK
            EndIf ! newCho
*                                                                      *
************************************************************************
*                                                                      *
*           Factorize G:
*
*               T
*           G=LL
*
*            write(6,*) 'Preconditioner i =',iB
*            Do i=1,min(nd,10)
*             write(6,'(10F12.8)') (Work(ipTemp3+(j-1)*(2*nd-j+2)/2+i-j),
*     &                             j=1,i)
*            End Do

            Call SQM(Work(ipTemp3),rpre(ip),nd)
C     write (*,*) "ib=",ib
C     call sqprt(rpre(ip),nd)

!            write(*,*)" ====== rpre ====== "
!            do i=1,nd*nd
!              write(*,*)i,"rpre",rpre(ip+i-1)
!            end do

            irc=0
            call c_f_pointer(c_loc(rpre(ip+nd**2)),ipre,[nd])
            call dgetrf_(nd,nd,rpre(ip),nd,ipre,irc)
            nullify(ipre)
            If (irc.ne.0) then
               Write(6,*) 'Error in DGETRF called from prec'
               Call Abend
            End If
            ip=ip+nD*(nd+1)
*
         End Do

 100     Continue
*                                                                      *
************************************************************************
*                                                                      *
         Do iB=1,nAsh(iS)
            ibb=nOrb(is)*(nish(is)+ib-1)+nish(is)+ib-2
            If (ib.le.nRs1(iS)+nRs2(is)+nRs3(is)) iR=3
            If (ib.le.nRs1(iS)+nRs2(is)) iR=2
            If (ib.le.nRs1(iS)) iR=1
            If (ActRot) Then
               nD=nOrb(js)
            Else
               If (ir.eq.1) nD=nOrb(js)-nRs1(js)
               If (ir.eq.2) nD=nOrb(js)-nRs2(js)
               If (ir.eq.3) nD=nOrb(js)-nRs3(js)
            End If
            If (nd.eq.0) Goto 110
            call dcopy_(nD**2,[0.0d0],0,Work(ipTemp3),1)
            ndtri=nd*(nd+1)/2
*
**  New Cholesky code
*
            If (newCho) Then
               Call Preca_cho(ib,is,js,nd,ir,Work(ipTemp3),
     &                        nOrb(is),nOrb(js),
     &                        Work(ipFIMO+ipCM(is)+ibb),
     &                        Work(ipFAMO+ipCM(is)+ibb),
     &                        Work(ipF0SqMO+ipCM(is)+ibb),
     &                        Work(ipFIMO+ipCM(js)-1),
     &                        Work(ipFAMO+ipCM(js)-1),
     &                        Work(ipF0SqMO+ipCM(js)-1),sign,
     &                        Work(ipJ),Work(ipK),Work(ipS),n2,
     &                        iAdr2)
            Else
            If (nish(js).gt.0)
     &         Call Precaii(ib,is,js,nd,ir,Work(ipTemp3),
     &                      nOrb(is),nOrb(js),
     &                      Work(ipFIMO+ipCM(is)+ibb),
     &                      Work(ipFAMO+ipCM(is)+ibb),
     &                      Work(ipF0SqMO+ipCM(is)+ibb),
     &                      Work(ipFIMO+ipCM(js)-1),
     &                      Work(ipFAMO+ipCM(js)-1),
     &                      Work(ipF0SqMO+ipCM(js)-1),sign,
     &                      Work(ipJ),Work(ipK),Work(ipS),n2) ! OK
*           Call Precaai(ib,nd,ir,rpre(ip))
*           Call Precaaa(ib,nd,ir,rpre(ip))
            If (nish(js)*nOrb(js).gt.0)
     &         Call Precabi(ib,is,js,ir,nd,Work(ipTemp3),nOrb(js),
     &                      Work(ipFIMO+ipCM(js)-1),
     &                      Work(ipFAMO+ipCM(js)-1),
     &                      Work(ipF0SQMO+ipCM(js)-1),sign,
     &                      Work(ipJ),Work(ipK),Work(ipScr),n2) !+/-?
*           Call Precaba(ib,nd,ir,rpre(ip))
            If (nOrb(js).gt.0)
     &            Call Precabb_2(ib,is,js,nd,nbas(js),nOrb(js),
     &                           Work(ipTemp3),
     &                           Work(ipTemp1),ntemp,Work(ipScr),
     &                           Work(ipTemp2),
     &                           Work(ipF0SQMO+ipCM(is)+ibb),
     &                           Work(ipFiMo+ipCM(js)-1),
     &                           Work(ipFAMO+ipcm(js)-1) ,
     &                           Work(ipF0SQMO+ipCM(js)-1),sign)
*
            !! symmetry not yet
            !! Eq. (C.12e)
            If (ActRot.or.nRs1(iS).ne.0.or.nRs3(iS).ne.0)
     &         Call Precaaa(ib,is,js,nd,ir,Work(ipTemp3),
     &                      nOrb(is),nOrb(js),
     &                      Work(ipFIMO+ipCM(js)-1),
     &                      Work(ipF0SqMO+ipCM(js)-1),sign,
     &                      Work(ipJ),Work(ipK),Work(ipS),n2,
     &                      Work(ipActInt)) ! OK
            EndIf ! newCho

            Call SQM(Work(ipTemp3),rpre(ip),nD)
C      write (*,*) "ib = ", ib,nd
C      call sqprt(rpre(ip),nd)
            irc=0
            call c_f_pointer(c_loc(rpre(ip+nd**2)),ipre,[nd])
            call dgetrf_(nd,nd,rpre(ip),nd,ipre,irc)
C      write (*,*) "after LU"
C      call sqprt(rpre(ip),nd)
            nullify(ipre)
            If (irc.ne.0) then
               Write(6,*) 'Error in DGETRF called from prec'
               Call Abend
            End If
            ip=ip+nD*(nd+1)
         End Do
110      Continue
         Call GetMem('1TEMP1','FREE','REAL',ipTemp1,nOrb(js)**2)
         Call GetMem('2TEMP2','FREE','REAL',ipTemp2,nOrb(js)**2)
         Call GetMem('3TEMP3','FREE','REAL',ipTemp3,nOrb(js)**2)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (ActRot.or.nRs1(iDSym).ne.0.or.nRs3(iDSym).ne.0) Then
        Call GetMem('AINT','FREE','REAL',ipActInt,ntAsh**4)
      End If
      Call GetMem('Scr ','Free','Real',ipS,n2)
      Call GetMem('KInt','Free','Real',ipK,n2)
      Call GetMem('JInt','Free','Real',ipJ,n2)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End Subroutine Prec_internal
*
      End
