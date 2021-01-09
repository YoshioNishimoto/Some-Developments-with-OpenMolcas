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
* Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
*               1990, IBM                                              *
************************************************************************
      SubRoutine OneEl(Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,
     &                 nOrdOp,rNuc,rHrmt,iChO,
     &                 opmol,ipad,opnuc,iopadr,idirect,isyop,
     &                 PtChrg,nGrid,iAddPot)
      use Basis_Info, only: nBas
      use PrpPnt
      use Temporary_Parameters, only: PrPrt, Short, IfAllOrb
      use Sizes_of_Seward, only: S
      use Real_Info, only: Thrs
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
      External Kernel, KrnlMm
#include "stdalloc.fh"
#include "real.fh"
      Real*8, Dimension(:), Allocatable :: Out, Nuc, TMat, Temp, El,
     &                                     Array
      Character*16, Dimension(:), Allocatable :: plabs
      Character Label*8, LBL*4
      Character L_Temp*8

C-tbp-npy-start
      character filename*16
C-tbp-npy-end

      Real*8 CCoor(3,nComp), rNuc(nComp), PtChrg(nGrid)
      Integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7)
      dimension opmol(*),opnuc(*),iopadr(ncomp,*)
      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) ' In OneEl: Label', Label
      Write (6,*) ' In OneEl: nComp'
      Write (6,'(1X,8I5)') nComp
      Write (6,*) ' In OneEl: lOper'
      Write (6,'(1X,8I5)') lOper
      Write (6,*) ' In OneEl: n2Tri'
      Do iComp = 1, nComp
         ip(iComp) = n2Tri(lOper(iComp))
      End Do
      Write (6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
      Call RecPrt(' CCoor',' ',CCoor,3,nComp)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Compute the number of blocks from each component of the operator
*     and the irreps it will span.
*
      nIC = 0
      llOper = 0
      Do iComp = 1, nComp
         llOper = iOr(llOper,lOper(iComp))
         Do iIrrep = 0, nIrrep-1
            If (iAnd(lOper(iComp),iTwoj(iIrrep)).ne.0) nIC = nIC + 1
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*) ' nIC =',nIC
#endif
*
      If (nIC.eq.0) Return
*
      Call SOS(iStabO,nStabO,llOper)
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate memory for symmetry adapted one electron integrals.
*     Will just store the unique elements, i.e. low triangular blocks
*     and lower triangular elements in the diagonal blocks.
*
      Call ICopy(nComp,[-1],0,ip,1)
      LenTot=0
      Do iComp = 1, nComp
         LenInt=n2Tri(lOper(iComp))
         LenTot=LenTot+LenInt+4
      End Do
      Call mma_allocate(Array,LenTot,label='Array')
      ip(1)=1
      call dcopy_(LenTot,[Zero],0,Array(ip(1)),1)
      iadr=ip(1)
      do iComp = 1, nComp
         LenInt=n2Tri(lOper(iComp))
         ip(icomp)=iadr
         iadr=iadr+LenInt+4
*        Copy center of operator to work area.
         call dcopy_(3,Ccoor(1,iComp),1,Array(ip(iComp)+LenInt),1)
*        Copy nuclear contribution to work area.
         Array(ip(iComp)+LenInt+3) = rNuc(iComp)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute all SO integrals for all components of the operator.
*
      Call OneEl_Internal
     &           (Kernel,KrnlMm,Label,ip,lOper,nComp,CCoor,
     &            nOrdOp,rHrmt,iChO,
     &            opmol,opnuc,ipad,iopadr,idirect,isyop,
     &            iStabO,nStabO,nIC,
     &            PtChrg,nGrid,iAddPot,Array,LenTot)
*                                                                      *
************************************************************************
*                                                                      *
*                    P O S T P R O C E S S I N G                       *
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Call PrMtrx(Label,lOper,nComp,ip,Array)
#endif

*                                                                      *
************************************************************************
*                                                                      *
*     Make a square sum on all the integrals for verification
      Call VrfMtrx(Label,lOper,nComp,ip,Array)
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute properties or write integrals to disc.
*
      Do iComp = 1, nComp
         iSmLbl = lOper(iComp)
         If (Prprt) Then
*                                                                      *
************************************************************************
*                                                                      *
*---------- Compute properties directly from integrals
*
*---------- Allocate some memory
*
            If (iComp.eq.1) Then
               If (short) Then
                  mDim = 1
               Else
                  mDim = S%nDim
               End If
               call mma_allocate(Out,mDim*nComp,label='Out')
               ipOut=1
               call dcopy_(mDim*nComp,[Zero],0,Out,1)
               call mma_allocate(Nuc,nComp,label='Nuc')
               ipNuc=1
               call dcopy_(nComp,[Zero],0,Nuc,1)
            End If
            nInt=n2Tri(iSmLbl)
            If (nInt.ne.0)
     &      Call CmpInt(Array(ip(iComp)),nInt,nBas,nIrrep,iSmLbl)
            Nuc(ipNuc+(iComp-1)) = Array(ip(iComp)+nInt+3)
            If (nInt.ne.0)
     &      Call XProp(Short,ifallorb,
     &                 nIrrep,nBas,nVec,Vec,nOcc,Occ,
     &                 Thrs,nDen,Array(ip(iComp)),
     &                 Out(ipOut+(iComp-1)*mDim))
*
            If (Label(1:3).eq.'PAM') Then
c               Open(unit=28,file='R_vect',access='append')
               EndFile(28)
               If(Short) Then
                  write(28,'(a8,2x,f20.14)')
     &                 Label, Out(ipOut+(iComp-1)*mDim)
               Else
                  Sum = 0.00d0
                  Do ii= 1 , nOcc
                    Sum = Sum + Out(ipOut+ii-1+(iComp-1)*mDim)
                  End Do
                  write(28,'(a8,2x,f20.14)') Label,-Sum
               End If
c               Close(28)
             End If
*
*---------- Once all components have been computed print them.
*
            If (iComp.eq.nComp) Then
               LBL=Label(1:4)
               Call UpCase(LBL)
               lpole = 0
               If (LBL.eq.'MLTP') Then
                   Read(Label,'(5X,I3)') lpole
               Else If (LBL.eq.'PAM ') Then
                   Read(Label,'(5X,I1)') lpole
               Else If (LBL.eq.'L_MP') Then
                   Read(Label,'(5X,I1)') lpole
               Else If (LBL(1:2).eq.'EF') Then
                   Read(Label,'(2X,I1)') lpole
               Else If (LBL.eq.'DMS ') Then
                   lpole = 3
               Else If (LBL.eq.'VELO') Then
                   lpole = 1
               End If
               Call mma_allocate(plabs,nComp,label='plabs')
               Call mma_allocate(TMat,nComp**2,label='TMat')
               Call mma_allocate(Temp,nComp,label='Temp')
               If (nComp.eq.1) Then
                  ipC2 = 1 ! dummy
               Else
                  ipC2 = 2 ! Used only for diamagnetic shielding.
               End If
               Call Prop(Short,Label,Ccoor(1,1),Ccoor(1,ipC2),
     &                   nIrrep,nBas,mDim,Occ,Thrs,
     &                   Out,Nuc,lpole,plabs,TMat,Temp,ifallorb)
*
* For a properties calculation, save the values of EF or CNT operators,
* they will be used to write the sum through Add_Info in Drv1El
*
               If (PrPrt.and.
     &             (LBL(1:2).eq.'EF'.or.LBL(1:3).eq.'CNT')) Then
                 Call mma_allocate(El,nComp,label='El')
                 ipEl=1
                 Call FZero(El,nComp)
*                Compute the sum of all orbital components
                 Do jComp=0,nComp-1
                   iInd1=ipEl+jComp
                   iInd2=ipOut+jComp*mDim
                   Do iOcc=0,mDim-1
                     El(iInd1)=El(iInd1)+Out(iInd2+iOcc)
                   End Do
                 End Do
*                Write electronic and nuclear components in temp file
                 LuTmp=10
                 Call DaName(LuTmp,'TMPPRP')
                 Read(Label(4:8),*) iEF
                 iDisk=(iEF-1)*2
                 Call dDaFile(LuTmp,1,El,nComp,iDisk)
                 Call dDaFile(LuTmp,1,Nuc,nComp,iDisk)
                 Call DaClos(LuTmp)
                 Call mma_deallocate(El)
               End If
*
               Call mma_deallocate(Temp)
               Call mma_deallocate(TMat)
               Call mma_deallocate(plabs)
               Call mma_deallocate(Nuc)
               Call mma_deallocate(Out)
            End If
         Else
*                                                                      *
************************************************************************
*                                                                      *
*---------- Write integrals to disc
*
            iOpt = 0
            iRC = -1
            If (Label(1:3).eq.'PAM') Then
               Write(L_Temp,'(A5,I3.3)') 'PAM  ',iPAMcount
               iPAMcount=iPAMcount+1
               iComp_=1
            Else If (Label(1:5).eq.'EMFR0') Then
               iComp_=1
               If (iComp.eq.1) Then
                  L_Temp='EMFR0  R'
               Else
                  L_Temp='EMFR0  I'
               End If
            Else If (Label(1:5).eq.'EMFR ') Then
               iComp_=MOD(iComp+2,3)+1
               If ((iComp+2)/3.eq.1) Then
                  L_Temp='EMFR  RS'
               Else If ((iComp+2)/3.eq.2) Then
                  L_Temp='EMFR  RA'
               Else If ((iComp+2)/3.eq.3) Then
                  L_Temp='EMFR  IS'
               Else If ((iComp+2)/3.eq.4) Then
                  L_Temp='EMFR  IA'
               End If
             Else If (Label(1:5).eq.'TMOM0') Then
               If (iComp.eq.1) Then
                  L_Temp='TMOM0  R'
                  iComp_=1
               Else
                  L_Temp='TMOM0  I'
                  iComp_=1
               End If
             Else If (Label(1:5).eq.'TMOM2') Then
               If (iComp.eq.1) Then
                  L_Temp='TMOM2  R'
                  iComp_=1
               Else
                  L_Temp='TMOM2  I'
                  iComp_=1
               End If
            Else If (Label(1:5).eq.'TMOM ') Then
               iComp_=MOD(iComp+2,3)+1
               If ((iComp+2)/3.eq.1) Then
                  L_Temp='TMOM  RS'
               Else If ((iComp+2)/3.eq.2) Then
                  L_Temp='TMOM  RA'
               Else If ((iComp+2)/3.eq.3) Then
                  L_Temp='TMOM  IS'
               Else If ((iComp+2)/3.eq.4) Then
                  L_Temp='TMOM  IA'
               End If
            Else
               L_Temp=Label
               iComp_=iComp
            End If
C-tbp-npy-start
            write(filename,'(A8,A1,I3,A4)') L_Temp,'_',iComp_,'.npy'
            do i = 1,12
               if (filename(i:i) .eq. ' ') then
                  filename(i:i) = '_'
               end if
            end do
            call npy_write_dbl_vec(filename, n2Tri(lOper(iComp)),
     *                             Array(ip(iComp)))
C-tbp-npy-end

            Call WrOne(iRC,iOpt,L_Temp,iComp_,Array(ip(iComp)),iSmLbl)

            If (iRC.ne.0) then
               Call WarningMessage(2,
     &               ' *** Error in subroutine ONEEL ***,'//
     &               '     Abend in subroutine WrOne')
               Call Abend()
            End If
         End If
      End Do  ! iComp
*
      if (Label.eq.'Attract ')
     &   Call Add_info('SEWARD_ATTRACT',Array(ip(1)),1,5)
      if (Label.eq.'Kinetic ')
     &   Call Add_info('SEWARD_KINETIC',Array(ip(1)),1,5)
      if (Label.eq.'Mltpl  1')
     &   Call Add_info('SEWARD_MLTPL1X',Array(ip(1)),1,5)
      If (Label(1:4).eq.'EMFR') Then
         CheckSum=Zero
         Do iComp = 1, nComp
            CheckSum=CheckSum
     &              +DDot_(n2Tri(lOper(iComp)),Array(ip(iComp)),1,One,0)
         End Do
         Call Add_info('SEWARD_EMFR',CheckSum,1,6)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Deallocate memory for integral
*
      Call mma_deallocate(Array)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
C-tbp-npy-start
      subroutine npy_write_dbl_vec(filename, n, vec)
C
C     Thomas Bondo Pedersen, Dec. 2020
C
C     Write double-precision vector to npy file
C
C     Note:
C     - uses stream io (f2003 standard)
C     - uses 'i0' format specifier (f95 standard)
C     - avoids gpopen/gpclose but uses the same
C       table to find a free unit number
C
C
      implicit none
      character*(*) filename
      integer n
      real*8 vec(n)

      integer*4 header_len
      integer lunit, i, l
      logical done

      integer maxlen_dict_str
      parameter (maxlen_dict_str = 85)
      character*(maxlen_dict_str) dict_str
      character*5 magic_str
      character magic_num, major, minor

      integer isfreeunit
      external isfreeunit

      if (n .le. 0) return ! nothing to write

      magic_str = 'NUMPY'
      magic_num = achar(147) ! x93
      major = achar(2) ! npy file format, major version
      minor = achar(0) ! npy file format, minor version


      ! note that the format specifier "I0" will only work properly
      ! with a fortran95-compliant compiler...
      write(dict_str,'(A,I0,A)')
     * "{'descr': '<f8', 'fortran_order': True, 'shape': (",
     * n,",), }"

      l = maxlen_dict_str
      done = .false.
      do while (.not. done .and. l .gt. 54)
         if (dict_str(l:l) .eq. '}') then
            done = .true.
         else
            l = l - 1
         end if
      end do
      if (.not. done) then
         write(6,*) 'npy_write_dbl_vec: unable to construct header'
         call abend()
      end if
      do while (l.lt.maxlen_dict_str .and. mod(l + 11, 16) .ne. 0)
          l = l + 1
          dict_str(l:l) = " "
      end do
      dict_str(l:l) = achar(10)
      header_len = INT(l,4)

      l = 7
      lunit = isfreeunit(l)
      open(unit=lunit, file=filename, form='unformatted',
     *     access='stream')
      write(lunit) magic_num, magic_str, major, minor
      write(lunit) header_len
      write(lunit) dict_str(1:header_len)
      write(lunit) (vec(i), i=1,n)
      close(unit=lunit)

      return
      end
C-tbp-npy-end
