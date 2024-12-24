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
        subroutine makeT2pHlp1 (T2p,Tau,aGrp,bGrp,aSGrp,bSGrp,key,
     c                          dimi,dimij,dimapp,dimabpp,dimap,dimabp)
c
c       this routine do:
c       define T2(+-)(ij,(ab)") = Tau((ab)',i,j) +- Tau((ab)',j,i)
c       for the case: aGrp=bGrp and aSGrp=bSGrp
c
c       parameter description:
c       T2p     - array for T2+ (O)
c       Tau     - array for Tau (I)
c       xGrp    - Groups of a',b' (I)
c       xSGrp   - SubGroups of a",b" (I)
c       key     - 0 - T2+; 1 = T2- will be produced
c       dimx    - Dimension of i,(i>=j),a",(a>b)",a',(a>=b)' (I)
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aGrp,bGrp,aSGrp,bSGrp,key
        integer dimi,dimij,dimapp,dimabpp,dimap,dimabp
        real*8 T2p(1:dimij,1:dimabpp)
        real*8 Tau(1:dimabp,1:dimi,1:dimi)
c
c
c       help variables
        integer i,j,ij,app,bpp,abpp,ap,abp,appAdd,bppAdd
c
c
c1      def appAdd,bppAdd
c
        appAdd=0
        if (aSGrp.ne.Grpalow(aGrp)) then
          do i=Grpalow(aGrp),aSGrp-1
          appAdd=appAdd+DimSGrpa(i)
          end do
        end if
c
        bppAdd=0
        if (bSGrp.ne.Grpalow(bGrp)) then
          do i=Grpalow(bGrp),bSGrp-1
          bppAdd=bppAdd+DimSGrpa(i)
          end do
        end if
c
c
c2      define T2(+,-)
c
        if (key.eq.0) then
c
c2.1    define T2+
        abpp=0
        do app=2,dimapp
        ap=appAdd+app
        abp=ap*(ap-1)/2+bppAdd
        do bpp=1,app-1
        abp=abp+1
        abpp=abpp+1
          ij=0
          do i=1,dimi
          do j=1,i
            ij=ij+1
            T2p(ij,abpp)=Tau(abp,i,j)+Tau(abp,j,i)
          end do
          end do
        end do
        end do
c
        else
c
c2.2    define T2-
        abpp=0
        do app=2,dimapp
        ap=appAdd+app
        abp=ap*(ap-1)/2+bppAdd
        do bpp=1,app-1
        abp=abp+1
        abpp=abpp+1
          ij=0
          do i=2,dimi
          do j=1,i-1
            ij=ij+1
            T2p(ij,abpp)=Tau(abp,i,j)-Tau(abp,j,i)
c           T2p(ij,abpp)=T2c(ap,bpp,i,j)-T2c(ap,bpp,j,i)
          end do
          end do
        end do
        end do
c
        end if
c
        call mv0sv (dimij*dimabpp,dimij*dimabpp,
     c              T2p(1,1),0.5d0)
c
        return
c Avoid unused argument warnings
        if (.false.) call Unused_integer(dimap)
        end
