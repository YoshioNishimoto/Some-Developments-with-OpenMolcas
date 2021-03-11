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
      subroutine two2mean34b(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
      implicit real*8 (a-h,o-z)
#include "para.fh"
      logical sameorb
      dimension
     *carteSO(norbsum,ncontmf,norbsum,ncontmf),
     *carteOO(norbsum,ncontmf,norbsum,ncontmf),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) then
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo
      coeff=0.5D0*coeff
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *carteSO(icartleft,jrun,icartright,irun)
      enddo
      enddo
      enddo
      enddo
      else
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo
      coeff=0.5D0*coeff
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *(carteSO(icartleft,jrun,icartright,irun)+
     *2d0*carteOO(icartleft,jrun,icartright,irun))
      enddo
      enddo
      enddo
      enddo
      endif
      return
      end
