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
* Copyright (C) 2020, Roland Lindh                                     *
************************************************************************
      Subroutine Dispersion_Kriging_Layer(qInt,E_Disp,nInter)
      Use Limbo
      Implicit None
#include "stdalloc.fh"
      Integer nInter
      Real*8 qInt(nInter), E_Disp
      Real*8, Allocatable:: qInt_s(:)
*
      Call mma_allocate(qInt_s,nInter,Label='qInt_s')
*
      Call Trans_K(U,qInt,qInt_s,nInter,1)
*#define _DEBUG_
#ifdef _DEBUG_
      Call RecPrt('Dispersion_Kriging_Layer: qInt',' ',qInt,nInter,1)
      Call RecPrt('Dispersion_Kriging_Layer: qInt_s',' ',qInt_s,
     &                                                   nInter,1)
#endif
      Call Dispersion_Kriging(qInt_s,E_Disp,nInter)
*
      Call mma_deallocate(qInt_s)
*
      End Subroutine Dispersion_Kriging_Layer
