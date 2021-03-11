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
* Copyright (C) Roland Lindh                                           *
************************************************************************
*  Get_AnalHess
*
*> @brief
*>   Read the the symmetry blocked nuclear Hessian from the runfile and return a
*>   pointer to the array's location in \c Work and the length of the array
*> @author R. Lindh
*>
*> @details
*> The utility will read the symmetry blocked nuclear Hessian from the runfile and
*> return a pointer and the length of the array.
*>
*> @param[out] ipAnalHess pointer to array with the symmetry blocked nuclear Hessian
*>                        in Cartesian coordinates
*> @param[out] nAnalHess  size of the array of the symmetry blocked nuclear Hessian
************************************************************************
      Subroutine Get_AnalHess(ipAnalHess,nAnalHess)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character*24 Label
      Logical      Found

      Label='Analytic Hessian'
      Call qpg_dArray(Label,Found,nAnalHess)
      If(.not.Found .or. nAnalHess.eq.0) Return
      Call GetMem('AnalHess','Allo','Real',ipAnalHess,nAnalHess)
      Call Get_dArray(Label,Work(ipAnalHess),nAnalHess)

      Return
      End
