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
       Subroutine Get_Temp(Label,Temp,nTemp)
       Implicit Real*8 (a-h,o-z)

       Real*8        Temp(nTemp)
       Character*(*) Label

       Call NameRun('.RUNTMP')
       Call Get_dArray(Label,Temp,nTemp)
       Call NameRun('#Pop')

       Return
       End
