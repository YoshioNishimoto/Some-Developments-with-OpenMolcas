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
* Copyright (C) Yannick Carissan                                       *
************************************************************************
      subroutine RotateGeoms(Q)
      implicit none
#include "geoms.fh"
#include "debug.fh"
#include "WrkSpc.fh"
#include "options.fh"
      Real*8 Q(0:3)
      integer igeom

      Do igeom=3,ngeoms+2
        if (debug) then
          Write(6,*) "Before rotation"
          call PrintGeom(6,nat(igeom),title(igeom),
     &           Work(ipgeo(igeom)),geolbl(1,igeom))
        End If
        Call RotateGeom(Q,nat(igeom),
     &           Work(ipgeo(igeom)),Work(ipgeo(igeom)))
        if (debug) then
          Write(6,*) "After rotation"
          call PrintGeom(6,nat(igeom),title(igeom),
     &           Work(ipgeo(igeom)),geolbl(1,igeom))
        End If
      End Do
      Return
      End
