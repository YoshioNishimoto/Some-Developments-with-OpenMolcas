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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************

      Subroutine Start_deepl(nPoints,nInter,x_,dy_,y_)
      use globvar
#include "stdalloc.fh"
!
        Integer nInter,nPoints
        Real*8 x_(nInter,nPoints),dy_(nInter,nPoints),y_(nPoints)
!
        Call mma_Allocate(x,nInter,nPoints,Label="x")
        Call mma_Allocate(dy,nInter*nPoints,Label="dy")
        Call mma_Allocate(y,nPoints,Label="y")

!
        return
      end
