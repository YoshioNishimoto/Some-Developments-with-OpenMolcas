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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************
*  Name_to_lm
*
*> @brief
*>   Get \f$ l \f$ and \f$ m \f$ numbers from a basis function name.
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Given a basis function name as a character string, extract the \f$ l
*> \f$ and \f$ m \f$ quantum numbers.
*>
*> Spherical harmonics functions are expected with the format `nnLmms`,
*> where `nn` is the shell number, `L` is a letter denoting angular
*> momentum, `mm` is the absolute value of the \f$ m \f$ number, and `s`
*> is the sign of \f$ m \f$.
*>
*> Cartesian functions are expected with the format `Lxxyyzz`, where `L`
*> is a letter denoting angular momentum, and `xx`, `yy`, `zz` are the
*> powers of \f$ x \f$, \f$ y \f$ and \f$ z \f$ (\f$ m_x,m_y,m_z \f$).
*> in this case, the numbers returned are \f$ -l \f$ and
*> \f$ ((m_y+m_z)^2+m_z-m_y)/2-m_x \f$.
*>
*> @param[in]  BName Basis function name
*> @param[out] l     \f$ l \f$ number (\f$ -l \f$ for Cartesians)
*> @param[out] m     \f$ m \f$ number (see details for Cartesians)
************************************************************************
      Subroutine Name_to_lm(BName,l,m)
      Implicit None
      Character(Len=*), Intent(In) :: BName
      Integer, Intent(Out) :: l, m
      Character :: Letter
      Integer :: i, lx, ly, lz
#include "angtp.fh"
*
      Letter = BName(3:3)
      Call LoCase(Letter)
      l = 0
      m = 0
      If (Letter .eq. 's') Then
*     Default is s
        Return
      Else if (Letter .eq. 'p') Then
*       p usually appear as px, py, pz, except when they are contaminants
        l = 1
        If (BName(4:4) .ne. '0') Then
          Letter = BName(4:4)
          Call LoCase(Letter)
          If (Letter .eq. 'x') Then
            m = 1
          Else If (Letter .eq. 'y') Then
            m = -1
          Else If (Letter .eq. 'z') Then
            m = 0
          End If
          Return
        End If
      End If
*     Parse the label for other cases
      l = -1
*     Find if there is an angular label
      Do i=Sum(LBound(AngTp)),Sum(UBound(AngTp))
        If (Letter .eq. AngTp(i)) Then
          l = i
          Exit
        End If
      End Do
      If (l .ge. 0) Then
*       If a label is found it is a spherical shell, just read m
        Read(BName(4:5),*) m
        If (BName(6:6) .eq. '-') m = -m
      Else
*       If no label, this is a Cartesian shell, return -l and some convention for m.
*       We use m=T(ly+lz)-(lx+ly), where T(n) is the nth triangular number: n*(n+1)/2).
*       From here, ly+lz can be recovered as the (int) triangular root of m+l: (sqrt(8*(m+l)+1)-1)/2,
*       lz is m+l-T(ly+lz) and lx is l-(ly+lz). This has the property that all
*       possible combinations of lx,ly,lz are encoded in -l plus a number from -l to l*(l+1)/2,
*       in descending order with priority lx>ly>lz: (3,0,0), (2,1,0), (2,0,1), (1,2,0),
*       (1,1,1), (1,0,2), (0,3,0), (0,2,1), (0,1,2), (0,0,3)
        Read(BName(2:3),*) lx
        Read(BName(4:5),*) ly
        Read(BName(6:7),*) lz
        l = -lx-ly-lz
        m = (ly+lz)*(ly+lz+1)/2-(lx+ly)
      End If
      Return
*
      End Subroutine Name_to_lm
