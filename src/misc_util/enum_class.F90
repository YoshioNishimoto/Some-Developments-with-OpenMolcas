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
! Copyright (C) 2021, Oskar Weser                                      *
!***********************************************************************

module enum_class_mod
    implicit none
    private
    public :: EnumBase_t

    !> @brief An enum class in Fortran.
    !>
    !> @details
    !> It is typesafe and scope bound.
    !> (No comparison between integers and enum values for example.)
    !>
    !>   type, extends(EnumBase_t) :: Color_t
    !>   end type
    !>
    !>   type :: PossibleColors_t
    !>       type(Color_t) :: Red = Color_t(1), Blue = Color_t(2)
    !>   end type
    !>
    !>   write(stdout, *) possible_colors%Red == possible_colors%Blue

    type, abstract :: EnumBase_t
        integer :: val
    contains
        private
        procedure :: eq_EnumBase_t
        procedure :: neq_EnumBase_t
        generic, public :: operator(==) => eq_EnumBase_t
        generic, public :: operator(/=) => neq_EnumBase_t
    end type

contains

    logical elemental function eq_EnumBase_t(this, other)
        class(EnumBase_t), intent(in) :: this, other
        if (.not. SAME_TYPE_AS(this, other)) error stop 'Can only compare objects of same type'
        eq_EnumBase_t = this%val == other%val
    end function

    logical elemental function neq_EnumBase_t(this, other)
        class(EnumBase_t), intent(in) :: this, other
        if (.not. SAME_TYPE_AS(this, other)) error stop 'Can only compare objects of same type'
        neq_EnumBase_t = this%val /= other%val
    end function
end module
