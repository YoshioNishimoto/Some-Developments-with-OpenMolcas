#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2020, Ignacio Fdez. Galván                             *
#***********************************************************************

# Template for building an OpenMolcas utility library,
# that will be part of libmolcas
# * Creates:
#   - An object library ${prog} with all sources
#     (except those that define Fortran modules)
# * Supports:
#   - ${util}_defs with utility-specific compile definitions
#   - ${util}_incs with utility-specific include directories
#   - ${util}_deps with utility-specific dependencies
# * Defines:
#   - ${util}_src with the path for the source directory
#   - ${util}_sources with all the non-module source files
#   - ${util}_mods with all the module source files

get_filename_component (util ${CMAKE_CURRENT_SOURCE_DIR} NAME)

if (DEFINED sources)
  set_absolute_paths (sources ${CMAKE_CURRENT_SOURCE_DIR} ${sources})
else ()
  file (GLOB sources *.f *.f90 *.F *.F90 *.c)
endif ()

# Ignore an empty directory
if (NOT sources)
  return ()
endif ()

# Find source files that define Fortran modules.
# These must be moved to another target to ensure
# that they are compiled first.
# If modfile_list is predefined, it may contain
# only the module sources that need to be available
# to other directories.
#--------------------------------------------------
if (DEFINED modfile_list)
  set_absolute_paths (modfile_list ${CMAKE_CURRENT_SOURCE_DIR} ${modfile_list})
else ()
  foreach (fname ${sources})
    file (READ ${fname} contents)
    if ("${contents}" MATCHES "\n *[mM][oO][dD][uU][lL][eE] ")
      list (APPEND modfile_list ${fname})
    endif ()
  endforeach ()
endif ()
foreach (fname ${modfile_list})
  list (REMOVE_ITEM sources ${fname})
endforeach ()

# Create an object library
#-------------------------
add_Molcas_library (${util}_obj OBJECT ${sources})

# dependencies
if (DEFINED ${util}_deps)
  add_dependencies(${util}_obj ${${util}_deps})
endif()
# utility-specific compile definitions
if (DEFINED ${util}_defs)
  target_compile_definitions (${util}_obj PRIVATE "${${util}_defs}")
endif ()
# utility-specific include directories
list (APPEND ${util}_incs ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (${util}_obj PRIVATE "${${util}_incs}")
# add a dummy source if there are none left,
# (CMake does not like empty targets)
if (NOT sources)
  set (dummy_file ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${util}/dummy.f90)
  if (NOT EXISTS ${dummy_file})
    file (WRITE ${dummy_file} "subroutine empty_dummy() ; end subroutine")
  endif ()
  target_sources (${util}_obj PRIVATE ${dummy_file})
endif ()

# Info for the main project
set (${util}_src     ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)
set (${util}_sources ${sources}                  PARENT_SCOPE)
set (${util}_mods    ${modfile_list}             PARENT_SCOPE)
