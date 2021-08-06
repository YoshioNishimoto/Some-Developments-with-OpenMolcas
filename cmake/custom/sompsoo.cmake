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
# Copyright (C) 2021, Stefan Knecht                                    *
#***********************************************************************
#                                                                      *
# Based on NEVPT2 Makefile by Leon Freitag                             *
#                                                                      *
#***********************************************************************
# CMakeLists.txt for soMPSoo                                           *
#***********************************************************************

# load External Project macro
include(ExternalProject)
# Set up compilation of soMPSoo components
set(CUSTOM_soMPSoo_LOCATION ${PROJECT_BINARY_DIR}/External/soMPSoo)

# soMPSoo (and QCMaquis) do not know profile
if(CMAKE_BUILD_TYPE MATCHES "profile")
  set(soMPSoo_BUILD_TYPE "release")
else()
  set(soMPSoo_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_ROOT})
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/custom)

if (LINALG_LIBRARIES)
  target_files(LINALG_LIBRARIES_FILES ${LINALG_LIBRARIES})
elseif (LINALG STREQUAL "Internal")
  set (LINALG_LIBRARIES_FILES $<TARGET_FILE:blas> $<TARGET_FILE:lapack>)
endif()

# CMake does not support passing lists inside lists, so we need to
# replace the semicolons in the lists and pass them as normal strings
# and then replace the new separators with semicolons on the other side

message(status "flags: ${CMAKE_Fortran_FLAGS}")

list(APPEND soMPSooCMakeArgs
  "-DENABLE_STANDALONE:BOOL=OFF"
  "-DHOST_PROGRAM:STRING=OpenMolcas"
  "-DCMAKE_BUILD_TYPE=${soMPSoo_BUILD_TYPE}"
  "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/External"
  "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}"
  "-DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}"
  "-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
  "-DCMAKE_C_FLAGS=${CMake_C_FLAGS}"
  "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
  "-DCMAKE_CXX_FLAGS=${CMake_CXX_FLAGS}"
  "-DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>"
  "-DENABLE_OPENMP:BOOL=OFF"
  "-DENABLE_BLAS:STRING=auto"
  "-DENABLE_LAPACK:STRING=auto"
  )
if(LINALG STREQUAL "Accelerate")
  set(HasVeclib "-DENABLE_ACCELERATE:BOOL=ON")
else()
  set(HasVeclib "-DENABLE_ACCELERATE:BOOL=OFF")
endif()
list(APPEND soMPSooCMakeArgs hasVeclib)

######################################
# git references for soMPSoo         #
######################################
set(reference_git_repo https://github.com/soMPSoo/beta-release.git)
set(reference_git_commit 0d59196)


set(EP_PROJECT soMPSoo)

# Enabling source changes to keep ExternalProject happy
set (CMAKE_DISABLE_SOURCE_CHANGES OFF)

set (last_hash "None")
set (hash_file ${CUSTOM_soMPSoo_LOCATION}/${EP_PROJECT}.hash)
if (EXISTS ${hash_file})
  file (READ ${hash_file} last_hash)
  string (REGEX REPLACE "\n$" "" last_hash "${last_hash}")
endif ()
if (last_hash STREQUAL ${reference_git_commit})
  set (EP_SkipUpdate ON)
else ()
  set (EP_SkipUpdate OFF)
endif ()

ExternalProject_Add(${EP_PROJECT}
                    PREFIX ${CUSTOM_soMPSoo_LOCATION}
                    GIT_REPOSITORY ${reference_git_repo}
                    GIT_TAG ${reference_git_commit}
                    UPDATE_DISCONNECTED ${EP_SkipUpdate}
                    CMAKE_ARGS "${soMPSooCMakeArgs}"
                    INSTALL_DIR "${PROJECT_BINARY_DIR}/soMPSoo"
                   )

ExternalProject_Add_Step (${EP_PROJECT} update_hash
                          COMMAND echo ${reference_git_commit} > ${hash_file}
                          DEPENDEES build
                         )

set (CMAKE_DISABLE_SOURCE_CHANGES ON)

# set variables for use in parent CMakeLists.txt
ExternalProject_Get_Property(${EP_PROJECT} BINARY_DIR)

set(SOMPSOO_INCLUDE ${BINARY_DIR}/modules PARENT_SCOPE)
set(SOMPSOO_LIBRARIES ${BINARY_DIR}/src/${CMAKE_FIND_LIBRARY_PREFIXES}soMPSoo.a PARENT_SCOPE)
