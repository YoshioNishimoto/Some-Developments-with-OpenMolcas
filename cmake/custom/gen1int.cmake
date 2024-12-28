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
# Copyright (C) 2017, Stefan Knecht                                    *
#               2021, Ignacio Fdez. Galván                             *
#***********************************************************************
#                                                                      *
#***********************************************************************
# CMakeLists.txt for GEN1INT in Molcas                                 *
#***********************************************************************

# load External Project macro
include(ExternalProject)
# Set up compilation of GEN1INT components
set(CUSTOM_GEN1INT_LOCATION ${PROJECT_BINARY_DIR}/External/gen1int)

# GEN1INT does not know profile
if(CMAKE_BUILD_TYPE MATCHES "profile")
  set(GEN1INT_BUILD_TYPE "release")
else()
  set(GEN1INT_BUILD_TYPE ${CMAKE_BUILD_TYPE})
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_ROOT})
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/custom)

if(SINGLE_MOD_DIR)
  set(mod_dir ${MAIN_MOD_DIR}/_single)
else()
  set(mod_dir ${MAIN_MOD_DIR}/gen1int_util)
endif()

list(APPEND GEN1INTCMakeArgs
  -DCMAKE_BUILD_TYPE=${GEN1INT_BUILD_TYPE}
  -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/External
  -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS}
  -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
  -DCMAKE_Fortran_MODULE_DIRECTORY=${mod_dir}
  -DEXTRA_INCLUDE=${OPENMOLCAS_DIR}/src/Include
  )

#####################################
# git references for GEN1INT module #
#####################################
set(reference_git_repo https://gitlab.com/Molcas/gen1int-molcas.git)
set(reference_git_commit 75353eea270b4cccb746efbb919e7f96dd605d7c)
set(EP_PROJECT gen1int)

# Enabling source changes to keep ExternalProject happy
set (CMAKE_DISABLE_SOURCE_CHANGES OFF)

set (last_hash "None")
set (hash_file ${CUSTOM_GEN1INT_LOCATION}/${EP_PROJECT}.hash)
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
                    PREFIX ${CUSTOM_GEN1INT_LOCATION}
                    CMAKE_ARGS "${GEN1INTCMakeArgs}"
                    GIT_REPOSITORY ${reference_git_repo}
                    GIT_TAG ${reference_git_commit}
                    UPDATE_DISCONNECTED ${EP_SkipUpdate}
                    INSTALL_DIR "${PROJECT_BINARY_DIR}"
                   )

ExternalProject_Add_Step (${EP_PROJECT} update_hash
                          COMMAND echo ${reference_git_commit} > ${hash_file}
                          DEPENDEES build
                         )

set (CMAKE_DISABLE_SOURCE_CHANGES ON)

# set variables for use in parent CMakeLists.txt
ExternalProject_Get_Property(${EP_PROJECT} install_dir)
ExternalProject_Get_Property(${EP_PROJECT} source_dir)
set(GEN1INT_INCLUDE ${mod_dir} ${source_dir}/src PARENT_SCOPE)
set(GEN1INT_LIBRARIES ${install_dir}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}gen1int-molcaslib.a PARENT_SCOPE)
