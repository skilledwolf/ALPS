#  Copyright Synge Todo 2011.
#   SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 4.3)
project(alps_fortran_examples LANGUAGES C CXX Fortran)
find_package(ALPS CONFIG REQUIRED COMPONENTS fortran)
option(ALPS_BUILD_TESTING "Build and register ALPS tests" ON)
if(ALPS_BUILD_TESTING)
  enable_testing()
endif()
include("${CMAKE_CURRENT_LIST_DIR}/../../../cmake/ALPSTesting.cmake")
add_subdirectory(hello)
add_subdirectory(ising)
