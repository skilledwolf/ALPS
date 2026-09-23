#  Copyright Matthias Troyer and Synge Todo 2010-2011.
#   SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 4.3)
project(alps_examples LANGUAGES C CXX)
find_package(ALPS CONFIG REQUIRED)
set(ALPS_ENABLE_MPI "${ALPS_MPI_ENABLED}")
option(ALPS_BUILD_TESTING "Build and register ALPS tests" ON)
if(ALPS_BUILD_TESTING)
  enable_testing()
endif()
include("${CMAKE_CURRENT_LIST_DIR}/../../cmake/ALPSTesting.cmake")

add_subdirectory(alea)
add_subdirectory(parapack)
add_subdirectory(scheduler)
add_subdirectory(sampling)
add_subdirectory(model)
add_subdirectory(ietl)
add_subdirectory(hdf5)
