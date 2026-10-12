#  Copyright Matthias Troyer 2009 - 2010.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

add_library(spinmc_impl STATIC abstract_fitter.C base_incr_fitter.C dummy_fitter.C
                               est_grad_fitter.C fitting_scheduler.C spinmc_factory.C)
target_link_libraries(spinmc_impl PUBLIC ALPS::alps)
add_executable(spinmc spinmc.C)
add_executable(spinmc_evaluate spinmc_evaluate.C)
target_link_libraries(spinmc PRIVATE spinmc_impl)
target_link_libraries(spinmc_evaluate PRIVATE spinmc_impl)
alps_install_applications(spinmc spinmc_evaluate)

if(ALPS_BUILD_TESTING)
  add_executable(spinmc_scientific_corrections tests/scientific_corrections.cpp)
  target_include_directories(spinmc_scientific_corrections PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_definitions(spinmc_scientific_corrections PRIVATE
    ALPS_SPINMC_LATTICES="${PROJECT_SOURCE_DIR}/src/alps/resources/lattices.xml")
  target_link_libraries(spinmc_scientific_corrections PRIVATE spinmc_impl)
  alps_add_gtest(spinmc_scientific_corrections COMPONENT applications LABELS regression scientific TIMEOUT 180)
endif()

if(ALPS_BUILD_TESTING)
  add_executable(spinmc_legacy_measurements tests/legacy_measurements.cpp)
  target_include_directories(spinmc_legacy_measurements PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_definitions(spinmc_legacy_measurements PRIVATE ALPS_SPINMC_LATTICES="${PROJECT_SOURCE_DIR}/src/alps/resources/lattices.xml")
  target_link_libraries(spinmc_legacy_measurements PRIVATE spinmc_impl)
  alps_add_gtest(spinmc_legacy_measurements COMPONENT applications LABELS regression scientific)
endif()
