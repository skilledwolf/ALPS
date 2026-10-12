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
add_library(worm_impl STATIC WModel.C WRun.C Wcheck.C Wdostep.C Winit.C Wmeas.C)
target_link_libraries(worm_impl PUBLIC ALPS::alps)

add_executable(worm main.C)
add_executable(worm_evaluate evaluate.C)

target_link_libraries(worm PRIVATE worm_impl)
target_link_libraries(worm_evaluate PRIVATE worm_impl)
alps_install_applications(worm worm_evaluate)

if(ALPS_BUILD_TESTING)
  add_executable(worm_time_endpoints tests/time_endpoints.cpp)
  target_include_directories(worm_time_endpoints PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(worm_time_endpoints PRIVATE worm_impl)
  alps_add_gtest(worm_time_endpoints COMPONENT applications LABELS regression scientific)
  add_executable(worm_initial_state tests/initial_state.cpp)
  target_include_directories(worm_initial_state PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_definitions(worm_initial_state PRIVATE
    ALPS_WORM_LATTICES="${PROJECT_SOURCE_DIR}/src/alps/resources/lattices.xml"
    ALPS_WORM_MODELS="${PROJECT_SOURCE_DIR}/src/alps/resources/models.xml")
  target_link_libraries(worm_initial_state PRIVATE worm_impl)
  alps_add_gtest(worm_initial_state COMPONENT applications LABELS regression scientific)
endif()

if(ALPS_BUILD_TESTING)
  add_executable(worm_legacy_geometry tests/legacy_geometry.cpp)
  target_include_directories(worm_legacy_geometry PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_definitions(worm_legacy_geometry PRIVATE
    ALPS_WORM_LATTICES="${PROJECT_SOURCE_DIR}/src/alps/resources/lattices.xml"
    ALPS_WORM_MODELS="${PROJECT_SOURCE_DIR}/src/alps/resources/models.xml")
  target_link_libraries(worm_legacy_geometry PRIVATE worm_impl BLAS::BLAS LAPACK::LAPACK)
  alps_add_gtest(worm_legacy_geometry COMPONENT applications LABELS regression scientific TIMEOUT 180)
endif()
