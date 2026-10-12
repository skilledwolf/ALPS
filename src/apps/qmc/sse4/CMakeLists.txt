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

add_executable(dirloop_sse main.cc lp_sse.cpp)
target_link_libraries(dirloop_sse PRIVATE ALPS::alps)
alps_install_applications(dirloop_sse)

if(ALPS_BUILD_TESTING)
  add_executable(sse_legacy_green tests/legacy_green.cpp lp_sse.cpp)
  target_include_directories(sse_legacy_green PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_definitions(sse_legacy_green PRIVATE
    ALPS_DIRLOOP_LATTICES="${PROJECT_SOURCE_DIR}/src/alps/resources/lattices.xml"
    ALPS_DIRLOOP_MODELS="${PROJECT_SOURCE_DIR}/src/alps/resources/models.xml")
  target_link_libraries(sse_legacy_green PRIVATE ALPS::alps BLAS::BLAS LAPACK::LAPACK)
  alps_add_gtest(sse_legacy_green COMPONENT applications LABELS regression scientific TIMEOUT 180)
endif()
