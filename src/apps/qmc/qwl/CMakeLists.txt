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
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml" QWL_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml")
configure_file(schema.hpp.in schema.hpp @ONLY)
add_executable(qwl qwl.C)
target_include_directories(qwl PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(qwl_evaluate qwl_evaluate.C)

target_link_libraries(qwl PRIVATE ALPS::alps ALPS::run_config ALPS::statistics)
target_link_libraries(qwl_evaluate PRIVATE ALPS::alps ALPS::statistics)
alps_install_applications(qwl qwl_evaluate)

install(FILES schema.toml DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas RENAME qwl.toml COMPONENT applications)

if(ALPS_BUILD_TESTING)
  add_executable(qwl_native_contract tests/native.cpp)
  target_link_libraries(qwl_native_contract PRIVATE ALPS::alps ALPS::statistics)
  alps_add_test(qwl_native_contract)
  set_tests_properties(qwl_native_contract PROPERTIES LABELS "qwl;statistics;checkpoint" TIMEOUT 120
    ENVIRONMENT_MODIFICATION "ALPS_XML_PATH=set:${PROJECT_SOURCE_DIR}/src/alps/resources")
endif()
