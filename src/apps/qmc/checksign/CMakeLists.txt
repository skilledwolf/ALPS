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

add_executable(checksign checksign.C)
target_include_directories(checksign PRIVATE "${PROJECT_SOURCE_DIR}/src/apps/lattice_model")
target_link_libraries(checksign PRIVATE ALPS::alps ALPS::run_config)
alps_install_applications(checksign)

if(ALPS_BUILD_TESTING)
  add_test(NAME checksign COMMAND checksign unfrustrated.toml frustrated.toml
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  set_tests_properties(checksign PROPERTIES LABELS "qmc;cli;run_config" TIMEOUT 60
    PASS_REGULAR_EXPRESSION "unfrustrated.toml: OK\nfrustrated.toml: SIGN PROBLEM\n"
    ENVIRONMENT_MODIFICATION "ALPS_XML_PATH=set:${PROJECT_SOURCE_DIR}/src/alps/resources")
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      set_property(TEST checksign APPEND PROPERTY ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
endif()
