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

function(alps_diag_application target)
  set(schemas "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/schema.toml")
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml")
    list(APPEND schemas "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml")
  endif()
  set(DIAG_SCHEMA "application = \"${target}\"\n")
  foreach(schema IN LISTS schemas)
    file(READ "${schema}" text)
    string(APPEND DIAG_SCHEMA "${text}")
  endforeach()
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${schemas})
  configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/schema.hpp.in" schema.hpp @ONLY)
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${target}.toml" "${DIAG_SCHEMA}")

  add_executable(${target} ${target}.C)
  target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
  target_link_libraries(${target} PRIVATE ALPS::alps ALPS::run_config)
  alps_install_applications(${target})
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${target}.toml"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas COMPONENT applications)
endfunction()

add_subdirectory(sparsediag)
add_subdirectory(fulldiag)
