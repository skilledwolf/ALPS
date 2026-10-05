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

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml" ISING_BASE_SCHEMA)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/lattice-schema.toml" ISING_LATTICE_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml" "${CMAKE_CURRENT_SOURCE_DIR}/lattice-schema.toml")
foreach(variant RANGE 1 3)
  set(ISING_SCHEMA "${ISING_BASE_SCHEMA}")
  if(variant EQUAL 2)
    string(APPEND ISING_SCHEMA "${ISING_LATTICE_SCHEMA}")
  endif()
  configure_file(ising_schema.hpp.in ${variant}/ising_schema.hpp @ONLY)
  add_executable(ising${variant} main.C)
  target_link_libraries(ising${variant} PRIVATE ALPS::alps)
  target_include_directories(ising${variant} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/${variant}")
  target_compile_definitions(ising${variant} PRIVATE ISING_VARIANT=${variant} ISING_COMMAND="ising${variant}")
endforeach()

foreach(command IN ITEMS evaluate evaluate2)
  add_executable(${command} evaluate.C)
  target_link_libraries(${command} PRIVATE ALPS::statistics ALPS::hdf5 ALPS::run_config)
endforeach()
target_compile_definitions(evaluate PRIVATE BINDER_FACTOR=1)
target_compile_definitions(evaluate2 PRIVATE BINDER_FACTOR=3)
