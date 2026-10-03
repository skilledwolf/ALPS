#  Copyright Matthias Troyer 2009 - 2010
#            Lukas Gamper 2012
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

foreach(name contract checkpoint)
  add_executable(param_${name} ${name}.cpp)

  target_link_libraries(param_${name} PRIVATE ALPS::params)
  alps_add_test(param_${name})
  set_property(TEST param_${name} PROPERTY LABELS params)
endforeach()

if(ALPS_ENABLE_MPI)
  add_executable(param_mpi mpi.cpp)
  target_link_libraries(param_mpi PRIVATE ALPS::params)
  add_test(NAME param_mpi
    COMMAND "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 2
      ${MPIEXEC_PREFLAGS} "$<TARGET_FILE:param_mpi>" ${MPIEXEC_POSTFLAGS})
  set_tests_properties(param_mpi PROPERTIES LABELS "params;mpi" PROCESSORS 2 TIMEOUT 60)
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      set_property(TEST param_mpi APPEND PROPERTY ENVIRONMENT_MODIFICATION
        "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
endif()
