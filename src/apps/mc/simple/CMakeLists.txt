#  Copyright Synge Todo 2010 - 2015.
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
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml" SIMPLEMC_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml")
configure_file(schema.hpp.in schema.hpp @ONLY)
add_executable(simplemc main.C)
target_include_directories(simplemc PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

target_link_libraries(simplemc PRIVATE ALPS::alps ALPS::run_config ALPS::statistics)
alps_install_applications(simplemc)
install(FILES schema.toml DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas RENAME simplemc.toml COMPONENT applications)

if(ALPS_BUILD_TESTING)
  add_executable(simplemc_native_contract tests/simulation.cpp)

  target_link_libraries(simplemc_native_contract PRIVATE ALPS::alps ALPS::statistics ALPS::run_config)
  alps_add_test(simplemc_native_contract)
  set_tests_properties(simplemc_native_contract PROPERTIES LABELS "simplemc;statistics;checkpoint" TIMEOUT 60)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter)
  add_test(NAME simplemc_cli_contract
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_contract.py" "$<TARGET_FILE:simplemc>")
  set_tests_properties(simplemc_cli_contract PROPERTIES LABELS "simplemc;cli;run_config" TIMEOUT 60
    ENVIRONMENT_MODIFICATION "ALPS_XML_PATH=set:${PROJECT_SOURCE_DIR}/src/alps/resources")
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      set_property(TEST simplemc_cli_contract APPEND PROPERTY ENVIRONMENT_MODIFICATION
        "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
endif()
