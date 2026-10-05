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

find_package(tomlplusplus 3.4 CONFIG REQUIRED)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml" SPINMC_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml")
configure_file(schema.hpp.in schema.hpp @ONLY)
add_executable(spinmc spinmc.C)
target_include_directories(spinmc PRIVATE "${CMAKE_CURRENT_BINARY_DIR}"
  "$<TARGET_PROPERTY:tomlplusplus::tomlplusplus,INTERFACE_INCLUDE_DIRECTORIES>")
target_compile_definitions(spinmc PRIVATE TOML_HEADER_ONLY=1)
target_link_libraries(spinmc PRIVATE ALPS::alps ALPS::run_config ALPS::statistics)
alps_install_applications(spinmc)
install(FILES schema.toml DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas RENAME spinmc.toml COMPONENT applications)

if(ALPS_BUILD_TESTING)
  add_executable(spinmc_native_contract tests/native.cpp)
  target_link_libraries(spinmc_native_contract PRIVATE ALPS::alps ALPS::statistics)
  alps_add_test(spinmc_native_contract)
  set_tests_properties(spinmc_native_contract PROPERTIES LABELS "spinmc;statistics;checkpoint" TIMEOUT 180)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter)
  add_test(NAME spinmc_cli_contract
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_contract.py" "$<TARGET_FILE:spinmc>")
  set_tests_properties(spinmc_cli_contract PROPERTIES LABELS "spinmc;cli;run_config" TIMEOUT 60
    ENVIRONMENT_MODIFICATION "ALPS_XML_PATH=set:${PROJECT_SOURCE_DIR}/src/alps/resources")
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      set_property(TEST spinmc_cli_contract APPEND PROPERTY ENVIRONMENT_MODIFICATION
        "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
endif()
