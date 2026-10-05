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

function(alps_native_qmc target)
  find_package(tomlplusplus 3.4 CONFIG REQUIRED)
  file(READ "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/schema.toml" QMC_COMMON_SCHEMA)
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml" QMC_APPLICATION_SCHEMA)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/schema.toml" "${CMAKE_CURRENT_SOURCE_DIR}/schema.toml")
  configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/schema.hpp.in" schema.hpp @ONLY)
  target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}"
    "$<TARGET_PROPERTY:tomlplusplus::tomlplusplus,INTERFACE_INCLUDE_DIRECTORIES>")
  target_compile_definitions(${target} PRIVATE TOML_HEADER_ONLY=1)
  target_link_libraries(${target} PRIVATE ALPS::run_config ALPS::statistics)
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${target}.toml"
    "application = \"${target}\"\n${QMC_COMMON_SCHEMA}\n${QMC_APPLICATION_SCHEMA}")
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${target}.toml"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas COMPONENT applications)
endfunction()

add_subdirectory(checksign)
add_subdirectory(looper)
add_subdirectory(qwl)
add_subdirectory(sse4)
add_subdirectory(worms)

if(ALPS_BUILD_TESTING)
  add_executable(qmc_measurement_contract tests/measurements.cpp)
  target_link_libraries(qmc_measurement_contract PRIVATE ALPS::alps ALPS::statistics)
  alps_add_test(qmc_measurement_contract)
  set_tests_properties(qmc_measurement_contract PROPERTIES LABELS "qmc;statistics"
    ENVIRONMENT_MODIFICATION "ALPS_XML_PATH=set:${PROJECT_SOURCE_DIR}/src/alps/resources")
endif()
