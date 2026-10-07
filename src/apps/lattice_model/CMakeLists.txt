# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
# alps_lattice_model_application(<target> [SCHEMAS <file>...])
# Builds <target> from <target>.C with the run schema of the lattice and model
# libraries, extended by the listed schema files, and installs both.
function(alps_lattice_model_application target)
  cmake_parse_arguments(PARSE_ARGV 1 ARG "" "" "SCHEMAS")
  set(schemas "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/schema.toml")
  foreach(schema IN LISTS ARG_SCHEMAS)
    cmake_path(ABSOLUTE_PATH schema NORMALIZE)
    list(APPEND schemas "${schema}")
  endforeach()
  set(RUN_SCHEMA "application = \"${target}\"\n")
  foreach(schema IN LISTS schemas)
    file(READ "${schema}" text)
    string(APPEND RUN_SCHEMA "${text}")
  endforeach()
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${schemas})
  configure_file("${CMAKE_CURRENT_FUNCTION_LIST_DIR}/schema.hpp.in" schema.hpp @ONLY)
  file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/${target}.toml" "${RUN_SCHEMA}")

  add_executable(${target} ${target}.C)
  target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_FUNCTION_LIST_DIR}")
  target_link_libraries(${target} PRIVATE ALPS::alps ALPS::run_config)
  alps_install_applications(${target})
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${target}.toml"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas COMPONENT applications)
endfunction()
