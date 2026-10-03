# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema/cthyb.toml" CTHYB_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/schema/cthyb.toml")
configure_file(schema.hpp.in cthyb_schema.hpp @ONLY)
target_sources(alps_cthyb PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/run_config.cpp")
target_include_directories(alps_cthyb PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(alps_cthyb PUBLIC ALPS::run_config)
install(FILES schema/cthyb.toml DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas COMPONENT libraries)

if(ALPS_BUILD_TESTING)
  add_executable(cthyb_run_config tests/run_config.cpp)
  target_link_libraries(cthyb_run_config PRIVATE ALPS::cthyb)
  alps_add_test(cthyb_run_config)
  set_property(TEST cthyb_run_config PROPERTY LABELS dmft)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter)
  add_test(NAME cthyb_cli_contract
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/tests/cli_contract.py"
            "$<TARGET_FILE:hybridization>")
  set_property(TEST cthyb_cli_contract PROPERTY LABELS dmft)
  set_property(TEST cthyb_cli_contract PROPERTY TIMEOUT 60)
  if(WIN32)
    foreach(directory IN LISTS ALPS_RUNTIME_LIBRARY_DIRS)
      set_property(TEST cthyb_cli_contract APPEND PROPERTY ENVIRONMENT_MODIFICATION
        "PATH=path_list_prepend:${directory}")
    endforeach()
  endif()
endif()
