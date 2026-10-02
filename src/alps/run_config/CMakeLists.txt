# Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
find_package(tomlplusplus 3.4 CONFIG REQUIRED)
add_library(alps_run_config src/run_config.cpp)
add_library(ALPS::run_config ALIAS alps_run_config)
set_target_properties(alps_run_config PROPERTIES EXPORT_NAME run_config
  VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR}
  CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON
  INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
target_link_libraries(alps_run_config PUBLIC ALPS::params)
# Distribution packages may expose a compiled parser. Use its headers in
# header-only mode so static SDK consumers need no toml++ target or runtime.
target_include_directories(alps_run_config SYSTEM PRIVATE
  "$<TARGET_PROPERTY:tomlplusplus::tomlplusplus,INTERFACE_INCLUDE_DIRECTORIES>")
target_compile_definitions(alps_run_config PRIVATE TOML_HEADER_ONLY=1)
generate_export_header(alps_run_config EXPORT_MACRO_NAME ALPS_RUN_CONFIG_DECL
  EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/generated/include/alps/run_config_export.h")
if(NOT BUILD_SHARED_LIBS)
  target_compile_definitions(alps_run_config PUBLIC ALPS_RUN_CONFIG_STATIC_DEFINE)
endif()
target_sources(alps_run_config PUBLIC FILE_SET HEADERS
  BASE_DIRS "${PROJECT_BINARY_DIR}/generated/include"
  FILES "${PROJECT_BINARY_DIR}/generated/include/alps/run_config_export.h")
target_sources(alps_run_config PUBLIC FILE_SET run_config_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include" FILES include/alps/run_config.hpp)
alps_dependency_rpath(alps_run_config)
alps_copy_runtime_dependencies(alps_run_config)
alps_register_module(NAME run_config HEADER_SET run_config_headers
  HEADER_TARGETS alps_run_config TARGETS alps_run_config
  DEPENDS params hdf5 GENERATED_HEADERS alps/run_config_export.h)
if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()
