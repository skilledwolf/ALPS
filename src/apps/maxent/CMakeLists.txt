# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema/maxent.toml" MAXENT_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/schema/maxent.toml")
configure_file(src/schema.hpp.in maxent_schema.hpp @ONLY)

add_library(alps_maxent STATIC
  src/input.cpp src/maxent_run.cpp src/maxent_helper.cpp src/maxent_simulation.cpp src/maxent_parms.cpp)
add_library(ALPS::maxent ALIAS alps_maxent)
set_target_properties(alps_maxent PROPERTIES EXPORT_NAME maxent
  CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
target_include_directories(alps_maxent PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(alps_maxent PUBLIC ALPS::run_config ALPS::solver_headers ALPS::params ALPS::hdf5 ALPS::utilities ALPS::osiris
  Boost::filesystem Boost::date_time ${ALPS_BLAS_TARGET} ${ALPS_LAPACK_TARGET})
install(TARGETS alps_maxent EXPORT ALPSTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)

add_executable(maxent cli/maxent.cpp)
target_link_libraries(maxent PRIVATE ALPS::maxent)
alps_register_module(NAME maxent ROOT "${CMAKE_CURRENT_SOURCE_DIR}/src"
  TARGETS alps_maxent
  PRIVATE_FILES src/default_model.hpp src/maxent.hpp src/maxent_parms.hpp
  DEPENDS configuration hdf5 osiris params run_config solvers)
alps_register_module(NAME maxent_cli ROOT "${CMAKE_CURRENT_SOURCE_DIR}/cli"
  TARGETS maxent DEPENDS hdf5 params run_config solvers)
install(TARGETS maxent RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

install(FILES schema/maxent.toml DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas COMPONENT libraries)
