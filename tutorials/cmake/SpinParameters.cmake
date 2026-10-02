# Embed the application schema; the executable needs no source-tree files.
file(READ "${CMAKE_CURRENT_LIST_DIR}/../spin-schema.toml" SPIN_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_LIST_DIR}/../spin-schema.toml")
configure_file("${CMAKE_CURRENT_LIST_DIR}/spin_config.hpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/spin_config.hpp" @ONLY)
include_directories("${CMAKE_CURRENT_BINARY_DIR}")
