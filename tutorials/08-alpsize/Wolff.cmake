# The two lessons share physics/statistics, while demonstrating different lattices.
file(READ "${CMAKE_CURRENT_LIST_DIR}/wolff-schema.toml" WOLFF_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_LIST_DIR}/wolff-schema.toml")
string(APPEND WOLFF_SCHEMA "${WOLFF_LATTICE_SCHEMA}")
configure_file("${CMAKE_CURRENT_LIST_DIR}/wolff_schema.hpp.in" wolff_schema.hpp @ONLY)
target_include_directories(wolff PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(wolff PRIVATE ALPS::statistics ALPS::hdf5 ALPS::run_config)
