# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

alps_add_test(lattice2xml)

add_test(NAME parameter_conversion
  COMMAND ${CMAKE_COMMAND} -Dxml_program=$<TARGET_FILE:parameter2xml>
    -Dhdf5_program=$<TARGET_FILE:parameter2hdf5>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/parameter_conversion.cmake)
set_property(TEST parameter_conversion PROPERTY LABELS cli)

add_executable(spin_simulation spin_simulation.cpp)
target_link_libraries(spin_simulation PRIVATE spinmc_impl)
target_include_directories(spin_simulation PRIVATE ${PROJECT_SOURCE_DIR}/src/apps/mc/spins)
target_compile_definitions(spin_simulation PRIVATE
  ALPS_TEST_LATTICES="${PROJECT_SOURCE_DIR}/src/alps/resources/lattices.xml")
add_test(NAME spin_simulation COMMAND spin_simulation)
set_tests_properties(spin_simulation PROPERTIES LABELS spinmc TIMEOUT 60)
