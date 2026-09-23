# Copyright Matthias Troyer and Synge Todo 2009-2013.
# SPDX-License-Identifier: MIT

alps_add_test(lattice2xml)

add_test(NAME parameter_conversion
  COMMAND ${CMAKE_COMMAND} -Dxml_program=$<TARGET_FILE:parameter2xml>
    -Dhdf5_program=$<TARGET_FILE:parameter2hdf5>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/parameter_conversion.cmake)
set_property(TEST parameter_conversion PROPERTY LABELS cli)
