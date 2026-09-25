#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_executable(dmrg dmrg.C factory.C)
target_include_directories(dmrg PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(dmrg PRIVATE ALPS::alps)
alps_install_applications(dmrg)

if(ALPS_BUILD_TESTING)
  foreach(test product state_slice system_layout block_matrix seed)
    add_executable(dmrg_${test}_test ${test}_test.cpp)
    target_include_directories(dmrg_${test}_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(dmrg_${test}_test PRIVATE ALPS::alps)
    add_test(NAME dmrg_${test} COMMAND dmrg_${test}_test)
    set_property(TEST dmrg_${test} PROPERTY LABELS dmrg)
  endforeach()
endif()
