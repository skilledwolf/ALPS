#  Copyright Matthias Troyer 2009 - 2010.
#   SPDX-License-Identifier: MIT

add_executable(dmrg dmrg.C factory.C)
target_include_directories(dmrg PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(dmrg PRIVATE ALPS::alps)
alps_install_applications(dmrg)

if(ALPS_BUILD_TESTING)
  add_executable(dmrg_product_test product_test.cpp)
  target_include_directories(dmrg_product_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
  target_link_libraries(dmrg_product_test PRIVATE ALPS::alps)
  add_test(NAME dmrg_product COMMAND dmrg_product_test)
  set_property(TEST dmrg_product PROPERTY LABELS dmrg)
endif()
