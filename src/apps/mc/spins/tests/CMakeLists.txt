# SPDX-License-Identifier: MIT

add_executable(spinmc_characterization characterization.cpp)
target_include_directories(spinmc_characterization PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(spinmc_characterization PRIVATE spinmc_impl)
alps_add_gtest(spinmc_characterization COMPONENT spinmc LABELS scientific TIMEOUT 120)
