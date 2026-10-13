# SPDX-License-Identifier: MIT
add_executable(qmc_worm_characterization worm_characterization.cpp)
target_link_libraries(qmc_worm_characterization PRIVATE worm_impl)
alps_add_gtest(qmc_worm_characterization COMPONENT qmc LABELS characterization regression)

add_executable(qmc_sse_characterization sse_characterization.cpp ../sse4/lp_sse.cpp)
target_link_libraries(qmc_sse_characterization PRIVATE ALPS::alps)
alps_add_gtest(qmc_sse_characterization COMPONENT qmc LABELS characterization regression)

add_executable(qmc_qwl_characterization qwl_characterization.cpp)
target_link_libraries(qmc_qwl_characterization PRIVATE ALPS::alps)
alps_add_gtest(qmc_qwl_characterization COMPONENT qmc LABELS characterization regression)

add_executable(qmc_loop_characterization loop_characterization.cpp)
target_link_libraries(qmc_loop_characterization PRIVATE ALPS::alps)
alps_add_gtest(qmc_loop_characterization COMPONENT qmc LABELS characterization unit)

foreach(solver IN ITEMS worm sse qwl)
  add_executable(qmc_${solver}_settings ${solver}_settings.cpp)
  target_link_libraries(qmc_${solver}_settings PRIVATE ALPS::alps)
  alps_add_gtest(qmc_${solver}_settings COMPONENT qmc LABELS unit)
endforeach()
