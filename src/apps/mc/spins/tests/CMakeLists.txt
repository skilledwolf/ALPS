# SPDX-License-Identifier: MIT

add_executable(spinmc_characterization characterization.cpp)
target_include_directories(spinmc_characterization PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(spinmc_characterization PRIVATE spinmc_impl)
alps_add_gtest(spinmc_characterization COMPONENT spinmc LABELS scientific TIMEOUT 120)

add_executable(spinmc_settings settings.cpp)
target_include_directories(spinmc_settings PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(spinmc_settings PRIVATE ALPS::alps)
alps_add_gtest(spinmc_settings COMPONENT spinmc)

add_executable(spinmc_couplings couplings.cpp)
target_include_directories(spinmc_couplings PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(spinmc_couplings PRIVATE ALPS::alps)
alps_add_gtest(spinmc_couplings COMPONENT spinmc)

add_executable(spinmc_model_choice model_choice.cpp)
target_include_directories(spinmc_model_choice PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(spinmc_model_choice PRIVATE ALPS::alps)
alps_add_gtest(spinmc_model_choice COMPONENT spinmc)

add_executable(spinmc_sampler sampler.cpp)
target_include_directories(spinmc_sampler PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(spinmc_sampler PRIVATE spinmc_impl)
alps_add_gtest(spinmc_sampler COMPONENT spinmc LABELS scientific TIMEOUT 120)
