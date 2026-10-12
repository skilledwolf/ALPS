# SPDX-License-Identifier: MIT
foreach(example IN ITEMS ngs/1_accumulator_only ngs/2_single_core ngs/3_mpi
    ngs/4_mpi_pscan code-06-mcmain-c++ code-07-mcmain-mcbase)
  string(MAKE_C_IDENTIFIER "ising_${example}" target)
  set(directory "${PROJECT_SOURCE_DIR}/tutorials/${example}")
  add_executable(${target} ising_energy.cpp "${directory}/ising.cpp")
  target_include_directories(${target} PRIVATE "${directory}")
  target_link_libraries(${target} PRIVATE ALPS::alps)
  alps_add_gtest(${target} COMPONENT tutorials LABELS scientific)
endforeach()
