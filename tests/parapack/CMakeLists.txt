#  Copyright Matthias Troyer and Synge Todo 2009 - 2010.
#   SPDX-License-Identifier: MIT

foreach(
  name
  clone_info
  clone_phase
  clone_timer
  exmc_optimize
  exp_number
  footprint
  info_test
  integer_range
  integer_range_limits
  linear_regression
  merge
  percentage
  temperature_scan
  time
  version
  wl_weight)
  add_executable(${name} ${name}.C)

  target_link_libraries(${name} PRIVATE ALPS::alps)
  alps_add_test(${name})
  set_property(TEST ${name} PROPERTY LABELS parapack)
endforeach()

add_executable(scheduler_workflow scheduler_workflow.cpp)
target_link_libraries(scheduler_workflow PRIVATE ALPS::alps)
foreach(runtime classic ngs)
  add_test(NAME scheduler_workflow_${runtime}
    COMMAND ${CMAKE_COMMAND} -Dprogram=$<TARGET_FILE:scheduler_workflow>
      -Druntime=${runtime} -P ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_workflow.cmake)
  set_property(TEST scheduler_workflow_${runtime} PROPERTY LABELS parapack)
  if(ALPS_ENABLE_MPI)
    add_test(NAME scheduler_workflow_${runtime}_mpi
      COMMAND ${CMAKE_COMMAND} -Dprogram=$<TARGET_FILE:scheduler_workflow>
        -Druntime=${runtime}
        "-Dlauncher=${MPIEXEC_EXECUTABLE};${MPIEXEC_NUMPROC_FLAG};2;${MPIEXEC_PREFLAGS}"
        "-Dscheduler_options=--mpi;--total-threads;2"
        -P ${CMAKE_CURRENT_SOURCE_DIR}/scheduler_workflow.cmake)
    set_tests_properties(scheduler_workflow_${runtime}_mpi PROPERTIES
      LABELS "parapack;mpi" PROCESSORS 2 TIMEOUT 120)
  endif()
endforeach()

if(ALPS_ENABLE_MPI)
  foreach(name collect_mpi comm_mpi filelock_mpi halt_mpi info_test_mpi process_mpi)
    add_executable(${name} ${name}.C)

    target_link_libraries(${name} PRIVATE ALPS::alps)
    # alps_add_test(${name})
    # set_property(TEST ${name} PROPERTY LABELS parapack)
  endforeach()
endif()
