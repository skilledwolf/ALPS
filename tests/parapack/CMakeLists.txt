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
add_executable(exchange_checkpoint exchange_checkpoint.cpp)
target_link_libraries(exchange_checkpoint PRIVATE ALPS::alps)
if(ALPS_ENABLE_MPI)
  add_test(NAME exchange_checkpoint
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:exchange_checkpoint> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(exchange_checkpoint PROPERTIES LABELS "parapack;mpi" PROCESSORS 2 TIMEOUT 60)
  add_test(NAME exchange_checkpoint_grouped
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:exchange_checkpoint> 2 ${MPIEXEC_POSTFLAGS})
  set_tests_properties(exchange_checkpoint_grouped PROPERTIES LABELS "parapack;mpi" PROCESSORS 2 TIMEOUT 60)
else()
  add_test(NAME exchange_checkpoint COMMAND exchange_checkpoint)
  set_property(TEST exchange_checkpoint PROPERTY LABELS parapack)
endif()
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

add_executable(task_state task_state.cpp)
target_link_libraries(task_state PRIVATE ALPS::alps)
add_test(NAME task_state COMMAND task_state)
set_property(TEST task_state PROPERTY LABELS parapack)

add_executable(clone_lifecycle clone_lifecycle.cpp)
target_link_libraries(clone_lifecycle PRIVATE ALPS::alps)
add_test(NAME clone_lifecycle COMMAND clone_lifecycle)
set_tests_properties(clone_lifecycle PROPERTIES LABELS parapack TIMEOUT 60)
if(ALPS_ENABLE_MPI)
  add_test(NAME clone_lifecycle_mpi
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:clone_lifecycle> mpi ${MPIEXEC_POSTFLAGS})
  set_tests_properties(clone_lifecycle_mpi PROPERTIES LABELS "parapack;mpi" PROCESSORS 2 TIMEOUT 60)
endif()
