add_subdirectory(ngs)

add_executable(scheduler_seeds_characterization seeds_characterization.cpp)
target_link_libraries(scheduler_seeds_characterization PRIVATE ALPS::alps)
alps_add_gtest(scheduler_seeds_characterization COMPONENT scheduler)

add_executable(scheduler_task_characterization task_characterization.cpp)
target_link_libraries(scheduler_task_characterization PRIVATE ALPS::alps)
alps_add_gtest(scheduler_task_characterization COMPONENT scheduler)

# seed_vectors.txt is shared with the check of pyalps/seeds.py below.
add_executable(scheduler_seeds seeds.cpp)
target_compile_definitions(scheduler_seeds PRIVATE
  ALPS_SEED_VECTORS="${CMAKE_CURRENT_SOURCE_DIR}/seed_vectors.txt")
target_link_libraries(scheduler_seeds PRIVATE ALPS::alps)
alps_add_gtest(scheduler_seeds COMPONENT scheduler)

find_package(Python 3.11 QUIET COMPONENTS Interpreter)
if(Python_Interpreter_FOUND)
  add_test(NAME pyalps_task_seeds
    COMMAND "${Python_EXECUTABLE}" -I "${CMAKE_CURRENT_SOURCE_DIR}/check_task_seeds.py"
      "${PROJECT_SOURCE_DIR}/python/pyalps/src/pyalps/seeds.py"
      "${CMAKE_CURRENT_SOURCE_DIR}/seed_vectors.txt")
  set_tests_properties(pyalps_task_seeds PROPERTIES TIMEOUT 60 LABELS "scheduler;unit")
endif()
