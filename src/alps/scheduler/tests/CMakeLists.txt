add_subdirectory(ngs)

add_executable(scheduler_checkpoint checkpoint.cpp)
target_link_libraries(scheduler_checkpoint PRIVATE ALPS::alps)
alps_add_test(scheduler_checkpoint)
set_tests_properties(scheduler_checkpoint PROPERTIES LABELS "scheduler;hdf5")
