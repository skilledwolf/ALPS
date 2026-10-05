# Mc source ownership.
target_sources(alps PRIVATE
  src/mcbase.cpp
  src/ngs/api.cpp
  src/stop_callback.cpp)
target_link_libraries(alps PUBLIC ALPS::statistics ALPS::run_config)
target_link_libraries(alps_headers INTERFACE "$<COMPILE_ONLY:ALPS::statistics;ALPS::run_config>")

target_sources(alps_headers INTERFACE
  FILE_SET mc_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/check_schedule.hpp
    include/alps/mcbase.hpp
    include/alps/mc/driver.hpp
    include/alps/mc/measurements.hpp
    include/alps/mc/physical_moments.hpp
    include/alps/mc/parallel.hpp
    include/alps/mcmpiadapter.hpp
    include/alps/ngs.hpp
    include/alps/ngs/api.hpp
    include/alps/progress_callback.hpp
    include/alps/stop_callback.hpp)

alps_register_module(NAME mc HEADER_SET mc_headers
  TARGETS alps
  DEPENDS alea cli configuration hdf5 params random run_config utilities xml
)

if(ALPS_BUILD_TESTING)
  add_executable(mcbase_checkpoint tests/checkpoint.cpp)
  target_link_libraries(mcbase_checkpoint PRIVATE ALPS::alps)
  alps_add_test(mcbase_checkpoint)
  set_property(TEST mcbase_checkpoint PROPERTY LABELS "mc;hdf5")

  add_executable(mcbase_ising_checkpoint
    "${PROJECT_SOURCE_DIR}/tutorials/10-ngs/2_single_core/checkpoint.cpp"
    "${PROJECT_SOURCE_DIR}/tutorials/09-code/07-mcmain-mcbase/ising.cpp")
  target_include_directories(mcbase_ising_checkpoint PRIVATE
    "${PROJECT_SOURCE_DIR}/tutorials/09-code/07-mcmain-mcbase")
  target_link_libraries(mcbase_ising_checkpoint PRIVATE ALPS::alps)
  alps_add_test(mcbase_ising_checkpoint)
  set_property(TEST mcbase_ising_checkpoint PROPERTY LABELS "mc;alea;hdf5;tutorial")

  add_executable(mcbase_heisenberg_boltzmann
    "${PROJECT_SOURCE_DIR}/tutorials/09-code/07-mcmain-mcbase/heisenberg/1d_lattice/two_spins.cpp"
    "${PROJECT_SOURCE_DIR}/tutorials/09-code/07-mcmain-mcbase/heisenberg/1d_lattice/heisenberg.cpp")
  target_include_directories(mcbase_heisenberg_boltzmann PRIVATE
    "${PROJECT_SOURCE_DIR}/tutorials/09-code/07-mcmain-mcbase/heisenberg/1d_lattice")
  target_link_libraries(mcbase_heisenberg_boltzmann PRIVATE ALPS::alps)
  alps_add_test(mcbase_heisenberg_boltzmann)
  set_property(TEST mcbase_heisenberg_boltzmann PROPERTY LABELS "mc;alea;tutorial")

  set(heisenberg_tutorial "${PROJECT_SOURCE_DIR}/tutorials/09-code/07-mcmain-mcbase/heisenberg")
  foreach(variant IN ITEMS 1d_lattice nd_lattice o_n_model)
    set(test mcbase_heisenberg_${variant}_checkpoint)
    add_executable(${test} "${heisenberg_tutorial}/checkpoint.cpp")
    if(variant STREQUAL "o_n_model")
      target_compile_definitions(${test} PRIVATE ALPS_NDIM_HEISENBERG)
    else()
      target_sources(${test} PRIVATE "${heisenberg_tutorial}/${variant}/heisenberg.cpp")
      target_include_directories(${test} PRIVATE "${heisenberg_tutorial}/${variant}")
      if(variant STREQUAL "nd_lattice")
        target_compile_definitions(${test} PRIVATE ALPS_LATTICE_HEISENBERG)
      endif()
    endif()
    target_link_libraries(${test} PRIVATE ALPS::alps)
    alps_add_test(${test})
    set_property(TEST ${test} PROPERTY LABELS "mc;alea;hdf5;tutorial")
  endforeach()

  if(ALPS_ENABLE_MPI)
    add_executable(mcbase_mpi tests/mpi.cpp)
    target_link_libraries(mcbase_mpi PRIVATE ALPS::alps MPI::MPI_CXX)
    add_test(NAME mcbase_mpi COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:mcbase_mpi> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(mcbase_mpi PROPERTIES LABELS "mc;alea;mpi" PROCESSORS 2 TIMEOUT 60)
  endif()
endif()
