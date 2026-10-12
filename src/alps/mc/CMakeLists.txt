# Mc source ownership.
target_sources(alps PRIVATE
  src/mcbase.cpp
  src/ngs/api.cpp
  src/stop_callback.cpp)

target_sources(alps_headers INTERFACE
  FILE_SET mc_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/check_schedule.hpp
    include/alps/mcbase.hpp
    include/alps/mcmpiadapter.hpp
    include/alps/ngs/detail/mpi_agreement.hpp
    include/alps/ngs.hpp
    include/alps/ngs/api.hpp
    include/alps/progress_callback.hpp
    include/alps/stop_callback.hpp)

# Public Python integration spellings survive the move to Python ownership.
target_sources(alps_headers INTERFACE FILE_SET python_compatibility_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/ngs/detail/export_sim_to_python.hpp
        include/alps/python/save_observable_to_hdf5.hpp)

if(ALPS_BUILD_TESTING AND ALPS_ENABLE_MPI)
  add_executable(mc_legacy_collect tests/legacy_collect.cpp)
  target_link_libraries(mc_legacy_collect PRIVATE ALPS::alps)
  alps_add_mpi_gtest(mc_legacy_collect COMPONENT mc RANKS 2 LABELS regression)
  alps_add_mpi_gtest(mc_legacy_collect COMPONENT mc RANKS 3 LABELS regression)
endif()
