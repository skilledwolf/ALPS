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
    include/alps/ngs.hpp
    include/alps/ngs/api.hpp
    include/alps/progress_callback.hpp
    include/alps/stop_callback.hpp)

alps_register_module(NAME mc HEADER_SET mc_headers
  TARGETS alps
  DEPENDS accumulators alea cli hdf5 params random utilities
)

if(ALPS_BUILD_TESTING)
  add_executable(mcbase_checkpoint tests/checkpoint.cpp)
  target_link_libraries(mcbase_checkpoint PRIVATE ALPS::alps)
  alps_add_test(mcbase_checkpoint)
  set_property(TEST mcbase_checkpoint PROPERTY LABELS "mc;hdf5")
endif()
