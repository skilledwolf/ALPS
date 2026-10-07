# Scheduler source ownership.
target_sources(alps PRIVATE
  src/scheduler/abstract_task.C
  src/scheduler/convertxdr.C
  src/scheduler/factory.C
  src/scheduler/info.C
  src/scheduler/master_scheduler.C
  src/scheduler/mpp_scheduler.C
  src/scheduler/options.C
  src/scheduler/remote_task.C
  src/scheduler/remote_worker.C
  src/scheduler/scheduler.C
  src/scheduler/serial_scheduler.C
  src/scheduler/signal.C
  src/scheduler/single_scheduler.C
  src/scheduler/slave_task.C
  src/scheduler/task.C
  src/scheduler/worker.C
  src/scheduler/workertask.C)

target_sources(alps_headers INTERFACE
  FILE_SET scheduler_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/scheduler.h
    include/alps/scheduler/convert.h
    include/alps/scheduler/diag.hpp
    include/alps/scheduler/factory.h
    include/alps/scheduler/info.h
    include/alps/scheduler/options.h
    include/alps/scheduler/scheduler.h
    include/alps/scheduler/signal.hpp
    include/alps/scheduler/task.h
    include/alps/scheduler/types.h
    include/alps/scheduler/worker.h)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME scheduler HEADER_SET scheduler_headers
  TARGETS alps
  DEPENDS numeric_io configuration expression hdf5 lattice legacy_parameters model numerics osiris params random utilities xml
)
