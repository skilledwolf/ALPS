# DMFT's in-process scheduler solvers and lattice tutorials still use the legacy Parameters API.
target_sources(alps PRIVATE make_deprecated_parameters.cpp)
target_sources(alps_headers INTERFACE
  FILE_SET params_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES include/alps/ngs/make_deprecated_parameters.hpp)
