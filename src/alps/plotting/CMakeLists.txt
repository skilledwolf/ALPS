# Plot output writes typed parameters and data sets as ALPS XML.
target_sources(alps_headers INTERFACE
  FILE_SET plotting_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/plot.h)

alps_register_module(NAME plotting HEADER_SET plotting_headers
  DEPENDS params xml
)
