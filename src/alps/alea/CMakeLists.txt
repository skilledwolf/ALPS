# Alea source ownership.
target_sources(alps PRIVATE
  src/alea/nan.C
  src/alea/observable.C
  src/alea/observablefactory.C
  src/alea/observableset.C
  src/ngs/mcobservable.cpp
  src/ngs/mcobservables.cpp
  src/ngs/mcresult.cpp
  src/ngs/mcresults.cpp
  src/ngs/observablewrappers.cpp)

target_sources(alps_headers INTERFACE
  FILE_SET alea_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/alea.h
    include/alps/alea/abstractbinning.h
    include/alps/alea/abstractsimpleobservable.h
    include/alps/alea/abstractsimpleobservable.ipp
    include/alps/alea/convergence.hpp
    include/alps/alea/detailedbinning.h
    include/alps/alea/histogram.h
    include/alps/alea/histogramdata.h
    include/alps/alea/histogrameval.h
    include/alps/alea/mcanalyze.hpp
    include/alps/alea/mcdata.hpp
    include/alps/alea/nan.h
    include/alps/alea/nobinning.h
    include/alps/alea/observable.h
    include/alps/alea/observable_fwd.hpp
    include/alps/alea/observablefactory.h
    include/alps/alea/observableset.h
    include/alps/alea/observableset_p.h
    include/alps/alea/output_helper.h
    include/alps/alea/recordableobservable.h
    include/alps/alea/signedobservable.h
    include/alps/alea/simplebinning.h
    include/alps/alea/simpleobsdata.h
    include/alps/alea/simpleobservable.h
    include/alps/alea/simpleobservable.ipp
    include/alps/alea/simpleobseval.h
    include/alps/alea/simpleobseval.ipp
    include/alps/alea/type_tag.hpp
    include/alps/alea/value_with_error.hpp
    include/alps/ngs/lib/mcresult_impl_base.ipp
    include/alps/ngs/lib/mcresult_impl_derived.ipp
    include/alps/ngs/mcobservable.hpp
    include/alps/ngs/mcobservables.hpp
    include/alps/ngs/mcresult.hpp
    include/alps/ngs/mcresults.hpp
    include/alps/ngs/observablewrappers.hpp)

# Modern Eigen-based ALEA from pinned ALPSCore; a separate scientific runtime.
add_library(alps_statistics
  src/statistics/autocorr.cpp
  src/statistics/batch.cpp
  src/statistics/covariance.cpp
  src/statistics/fisher_f.cpp
  src/statistics/galois.cpp
  src/statistics/mean.cpp
  src/statistics/model.cpp
  src/statistics/propagation.cpp
  src/statistics/result.cpp
  src/statistics/testing.cpp
  src/statistics/transform.cpp
  src/statistics/util.cpp
  src/statistics/variance.cpp
)
add_library(ALPS::statistics ALIAS alps_statistics)
target_compile_features(alps_statistics PUBLIC cxx_std_17)
set_target_properties(alps_statistics PROPERTIES
  EXPORT_NAME statistics
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  WINDOWS_EXPORT_ALL_SYMBOLS ON
  INSTALL_RPATH "${ALPS_LOADER_ORIGIN}")
target_link_libraries(alps_statistics PUBLIC Eigen3::Eigen Boost::headers)
alps_copy_runtime_dependencies(alps_statistics)
target_sources(alps_statistics PUBLIC
  FILE_SET alea_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/alea.hpp
    include/alps/alea/autocorr.hpp
    include/alps/alea/batch.hpp
    include/alps/alea/bundle.hpp
    include/alps/alea/checkpoint.hpp
    include/alps/alea/complex_op.hpp
    include/alps/alea/computed.hpp
    include/alps/alea/convert.hpp
    include/alps/alea/core.hpp
    include/alps/alea/covariance.hpp
    include/alps/alea/fisher_f.hpp
    include/alps/alea/hdf5.hpp
    include/alps/alea/internal/format.hpp
    include/alps/alea/internal/galois.hpp
    include/alps/alea/internal/joined.hpp
    include/alps/alea/internal/outer.hpp
    include/alps/alea/internal/pooling.hpp
    include/alps/alea/internal/util.hpp
    include/alps/alea/mean.hpp
    include/alps/alea/propagation.hpp
    include/alps/alea/result.hpp
    include/alps/alea/serialize.hpp
    include/alps/alea/testing.hpp
    include/alps/alea/transform.hpp
    include/alps/alea/transformer.hpp
    include/alps/alea/util.hpp
    include/alps/alea/util/model.hpp
    include/alps/alea/util/serializer.hpp
    include/alps/alea/var_strategy.hpp
    include/alps/alea/variance.hpp
    include/alps/common/view.hpp
    include/alps/serialization/core.hpp
    include/alps/serialization/eigen.hpp
    include/alps/serialization/util.hpp
)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME alea HEADER_SET alea_headers
  HEADER_TARGETS alps_headers alps_statistics
  TARGETS alps alps_statistics
  DEPENDS configuration containers hdf5 numerics osiris utilities xml
)
