# Native SDK migration

The native SDK now exports components with CMake package targets. CMake 3.27, C++17, Boost 1.76, HDF5 1.10.5 and LP64 BLAS/LAPACK are required. MPI defaults to OFF; existing CMake caches preserve their configured value.
