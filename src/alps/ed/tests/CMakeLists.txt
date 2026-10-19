# Unit tests of the ed component. The numeric locks include the models test
# helper dense_e0.hpp (an independent Kronecker assembly over LAPACK's zheev).
set(models_tests "${PROJECT_SOURCE_DIR}/src/alps/models/tests")

foreach(name IN ITEMS
    canonical_csr
    local_transitions
    product_basis
    row_csr
    streamed_local_operator
    term_accumulation
    thermodynamics)
  add_executable(ed_${name} ${name}.cpp)
  target_link_libraries(ed_${name} PRIVATE ALPS::ed)
  alps_add_gtest(ed_${name} COMPONENT ed LABELS "ed;unit")
endforeach()

foreach(name IN ITEMS
    conserved_sectors
    dense_eigensystem
    finite_group
    ground_states
    lanczos_csr
    momentum_compiled
    product_problem
    projected_basis
    projected_problem
    sparse_from_local_lock)
  add_executable(ed_${name} ${name}.cpp)
  target_link_libraries(ed_${name} PRIVATE ALPS::ed ${ALPS_LAPACK_TARGET} ${ALPS_BLAS_TARGET})
  target_include_directories(ed_${name} PRIVATE "${models_tests}")
  alps_add_gtest(ed_${name} COMPONENT ed LABELS "ed;numeric" TIMEOUT 300)
endforeach()
