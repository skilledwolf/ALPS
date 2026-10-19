foreach(name IN ITEMS
    ab_initio
    compile_hamiltonian
    compile_on_physical_basis
    finite_graph
    hermiticity
    heterogeneous_physical_basis
    lattice_tiling
    lattice_xml
    operator_decomposition
    operator_names
    qubit_encoding
    sign
    term_sum
    translation_generators)
  add_executable(models_${name} ${name}.cpp)
  target_link_libraries(models_${name} PRIVATE ALPS::models)
  alps_add_gtest(models_${name} COMPONENT models LABELS "models;unit")
endforeach()

# The Pauli locks assemble dense Hamiltonians and call LAPACK directly.
foreach(name IN ITEMS pauli_export pauli_fermions)
  add_executable(models_${name} ${name}.cpp)
  target_link_libraries(models_${name} PRIVATE ALPS::models ${ALPS_LAPACK_TARGET} ${ALPS_BLAS_TARGET})
  alps_add_gtest(models_${name} COMPONENT models LABELS "models;numeric" TIMEOUT 120)
endforeach()
