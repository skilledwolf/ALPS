# Existing parameter grammar and HDF5 parameter preparation tools.
add_executable(parameter2hdf5 parameter2hdf5.C)
target_link_libraries(parameter2hdf5 PRIVATE ALPS::alps)
add_executable(p2h5 p2h5.cpp)
target_link_libraries(p2h5 PRIVATE ALPS::alps)
install(TARGETS parameter2hdf5 p2h5
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

alps_register_module(NAME tools_parameters
  TARGETS parameter2hdf5 p2h5
  DEPENDS hdf5 legacy_parameters xml)
