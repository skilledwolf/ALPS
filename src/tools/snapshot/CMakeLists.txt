# Offline converter for snapshot files written by released ALPS versions.
add_executable(snap2vtk snap2vtk.C)
target_link_libraries(snap2vtk PRIVATE ALPS::alps)
install(TARGETS snap2vtk
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

alps_register_module(NAME tools_snapshot
  TARGETS snap2vtk
  DEPENDS osiris utilities)
