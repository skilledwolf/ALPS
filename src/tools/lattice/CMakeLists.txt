foreach(program IN ITEMS lattice2xml printgraph)
  add_executable(${program} ${program}.C)
  target_link_libraries(${program} PRIVATE ALPS::alps ALPS::run_config)
endforeach()
install(TARGETS lattice2xml printgraph
  RUNTIME_DEPENDENCY_SET alps_runtime
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT tools)

# pltgraph.py remains source-only.
alps_register_module(NAME tools_lattice
  TARGETS lattice2xml printgraph
  PRIVATE_FILES lattice_run.hpp
  DEPENDS lattice legacy_parameters run_config xml)
