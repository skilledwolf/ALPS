# Typed parameters reach the expression, lattice and model libraries as text symbols.
add_executable(param_symbol_table symbol_table.cpp)
target_link_libraries(param_symbol_table PRIVATE ALPS::alps)
alps_add_test(param_symbol_table)
set_property(TEST param_symbol_table PROPERTY LABELS "params;integration")
