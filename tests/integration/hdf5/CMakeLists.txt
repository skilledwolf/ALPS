# The legacy parameter bridge still has model-library consumers.
# Native ALEA archive contracts live in src/alps/alea/tests/statistics.cpp.
add_executable(hdf5_parms hdf5_parms.cpp)
target_link_libraries(hdf5_parms PRIVATE ALPS::alps)
alps_add_test(hdf5_parms)
set_property(TEST hdf5_parms PROPERTY LABELS "hdf5;integration")
