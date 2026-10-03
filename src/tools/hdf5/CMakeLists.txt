install(PROGRAMS convert.py
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  RENAME alps-hdf5-convert
  COMPONENT tools)

alps_register_module(NAME tools_hdf5
  EXTRA_FILES convert.py)
