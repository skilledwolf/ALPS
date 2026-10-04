#  Copyright Matthias Troyer 2009 - 2013.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

set(DMFT_QMC_SOURCE
    run_config.C
    externalsolver.C
    fouriertransform.C
    hilberttransformer.C
    hirschfyesim.C
    hirschfye_run.C
    alps_solver.C
    bandstructure.C
    auxiliaryfunctions.C)
add_library(dmft_qmc_impl STATIC ${DMFT_QMC_SOURCE})

set(DMFT_QMC_INTERACTION_EXPANSION_SOURCE
    interaction_expansion/auxiliary.cpp
    interaction_expansion/observables.cpp
    interaction_expansion/fastupdate.cpp
    interaction_expansion/selfenergy.cpp
    interaction_expansion/solver.cpp
    interaction_expansion/io.cpp
    interaction_expansion/splines.cpp
    interaction_expansion/interaction_expansion.cpp
    interaction_expansion/measurements.cpp
    interaction_expansion/model.cpp)

add_library(alps_cthyb STATIC
    hybridization/run.cpp
    hybridization/hybsim.cpp
    hybridization/hyblocal.cpp
    hybridization/hybint.cpp
    hybridization/hybfun.cpp
    hybridization/hybretintfun.cpp
    hybridization/hybmatrix.cpp
    hybridization/hybmatrix_ft.cpp
    hybridization/hybconfig.cpp
    hybridization/hybupdates.cpp
    hybridization/hybevaluate.cpp
    hybridization/hybmeasurements.cpp)
add_library(alps_ctint STATIC
    interaction_expansion2/run.cpp
    fouriertransform.C
    interaction_expansion2/auxiliary.cpp
    interaction_expansion2/observables.cpp
    interaction_expansion2/fastupdate.cpp
    interaction_expansion2/selfenergy.cpp
    interaction_expansion2/solver.cpp
    interaction_expansion2/io.cpp
    interaction_expansion2/splines.cpp
    interaction_expansion2/interaction_expansion.cpp
    interaction_expansion2/measurements.cpp
    interaction_expansion2/model.cpp)
foreach(solver IN ITEMS cthyb ctint)
  add_library(ALPS::${solver} ALIAS alps_${solver})
  set_target_properties(alps_${solver} PROPERTIES EXPORT_NAME ${solver}
    CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
  target_link_libraries(alps_${solver} PUBLIC ALPS::alps)
  target_include_directories(alps_${solver} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
endforeach()
add_subdirectory(hybridization)
add_subdirectory(interaction_expansion2)
install(TARGETS alps_cthyb alps_ctint EXPORT ALPSTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libraries)

#Interaction Expansion
add_library(interaction_expansion_impl STATIC ${DMFT_QMC_INTERACTION_EXPANSION_SOURCE})
target_link_libraries(interaction_expansion_impl PUBLIC ALPS::alps ALPS::run_config)

find_package(tomlplusplus 3.4 CONFIG REQUIRED)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema/dmft.toml" ALPS_DMFT_SCHEMA)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/schema/hirschfye.toml" ALPS_HIRSCHFYE_SCHEMA)
configure_file(schema.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/dmft_schema.hpp" @ONLY)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/schema/dmft.toml" "${CMAKE_CURRENT_SOURCE_DIR}/schema/hirschfye.toml")
target_include_directories(dmft_qmc_impl
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}"
  PRIVATE "$<TARGET_PROPERTY:tomlplusplus::tomlplusplus,INTERFACE_INCLUDE_DIRECTORIES>")
target_link_libraries(dmft_qmc_impl PUBLIC ALPS::alps ALPS::cthyb ALPS::ctint ALPS::run_config ALPS::statistics)
if(ALPS_ENABLE_MPI)
  target_link_libraries(dmft_qmc_impl PRIVATE MPI::MPI_CXX)
endif()
install(FILES schema/dmft.toml schema/hirschfye.toml
  DESTINATION ${CMAKE_INSTALL_DATADIR}/alps/schemas COMPONENT applications)

target_include_directories(interaction_expansion_impl PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

add_executable(dmft main.C selfconsistency.C)
add_executable(hirschfye solver_main.C)

add_executable(hybridization hybridization/hybmain.cpp)
add_executable(interaction interaction_expansion2/main.cpp)

target_link_libraries(dmft PRIVATE dmft_qmc_impl interaction_expansion_impl)
target_link_libraries(hirschfye PRIVATE dmft_qmc_impl)

target_link_libraries(hybridization PRIVATE ALPS::cthyb)
target_link_libraries(interaction PRIVATE ALPS::ctint)

if(ALPS_BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_test(NAME dmft_cli_contract COMMAND ${Python3_EXECUTABLE}
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/dmft_cli_contract.py"
    $<TARGET_FILE:dmft> $<TARGET_FILE:hybridization> $<TARGET_FILE:dmft_externalsolver_protocol>)
  set_tests_properties(dmft_cli_contract PROPERTIES TIMEOUT 180 LABELS "dmft;run_config")
  add_executable(dmft_hirschfye_fixture tests/hirschfye_fixture.cpp)
  target_link_libraries(dmft_hirschfye_fixture PRIVATE dmft_qmc_impl)
  add_test(NAME hirschfye_cli_contract COMMAND ${Python3_EXECUTABLE}
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/hirschfye_cli_contract.py"
    $<TARGET_FILE:hirschfye> $<TARGET_FILE:dmft_hirschfye_fixture>)
  set_tests_properties(hirschfye_cli_contract PROPERTIES TIMEOUT 60 LABELS "dmft;run_config")
  add_executable(hirschfye_statistics tests/hirschfye_statistics.cpp)
  target_link_libraries(hirschfye_statistics PRIVATE dmft_qmc_impl)
  alps_add_test(hirschfye_statistics)
  set_property(TEST hirschfye_statistics PROPERTY LABELS "dmft;hirschfye")
  if(ALPS_ENABLE_MPI)
    add_executable(hirschfye_mpi_contract tests/hirschfye_mpi.cpp)
    target_link_libraries(hirschfye_mpi_contract PRIVATE dmft_qmc_impl MPI::MPI_CXX)
    add_test(NAME hirschfye_mpi_contract COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:hirschfye_mpi_contract> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(hirschfye_mpi_contract PROPERTIES LABELS "dmft;hirschfye;mpi" PROCESSORS 2 TIMEOUT 60)
  endif()
  # Regression tests for the shared DMFT QMC infrastructure.
  add_executable(dmft_fouriertransform_offdiag_numeric dmft_fouriertransform_offdiag_numeric.C)
  target_link_libraries(dmft_fouriertransform_offdiag_numeric PRIVATE dmft_qmc_impl)
  alps_add_test(dmft_fouriertransform_offdiag_numeric)
  set_property(TEST dmft_fouriertransform_offdiag_numeric PROPERTY LABELS dmft)

  add_executable(dmft_externalsolver_protocol dmft_externalsolver_protocol.C)
  target_link_libraries(dmft_externalsolver_protocol PRIVATE dmft_qmc_impl)
  alps_add_test(dmft_externalsolver_protocol)
  set_property(TEST dmft_externalsolver_protocol PROPERTY LABELS dmft)

  add_executable(dmft_hilberttransformer_afm_multiband_numeric
                 dmft_hilberttransformer_afm_multiband_numeric.C)
  target_link_libraries(dmft_hilberttransformer_afm_multiband_numeric PRIVATE dmft_qmc_impl)
  alps_add_test(dmft_hilberttransformer_afm_multiband_numeric)
  set_property(TEST dmft_hilberttransformer_afm_multiband_numeric PROPERTY LABELS dmft)

  add_executable(dmft_interaction_expansion_choice dmft_interaction_expansion_choice.C)
  alps_add_test(dmft_interaction_expansion_choice)
  set_property(TEST dmft_interaction_expansion_choice PROPERTY LABELS dmft)
endif()

alps_install_applications(dmft hirschfye hybridization interaction)
