cmake_minimum_required(VERSION 3.27...4.3)
project(alpsize LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)
include("${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/SpinParameters.cmake")

# rule for generating single core ising program
add_executable(single ising.cpp single.cpp)
target_link_libraries(single PRIVATE ALPS::alps ALPS::run_config)

# alps::mc::main spreads the chains over MPI ranks when ALPS is built with MPI
# and runs them in one process otherwise.
add_executable(mpi ising.cpp mpi.cpp)
target_link_libraries(mpi PRIVATE ALPS::alps)
