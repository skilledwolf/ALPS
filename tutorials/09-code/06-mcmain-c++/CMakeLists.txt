cmake_minimum_required(VERSION 3.27...4.3)
project(alps_ising_tutorial LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)
include("${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/SpinParameters.cmake")
set(ISING_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../10-ngs/1_accumulator_only")
add_executable(ising "${ISING_SOURCE}/ising.cpp" main.cpp)
target_include_directories(ising PRIVATE "${ISING_SOURCE}")
target_link_libraries(ising PRIVATE ALPS::alps ALPS::run_config)
