cmake_minimum_required(VERSION 3.27...4.3)
project(alps_ising_tutorial LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)
include("${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/SpinParameters.cmake")
add_executable(ising ising.cpp main.cpp)
target_link_libraries(ising PRIVATE ALPS::alps ALPS::run_config)
