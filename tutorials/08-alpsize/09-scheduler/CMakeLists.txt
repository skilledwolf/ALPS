cmake_minimum_required(VERSION 3.27...4.3)
project(alpsize LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating 'hello world' program
add_executable(hello hello.C)
target_link_libraries(hello PRIVATE ALPS::alps)

# rule for generating wolff program
add_executable(wolff wolff.C)
target_link_libraries(wolff PRIVATE ALPS::alps)

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/execution.toml" WOLFF_LATTICE_SCHEMA)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/execution.toml")
include("${CMAKE_CURRENT_SOURCE_DIR}/../Wolff.cmake")
