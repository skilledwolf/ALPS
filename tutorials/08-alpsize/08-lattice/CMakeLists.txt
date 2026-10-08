cmake_minimum_required(VERSION 3.27...4.3)
project(alpsize LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating lattice program
add_executable(lattice lattice.C)
target_link_libraries(lattice PRIVATE ALPS::alps ALPS::run_config)

# rule for generating wolff program
add_executable(wolff wolff.C)

set(WOLFF_LATTICE_SCHEMA "
[parameters.LATTICE]
type = \"string\"
[parameters.GRAPH]
type = \"string\"
[input.lattice_library]
type = \"path\"
")
target_link_libraries(wolff PRIVATE ALPS::alps)

include("${CMAKE_CURRENT_SOURCE_DIR}/../Wolff.cmake")
