cmake_minimum_required(VERSION 3.27...4.3)
project(alpsize LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating wolff program
add_executable(wolff wolff.C)

include("${CMAKE_CURRENT_SOURCE_DIR}/../Wolff.cmake")
