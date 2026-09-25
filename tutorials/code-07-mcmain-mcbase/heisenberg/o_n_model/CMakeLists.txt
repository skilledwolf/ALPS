cmake_minimum_required(VERSION 3.27...4.3)
project(ndim_spin LANGUAGES C CXX)

find_package(ALPS CONFIG REQUIRED)

# rule for generating the heisenberg example program
add_executable(heisenberg heisenberg.cpp)
target_link_libraries(heisenberg PRIVATE ALPS::alps)

option(ALPS_BUILD_TESTING "Build and register ALPS tests" ON)
if(ALPS_BUILD_TESTING)
  enable_testing()
  add_executable(tinyvector_test tinyvector_test.cpp)
  target_link_libraries(tinyvector_test PRIVATE ALPS::alps)
  add_test(NAME tinyvector COMMAND tinyvector_test)
  set_tests_properties(tinyvector PROPERTIES ENVIRONMENT_MODIFICATION
    "PATH=path_list_prepend:$<TARGET_FILE_DIR:ALPS::alps>")
endif()
