cmake_minimum_required(VERSION 3.5)
cmake_policy(VERSION 3.5...3.14)

project(tests)
find_package(GTest)

if (GTEST_FOUND)
  enable_testing()

  add_executable(cache cache.cxx)
  set_target_properties(cache PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
  target_link_libraries(cache mdl2 GTest::GTest)
  add_test(NAME cache COMMAND $<TARGET_FILE:cache> WORKING_DIRECTORY ${CMAKE_BINARY_DIR}) 
  add_test(NAME mpicache COMMAND mpirun -n 2 $<TARGET_FILE:cache> WORKING_DIRECTORY ${CMAKE_BINARY_DIR}) 

endif(GTEST_FOUND)
