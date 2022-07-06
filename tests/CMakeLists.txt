cmake_minimum_required(VERSION 3.14)
cmake_policy(VERSION 3.14...3.21)

project(tests)
set(CMAKE_CXX_STANDARD 17)
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        main
)
# Prevent overriding the parent project's compiler/linker settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

set(USE_GTEST "no" CACHE STRING "Should we provide the tests using gtest framework")
if(USE_GTEST)
  FetchContent_MakeAvailable(googletest)
  enable_testing()
  add_executable(cache cache.cxx)
  set_target_properties(cache PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
  target_link_libraries(cache mdl2 gtest_main)
  add_executable(swaplocal swaplocal.cxx)
  set_target_properties(swaplocal PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
  target_link_libraries(swaplocal mdl2 gtest_main)
  add_test(NAME cache COMMAND $<TARGET_FILE:cache> WORKING_DIRECTORY ${CMAKE_BINARY_DIR}) 
  add_test(NAME mpicache COMMAND mpirun -n 2 $<TARGET_FILE:cache> WORKING_DIRECTORY ${CMAKE_BINARY_DIR}) 
  add_test(NAME swaplocal COMMAND $<TARGET_FILE:swaplocal> WORKING_DIRECTORY ${CMAKE_BINARY_DIR}) 
endif()
