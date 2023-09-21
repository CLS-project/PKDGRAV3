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

  add_executable(imf imf.cxx)
  target_include_directories(imf PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../ ${CMAKE_CURRENT_SOURCE_DIR}/../)
  set_target_properties(imf PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
  target_link_libraries(imf mdl2 gtest_main)
  target_link_libraries(imf blitz fmt)
  add_test(NAME imf COMMAND $<TARGET_FILE:imf> WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable(hydro hydro.cxx)
   # check_symbol_exists does not work for intrinsics..
   include(CheckCSourceCompiles)
   set(CMAKE_REQUIRED_FLAGS -march=native)
   check_c_source_compiles("
       #include <immintrin.h>
       int main() {
         __m128d a, b;
         __m128d tmp = _mm_pow_pd(a, b);
         return 0;
       }
     "
     HAVE_MM_POW
   )
   if(HAVE_MM_POW)
      target_compile_options(hydro PRIVATE -DHAVE_MM_POW)
   endif()
   check_c_source_compiles("
       #include <immintrin.h>
       int main() {
         __m256d a, b;
         __m256d tmp = _mm256_pow_pd(a, b);
         return 0;
       }
     "
     HAVE_MM256_POW
   )
   if(HAVE_MM256_POW)
      target_compile_options(hydro PRIVATE -DHAVE_MM256_POW)
   endif()
   check_c_source_compiles("
       #include <immintrin.h>
       int main() {
         __m512d a, b;
         __m512d tmp = _mm512_pow_pd(a, b);
         return 0;
       }
     "
     HAVE_MM512_POW
   )
   if(HAVE_MM512_POW)
      target_compile_options(hydro PRIVATE -DHAVE_MM512_POW)
   endif()
  target_include_directories(hydro PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../ ${CMAKE_CURRENT_SOURCE_DIR}/../)
  set_target_properties(hydro PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
  set_source_files_properties(hydro.cxx PROPERTIES COMPILE_OPTIONS "-fno-finite-math-only")
  target_link_libraries(hydro mdl2 gtest_main)
  target_link_libraries(hydro blitz fmt)
  add_test(NAME hydro COMMAND $<TARGET_FILE:hydro> WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable(limiter limiters.cxx)
  target_include_directories(limiter PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../ ${CMAKE_CURRENT_SOURCE_DIR}/../)
  set_target_properties(limiter PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
  target_link_libraries(limiter mdl2 gtest_main)
  target_link_libraries(limiter blitz fmt)
  add_test(NAME limiter COMMAND $<TARGET_FILE:limiter> WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  add_executable(eEOS eEOS.cxx)
  target_include_directories(eEOS PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/../ ${CMAKE_CURRENT_SOURCE_DIR}/../)
  set_target_properties(eEOS PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES CXX_EXTENSIONS NO)
  target_link_libraries(eEOS mdl2 gtest_main)
  target_link_libraries(eEOS blitz fmt)
  add_test(NAME eEOS COMMAND $<TARGET_FILE:eEOS> WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
