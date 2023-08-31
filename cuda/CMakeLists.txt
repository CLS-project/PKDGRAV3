cmake_minimum_required(VERSION 3.10)
cmake_policy(VERSION 3.10...3.22)

include(CheckLanguage)

project(pkdcuda)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
  enable_language(CUDA)
  set(USE_CUDA TRUE)
  set(USE_CUDA ${USE_CUDA} PARENT_SCOPE)
  add_library(${PROJECT_NAME} STATIC cudaewald.cu cudapppc.cu cudautil.cu)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
  )
  target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(${PROJECT_NAME} mdl2 blitz)
else(CMAKE_CUDA_COMPILER)
  add_library(${PROJECT_NAME} INTERFACE)
endif(CMAKE_CUDA_COMPILER)

