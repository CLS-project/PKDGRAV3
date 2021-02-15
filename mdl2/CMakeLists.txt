cmake_minimum_required(VERSION 3.5)
cmake_policy(VERSION 3.5...3.14)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Release build with debug info selected")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
project(mdl2 LANGUAGES C CXX Fortran)
add_subdirectory(openpa)	# Lockless queue support
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
include(FindPkgConfig)
find_package(FFTW)
if(FFTW_FOUND)
  set(MDL_FFTW ${FFTW_FOUND})
endif()
pkg_check_modules(HWLOC hwloc)
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILES (time.h HAVE_TIME_H)
CHECK_INCLUDE_FILES (sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES (unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES (signal.h HAVE_SIGNAL_H)
CHECK_INCLUDE_FILES (inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES (stdint.h HAVE_STDINT_H)

INCLUDE(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(memkind hbw_posix_memalign "" HAVE_MEMKIND)

INCLUDE(CheckSymbolExists)
check_symbol_exists(backtrace execinfo.h USE_BT)

find_package(Threads REQUIRED)
find_package(MPI REQUIRED)	# MPI support
find_package(CUDA)
find_package(elfutils)

add_library(${PROJECT_NAME} STATIC "")
target_compile_features(${PROJECT_NAME} PRIVATE cxx_auto_type cxx_range_for cxx_lambdas cxx_strong_enums)
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

target_sources(${PROJECT_NAME}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mdl.cxx
          ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mdlmessages.cxx
          ${CMAKE_CURRENT_SOURCE_DIR}/mdlbase.cxx
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mdl.h
          ${CMAKE_CURRENT_SOURCE_DIR}/mdlbase.h
          ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h
)
if(USE_BT)
  target_sources(${PROJECT_NAME}
	  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/mdlbt.cxx
	  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mdlbt.h)
  set(USE_BT TRUE)
endif()

if (HWLOC_FOUND)
  find_library(HWLOC_LIBRARY hwloc HINTS ${HWLOC_LIBDIR})
  CHECK_LIBRARY_EXISTS(hwloc hwloc_topology_init ${HWLOC_LIBDIR} HAVE_HWLOC)
  if (HAVE_HWLOC)
    set(USE_HWLOC TRUE)
    target_include_directories(${PROJECT_NAME} PUBLIC ${HWLOC_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${HWLOC_LIBRARY})
  else()
    message(WARNING
	" Found hwloc library but cannot link to it, so we won't use it.\n"
	" On Cray you can try: export CRAYPE_LINK_TYPE=dynamic")
  endif()
endif()
if (elfutils_FOUND)
   set(USE_ELFUTILS TRUE)
   target_link_libraries(${PROJECT_NAME} elfutils::dw elfutils::elf)
endif()

if (APPLE)
  target_sources(${PROJECT_NAME}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/mac/pthread_barrier.c
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mac/pthread_barrier.h
  )
endif()
if (MSVC)
  target_sources(${PROJECT_NAME}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/windows/pthread.c
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/windows/pthread.h
  )
  target_link_libraries(${PROJECT_NAME} wsock32 ws2_32)
endif()
if (APPLE)
  target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mac)
endif()
if (MSVC)
  target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/windows)
endif()
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/openpa)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/boost)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mpi)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_C_INCLUDE_PATH})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..) ############################
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/..)

target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} openpa)

if(FFTW_FOUND)
target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC ${FFTW_INCLUDES})
endif()
if(HAVE_MEMKIND)
find_library(MEMKIND_LIBRARY memkind)
target_link_libraries(${PROJECT_NAME} ${MEMKIND_LIBRARY})
endif()

if (CUDA_FOUND)
  set(USE_CUDA TRUE)
  get_property(MDL_INCLUDES TARGET mdl2 PROPERTY INCLUDE_DIRECTORIES)
  get_property(OPENPA_INCLUDES TARGET openpa PROPERTY INCLUDE_DIRECTORIES)
  CUDA_INCLUDE_DIRECTORIES(${MDL_INCLUDES} ${OPENPA_INCLUDES})
#  CUDA_INCLUDE_DIRECTORIES(mdl2 openpa)
  target_link_libraries(${PROJECT_NAME} ${CUDA_LIBRARIES})
  CUDA_COMPILE(cuda_files ${CMAKE_CURRENT_SOURCE_DIR}/cuda/mdlcuda.cu OPTIONS -arch compute_35)
  target_sources(${PROJECT_NAME} PRIVATE ${cuda_files})
  target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/cuda ${CUDA_INCLUDE_DIRS})
endif(CUDA_FOUND)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/mdl_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h)
install(TARGETS ${PROJECT_NAME} DESTINATION "lib")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h mpi/mdl.h mdlbase.h DESTINATION "include")
