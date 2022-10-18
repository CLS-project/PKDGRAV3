cmake_minimum_required(VERSION 3.8)
cmake_policy(VERSION 3.8...3.22)
include(CheckLanguage)
include(FindPkgConfig)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Release build with debug info selected")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()
project(mdl2 LANGUAGES C CXX)
option(DEBUG_COUNT_CACHE "Count messages sent and received in the cache" OFF)
get_directory_property(hasParent PARENT_DIRECTORY)
add_subdirectory(cuda)
if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/metal")
  add_subdirectory(metal)
endif()
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
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
INCLUDE (CheckFunctionExists)
CHECK_FUNCTION_EXISTS(sched_yield HAVE_SCHED_YIELD)

INCLUDE(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(memkind hbw_posix_memalign "" HAVE_MEMKIND)

INCLUDE(CheckSymbolExists)
check_symbol_exists(backtrace execinfo.h USE_BT)

find_package(Threads REQUIRED)
find_package(MPI REQUIRED)	# MPI support
find_package(elfutils QUIET)
find_package(Boost REQUIRED)

add_library(${PROJECT_NAME} STATIC "")
target_compile_features(${PROJECT_NAME} PRIVATE cxx_auto_type cxx_range_for cxx_lambdas cxx_strong_enums)
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

target_link_libraries(${PROJECT_NAME} mdlcuda)
target_sources(${PROJECT_NAME}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mdl.cxx
          ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mpimessages.cxx
          ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mdlmessages.cxx
          ${CMAKE_CURRENT_SOURCE_DIR}/gpu/mdlgpu.cxx
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
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mpi)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..) ############################
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/..)
if(CMAKE_VERSION VERSION_LESS "3.10")
target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_C_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES})
else()
target_link_libraries(${PROJECT_NAME} MPI::MPI_C)
target_include_directories(${PROJECT_NAME} PUBLIC ${MPI_C_INCLUDE_DIRS})
endif()
target_link_libraries(${PROJECT_NAME} ${CMAKE_THREAD_LIBS_INIT})

if(FFTW_FOUND)
target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC ${FFTW_INCLUDES})
endif()
if(HAVE_MEMKIND)
find_library(MEMKIND_LIBRARY memkind)
target_link_libraries(${PROJECT_NAME} ${MEMKIND_LIBRARY})
endif()
target_link_libraries(${PROJECT_NAME} Boost::boost)

if (APPLE AND METAL AND FOUNDATION AND QUARTZCORE)
  target_link_libraries(${PROJECT_NAME} mdlmetal)
endif (APPLE AND METAL AND FOUNDATION AND QUARTZCORE)
if(hasParent)
  set(USE_CUDA ${USE_CUDA} PARENT_SCOPE)
  set(USE_METAL ${USE_METAL} PARENT_SCOPE)
endif(hasParent)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/mdl_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h)
install(TARGETS ${PROJECT_NAME} DESTINATION "lib")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h mpi/mdl.h mdlbase.h DESTINATION "include")
