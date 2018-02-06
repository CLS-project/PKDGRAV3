cmake_minimum_required(VERSION 3.1)
project(mdl2)
add_subdirectory(openpa)	# Lockless queue support

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/")
find_package(FFTW REQUIRED)
if(FFTW_FOUND)
#  INCLUDE_DIRECTORIES(${FFTW_INCLUDES})
  set(MDL_FFTW ${FFTW_FOUND})
endif()

INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILES (sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES (unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES (signal.h HAVE_SIGNAL_H)
CHECK_INCLUDE_FILES (inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES (stdint.h HAVE_STDINT_H)

INCLUDE(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(memkind hbw_posix_memalign "" HAVE_MEMKIND)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/mdl_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h)

add_library(${PROJECT_NAME} STATIC "")
target_sources(${PROJECT_NAME}
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mdl.c
          ${CMAKE_CURRENT_SOURCE_DIR}/mdlbase.c
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mpi/mdl.h
          ${CMAKE_CURRENT_SOURCE_DIR}/mdlbase.h
          ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h
)
if (APPLE)
  target_sources(${PROJECT_NAME}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/mac/pthread_barrier.c
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/mac/pthread_barrier.h
  )
endif()

if (APPLE)
  target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mac)
endif()
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mpi)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_C_INCLUDE_PATH})

find_package(MPI REQUIRED)	# MPI support
target_link_libraries(${PROJECT_NAME} ${MPI_C_LIBRARIES} openpa)

#if(FFTW_FOUND)
target_link_libraries(${PROJECT_NAME} ${FFTW_LIBRARIES})
target_include_directories(${PROJECT_NAME} PUBLIC ${FFTW_INCLUDES})
#endif()

install(TARGETS ${PROJECT_NAME} DESTINATION "lib")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h mpi/mdl.h mdlbase.h DESTINATION "include")
