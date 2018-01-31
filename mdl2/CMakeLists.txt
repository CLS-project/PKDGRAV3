cmake_minimum_required(VERSION 3.0)
project(mdl2)
find_package(MPI REQUIRED)

add_subdirectory(openpa)

INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILES (sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES (unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES (signal.h HAVE_SIGNAL_H)
CHECK_INCLUDE_FILES (inttypes.h HAVE_INTTYPES_H)

INCLUDE(CheckLibraryExists)
CHECK_LIBRARY_EXISTS(memkind hbw_posix_memalign "" HAVE_MEMKIND)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/mdl_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/mdl_config.h)

set(mdl_SOURCES mpi/mdl.c mdlbase.c)
if(APPLE)
set(mdl_SOURCES ${mdl_SOURCES} mac/pthread_barrier.c)
endif()

add_library(${PROJECT_NAME} STATIC ${mdl_SOURCES})
#target_link_libraries(${PROJECT_NAME} openpa)
if (APPLE)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mac)
endif()
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mpi)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${PROJECT_NAME} PRIVATE ${MPI_C_LIBRARIES} openpa)
target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_C_INCLUDE_PATH}) 
