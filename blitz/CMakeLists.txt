cmake_minimum_required(VERSION 3.12)
cmake_policy(VERSION 3.12...3.27)
project(blitz)


set(BLITZ_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(BLITZ_HEADERS
    ${BLITZ_INCLUDE_DIR}/blitz/array.h
)

add_library(blitz STATIC globals.cpp)
#target_include_directories(blitz INTERFACE $<BUILD_INTERFACE:${BLITZ_INCLUDE_DIR}>
#					   $<INSTALL_INTERFACE:include>)
target_include_directories(blitz INTERFACE ${BLITZ_INCLUDE_DIR})
target_include_directories(blitz PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/bzconfig.h.in ${CMAKE_CURRENT_BINARY_DIR}/bzconfig.h)
install(FILES ${BLITZ_HEADERS} DESTINATION "include/blitz")
