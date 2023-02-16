cmake_minimum_required(VERSION 3.0)
project(pkdmetal)

if (USE_METAL)
  add_library(${PROJECT_NAME} STATIC metal.mm)

  add_custom_command(
    MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/pppc.metal"
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/pppc.air"
    COMMAND ${metal_EXECUTABLE} -I "${pkdgrav3_SOURCE_DIR}" -I "${CMAKE_BINARY_DIR}"
      -c "${CMAKE_CURRENT_SOURCE_DIR}/pppc.metal" -o "${CMAKE_CURRENT_BINARY_DIR}/pppc.air"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Creating pppc.air"
    VERBATIM)
  add_custom_command(
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pppc.air
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pkdgrav3.metallib
    COMMAND ${metallib_EXECUTABLE}
      "${CMAKE_CURRENT_BINARY_DIR}/pppc.air"
      -o ${CMAKE_CURRENT_BINARY_DIR}/pkdgrav3.metallib
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    COMMENT "Creating pkdgrav3.metallib"
    VERBATIM)

  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pkdgrav3_metallib.h
    MAIN_DEPENDENCY ${CMAKE_CURRENT_BINARY_DIR}/pkdgrav3.metallib
    COMMAND ${PYTHON} ${CMAKE_CURRENT_SOURCE_DIR}/make_header.py ${CMAKE_CURRENT_BINARY_DIR}/pkdgrav3.metallib ${CMAKE_CURRENT_BINARY_DIR}/pkdgrav3_metallib.h
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    VERBATIM)
    add_custom_target(pkdgrav3_metallib ALL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pkdgrav3_metallib.h)
  add_dependencies(${PROJECT_NAME} pkdgrav3_metallib)
  target_link_libraries(${PROJECT_NAME} PUBLIC ${METAL_LIBRARIES} mdl2 blitz)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
  )
else(USE_METAL)
  add_library(${PROJECT_NAME} INTERFACE)
endif(USE_METAL)

#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/metal-cpp)
#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
#target_include_directories(fmt-header-only INTERFACE
#  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
#  $<INSTALL_INTERFACE:include>)


