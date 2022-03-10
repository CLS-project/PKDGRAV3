cmake_minimum_required(VERSION 3.0)
project(mdlmetal)

if (APPLE)
  find_library(METAL NAMES Metal)
  find_library(FOUNDATION NAMES Foundation)
  find_library(QUARTZCORE NAMES QuartzCore)
  if(METAL AND FOUNDATION AND QUARTZCORE)
    execute_process(
      COMMAND xcrun --find metal
      OUTPUT_VARIABLE _xcrun_out OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_VARIABLE _xcrun_err)
    if(_xcrun_out)
      set(metal_EXECUTABLE "${_xcrun_out}" CACHE FILEPATH "metal executable")
      mark_as_advanced(metal_EXECUTABLE)
    endif()
    execute_process(
      COMMAND xcrun --find metallib
      OUTPUT_VARIABLE _xcrun_out OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_VARIABLE _xcrun_err)
    if(_xcrun_out)
      set(metallib_EXECUTABLE "${_xcrun_out}" CACHE FILEPATH "metallib executable")
      mark_as_advanced(metallib_EXECUTABLE)
    endif()
    set(USE_METAL TRUE)
    set(USE_METAL ${USE_METAL} PARENT_SCOPE)
    add_library(${PROJECT_NAME} STATIC mdlmetal.mm)
    target_link_libraries(${PROJECT_NAME} ${METAL})
    target_link_libraries(${PROJECT_NAME} ${FOUNDATION})
    target_link_libraries(${PROJECT_NAME} ${QUARTZCORE})
  endif(METAL AND FOUNDATION AND QUARTZCORE)
endif(APPLE)
if(USE_METAL)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
  )
  target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/metal-cpp)
  target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
  target_link_libraries(${PROJECT_NAME} mdl2)
  install(TARGETS ${PROJECT_NAME} DESTINATION "lib")
else(USE_METAL)
  add_library(${PROJECT_NAME} INTERFACE)
endif(USE_METAL)


