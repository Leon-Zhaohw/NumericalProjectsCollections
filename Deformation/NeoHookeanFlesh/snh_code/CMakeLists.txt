cmake_minimum_required(VERSION 3.6.1)

execute_process(COMMAND ./eigen_found.sh RESULT_VARIABLE EIGEN_NOT_FOUND WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
if( EIGEN_NOT_FOUND )
  message(FATAL_ERROR "Failed to locate Eigen or Eigen's checksum is incorrect. Please run download_eigen.sh.")
endif()

if(NOT CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
else(NOT CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE STREQUAL "")

set(VALID_BUILD_TYPES None Debug Release RelWithDebInfo MinSizeRel)
if(NOT ";${VALID_BUILD_TYPES};" MATCHES ";${CMAKE_BUILD_TYPE};")
  message(FATAL_ERROR "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif()

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

option(USE_OPENMP "Use OpenMP for multi-threaded execution." OFF)

add_subdirectory(cubesim)
add_subdirectory(hexcli)
add_subdirectory(tetcli)
