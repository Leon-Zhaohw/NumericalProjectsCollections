project(kelvinlet)

cmake_minimum_required(VERSION 2.6)
set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -Wno-deprecated-declarations")

option(USE_DEBUG "Debug?" OFF)
option(USE_TBB   "TBB?"   OFF)

if(USE_DEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else()
    set(CMAKE_BUILD_TYPE Release)
endif()

set(LINK_FLAGS -m64)
add_definitions(-m64)

set(KELVINLET_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

set(SRCS 
    ${KELVINLET_DIR}/kelvinlet/dynaBase.cpp
    ${KELVINLET_DIR}/kelvinlet/dynaPulseAffine.cpp
    ${KELVINLET_DIR}/kelvinlet/dynaPulseBase.cpp
    ${KELVINLET_DIR}/kelvinlet/dynaPulseGrab.cpp
    ${KELVINLET_DIR}/kelvinlet/dynaPushAffine.cpp
    ${KELVINLET_DIR}/kelvinlet/dynaPushBase.cpp
    ${KELVINLET_DIR}/kelvinlet/dynaPushGrab.cpp
    main.cpp)

if(APPLE)
    set(CMAKE_SHARED_LIBRARY_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CXX_FLAGS} -framework Accelerate")
endif()

find_package(Eigen3)

if(USE_TBB)
    find_package(TBB)
endif()

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${KELVINLET_DIR})

if (USE_TBB AND TBB_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_TBB")
    include_directories(${TBB_INCLUDE_DIRS})
    include_directories(${TBB_MALLOC_INCLUDE_DIRS})
endif()

add_executable(${PROJECT_NAME} ${SRCS})

if (USE_TBB AND TBB_FOUND)
    target_link_libraries(${PROJECT_NAME} ${TBB_LIBRARIES})
endif()
