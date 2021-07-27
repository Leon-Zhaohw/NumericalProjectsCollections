cmake_minimum_required(VERSION 2.8.3)
project(2D_fluid_Semi)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# ==============================================================================
# Additional cmake find modules
# ==============================================================================
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#########################################################
#Find FFTW
#########################################################
find_library(FFTW_LIBRARY
     NAMES fftw3 fftw)
#if (FFTW_LIBRARY_FOUND)
#  message("-- Found FFTW: ${FFTW_INCLUDE_DIRS}")
#else (FFTW_LIBRARY_FOUND)
#  message(FATAL_ERROR "Can't find FFTW.")
#endif (FFTW_LIBRARY_FOUND)

#find_library(FFTWF_LIBRARY
#     NAMES fftw3f fftwf)
     
find_library(FFTW_OMP_LIBRARY
    NAMES fftw3_omp)
#if (FFTWF_LIBRARY_FOUND)
#  message("-- Found FFTWF: ${FFTWF_INCLUDE_DIRS}")
#else (FFTWF_LIBRARY_FOUND)
#  message(FATAL_ERROR "Can't find FFTWF")
#endif (FFTWF_LIBRARY_FOUND)


#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR " GLUT not found!")
endif(NOT GLUT_FOUND)
#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR " OPENGL not found!")
endif(NOT OPENGL_FOUND)
#########################################################
# FIND OpenMP
#########################################################
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
#########################################################
# FIND X11
#########################################################
find_package(X11)
if (X11_FOUND)
  message("--Found X11")
endif(X11_FOUND)
message("${X11_LIBRARIES}  ${X11_INCLUDE_DIR}")
#########################################################
# FIND libjpeg
#########################################################
find_package(JPEG)
if (JPEG_FOUND)
  message("--Found jpeg.")
endif(JPEG_FOUND)

add_subdirectory(libraries)

include_directories(
  src
  libraries/glvu
  libraries/Eigen/
  ${EIGEN_INCLUDE_DIRS}
  ${FFTW_INCLUDE_DIRS}
  #${FFTWF_INCLUDE_DIRS}
  ${GLUT_INCLUDE_DIRS}
  ${OpenGL_INCLUDE_DIRS}
  ${X11_INCLUDE_DIR}
  ${stlplus3_SOURCE_DIR}
  ${JPEG_INCLUDE_DIR}
)

add_library(glvu_libs libraries/glvu/glvu.cpp)

add_subdirectory(src)
add_subdirectory(applications)
