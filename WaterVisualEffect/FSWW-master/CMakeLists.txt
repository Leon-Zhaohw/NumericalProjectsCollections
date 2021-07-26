cmake_minimum_required(VERSION 3.8)
project(wawa LANGUAGES CXX CUDA)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules/")

set(CMAKE_CXX_STANDARD 11)

# find Eigen
find_package(Eigen REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(CUDA QUIET REQUIRED)
set(CUDA_VERSION 9.2)

# find OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "-g -O3 -DDEBUG -fPIC ${CMAKE_CXX_FLAGS}")
set(CMAKE_CUDA_FLAGS "-D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED ${CMAKE_CUDA_FLAGS}")
#set(CUDA_NVCC_FLAGS "-Wno-deprecated-gpu-targets --expt-relaxed-constexpr -arch=sm_30  -D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED -D__STRICT_ANSI__")

# find Magnum
find_package(Magnum REQUIRED
  MeshTools
  Primitives
  Shaders
  SceneGraph
  Sdl2Application)
include_directories(${SDL2_INCLUDE_DIR})
include_directories(build/external/include)

# find ImGui
find_library(IMGUI_LIB ImGui)
find_library(MAGNUM_IMGUI_LIB MagnumImGui)

# libraries to link
set(GRAPHICS_LIBS
  Magnum::MeshTools
  Magnum::Primitives
  Magnum::Application
  Magnum::Magnum
  Magnum::Shaders
  Magnum::SceneGraph
  ${IMGUI_LIB}
  ${MAGNUM_IMGUI_LIB})

add_subdirectory(src)

# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to 
# particles will also build with -std=c++11
target_compile_features(wawa PUBLIC cxx_std_11)
 
# We need to explicitly state that we need all CUDA files in the 
# particle library to be built with -dc as the member functions 
# could be called by other libraries and executables
#set_target_properties( wawa
#                       PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

