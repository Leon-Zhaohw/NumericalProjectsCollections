set(Sources
  CGInterface.cpp
  Cube.cpp
  CubeMesh.cpp
  CubeNewtonSolver.cpp
  EigenCG.cpp
  Material.cpp
  QuasistaticSolveResult.cpp
  StableNeoHookean.cpp
  TetMesh.cpp
)

set(Headers
  CGInterface.h
  CGSolverType.h
  Cube.h
  CubeMesh.h
  CubeNewtonSolver.h
  EigenCG.h
  InitialMesh.h
  LineSearches.h
  LocalHessianLocation.h
  Material.h
  QuasistaticSolveResult.h
  Settings.h
  StableNeoHookean.h
  TetMesh.h
  TetNewtonSolver.h
)

add_library(cubesim ${Headers} ${Sources})

find_package(Eigen3 3.3.5 REQUIRED)
target_include_directories(cubesim SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})

if(USE_OPENMP)
  find_package(OpenMP)
  if(NOT OPENMP_FOUND)
    message(FATAL_ERROR "Error, failed to find OpenMP.")
  endif()
  target_compile_definitions(cubesim PUBLIC OPENMP_FOUND)
  target_compile_options(cubesim PRIVATE ${OpenMP_CXX_FLAGS})
  target_link_libraries(cubesim INTERFACE ${OpenMP_CXX_FLAGS})
endif()
