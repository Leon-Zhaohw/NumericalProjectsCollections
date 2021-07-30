# Add headers
file(GLOB ALL_HDRS *.h)

set(ALL_SRC
    2D/FIELD2D.cpp
    2D/FLUID2D.cpp
    2D/obstacle_2d.cpp
    2D/VFIELD2D.cpp
    2D/laplacian_basis_2D.cpp
    2D/all_dirichlet_basis_2D.cpp
    2D/three_dirichlet_one_neumann.cpp
    2D/two_neumann_x_2D.cpp
    2D/basis_lookup_2D.cpp
    2D/laplacian_fluid_2D.cpp
    2D/drawer_2d.cpp
    2D/two_neumann_basis_2D.cpp
    3D/drawer_3d.cpp
    3D/laplacian_basis_set_3d.cpp
    3D/dirichlet_basis_3d.cpp
    3D/dirichlet_basis_set_3d.cpp
    3D/one_neumann_basis_3d.cpp
    3D/one_neumann_basis_set_3d.cpp
    3D/one_neumann_basis_set_3d_test.cpp
    3D/two_neumann_x_3d_basis.cpp
    3D/two_neumann_x_3d_basis_set.cpp
    3D/four_neumann_basis_3d.cpp
    3D/four_neumann_basis_set_3d.cpp
    3D/six_neumann_basis_3d.cpp
    3D/six_neumann_basis_set_3d.cpp
    3D/FIELD_3D.cpp
    3D/VECTOR3_FIELD_3D.cpp
    3D/laplacian_fluid_3d.cpp
    3D/FLUID_3D_MIC.cpp
    3D/trig_integral_3d.cpp
    3D/obstacle_box_3d.cpp
    3D/obstacle_wrapper_3d.cpp
    3D/obstacle_sphere_3d.cpp
    3D/obstacle_cylinder_3d.cpp
    3D/obstacle_cylinder_static_3d.cpp
    3D/obstacle_cylinder_interact.cpp
    3D/density_warpper.cpp
    3D/fluid_3D_DCT.cpp
    Alg/VECTOR.cpp
    Alg/MATRIX3.cpp
    solver/integrator_RK4.cpp
    solver/integrator_semi_implicit.cpp
    solver/obstacle_solver.cpp
    solver/integrator_2d.cpp
    solver/trapezoidal.cpp
    3D/particle_3d.cpp
    util/block_3d_dct.cpp
    util/trackball.cpp
    util/util.cpp
    util/timer.cpp
    util/read_write_tensor.cpp
    util/write_density_pbrt.cpp
    util/stringprintf.cc
    util/glPrintString.cpp
    util/solve_quadratic.cpp
    util/SIMPLE_PARSER.cpp)

set(ALL_SOURCE
  ${ALL_SRC}
  ${ALL_HDRS})

set(ALL_DEPENDICIES
  ${FFTW_OMP_LIBRARY}
  ${FFTW_LIBRARY}
  #${FFTWF_LIBRARY}
  ${OPENGL_LIBRARIES}
  ${GLUT_LIBRARY}
  ${X11_LIBRARIES}
  stlplus3
  ${JPEG_LIBRARIES}
  glvu_libs
  )

add_library(fluid ${ALL_SOURCE})
target_link_libraries(fluid ${ALL_DEPENDICIES})

#TEST CASES.
#add_executable(all_dirichlet_basis_2D_test 2D/all_dirichlet_basis_2D_test.cpp)
#target_link_libraries(all_dirichlet_basis_2D_test fluid ${ALL_DEPENDICIES})

#add_executable(three_dirichlet_one_neumann_test 2D/three_dirichlet_one_neumann_test.cpp)
#target_link_libraries(three_dirichlet_one_neumann_test fluid ${ALL_DEPENDICIES})

add_executable(two_neumann_x_2D_test 2D/two_neumann_x_2D_test.cpp)
target_link_libraries(two_neumann_x_2D_test fluid ${ALL_DEPENDICIES})

#add_executable(obstacle_2d_test 2D/obstacle_2d_test.cpp)
#target_link_libraries(obstacle_2d_test fluid ${ALL_DEPENDICIES})

add_executable(dirichlet_basis_set_3d_test 3D/dirichlet_basis_set_3d_test.cpp)
target_link_libraries(dirichlet_basis_set_3d_test fluid ${ALL_DEPENDICIES})

add_executable(one_neumann_basis_set_3d_test 3D/one_neumann_basis_set_3d_test.cpp)
target_link_libraries(one_neumann_basis_set_3d_test fluid ${ALL_DEPENDICIES})

add_executable(two_neumann_x_3d_basis_set_test 3D/two_neumann_x_3d_basis_set_test.cpp)
target_link_libraries(two_neumann_x_3d_basis_set_test fluid ${ALL_DEPENDICIES})

add_executable(four_neumann_basis_set_3d_test 3D/four_neumann_basis_set_3d_test.cpp)
target_link_libraries(four_neumann_basis_set_3d_test fluid ${ALL_DEPENDICIES})

add_executable(six_neumann_basis_set_3d_test 3D/six_neumann_basis_set_3d_test.cpp)
target_link_libraries(six_neumann_basis_set_3d_test fluid ${ALL_DEPENDICIES})

#add_executable(trig_integral_3d_test 3D/trig_integral_3d_test.cpp)
#target_link_libraries(trig_integral_3d_test fluid ${ALL_DEPENDICIES})
