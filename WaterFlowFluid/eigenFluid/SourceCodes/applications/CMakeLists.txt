
add_executable(laplacian_fluid_sim2D laplacian_fluid_sim2D.cpp)
target_link_libraries(laplacian_fluid_sim2D fluid ${ALL_DEPENDICIES})

add_executable(laplacian_fluid_sim3D laplacian_fluid_sim3D.cpp)
target_link_libraries(laplacian_fluid_sim3D fluid ${ALL_DEPENDICIES})

add_executable(precompute_2D_tensor precompute_2D_tensor.cpp)
target_link_libraries(precompute_2D_tensor fluid ${ALL_DEPENDICIES})

add_executable(precompute_3D_tensor precompute_3D_tensor.cpp)
target_link_libraries(precompute_3D_tensor fluid ${ALL_DEPENDICIES})

add_executable(drop_tensor_entries drop_tensor_entries.cpp)
target_link_libraries(drop_tensor_entries fluid ${ALL_DEPENDICIES})

add_executable(fluidStam3D fluidStam3D.cpp)
target_link_libraries(fluidStam3D fluid ${ALL_DEPENDICIES})

add_executable(fluid_DCT_sim3D fluid_DCT_sim3D.cpp)
target_link_libraries(fluid_DCT_sim3D fluid ${ALL_DEPENDICIES})
