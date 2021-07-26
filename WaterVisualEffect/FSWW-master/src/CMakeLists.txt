find_package(Magnum REQUIRED
  MeshTools
  Primitives
  Shaders
  SceneGraph
  Sdl2Application
  TextureTools
  Text)

add_executable(wawa main.cpp settings.cpp ui_parameters.cpp DrawingPrimitives.cpp Viewer.cpp definitions.cpp WaterSurface.cpp Times.cpp EquivalentSource.cpp LinearWave.cpp Obstacle.cpp InputPoint.cpp Grid.cpp CircularObstacle.cpp SquareObstacle.cpp TextureObstacle.cpp MovingEquivalentSource.cpp CudaWaterSurface.cu error.cpp Wave.cpp WavyObstacle.cpp ProjectedGrid.cpp)

target_link_libraries(wawa ${GRAPHICS_LIBS} Magnum::TextureTools Magnum::Text SDL2_image)

