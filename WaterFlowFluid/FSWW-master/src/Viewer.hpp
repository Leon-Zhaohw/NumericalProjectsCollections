#ifndef VIEWER_HPP
#define VIEWER_HPP

#include <Magnum/Buffer.h>
#include <Magnum/DefaultFramebuffer.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Mesh.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/MeshTools/CompressIndices.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Primitives/Plane.h>
#include <Magnum/Primitives/UVSphere.h>
#include <Magnum/Renderer.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/SceneGraph/Drawable.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Scene.h>
#include <Magnum/Shader.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Trade/MeshData3D.h>

#include <MagnumImGui.h>
#include <imgui.h>

#include "DrawingPrimitives.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

#include "definitions.hpp"
#include "WaterSurface.hpp"

namespace Magnum {

  using namespace Magnum::Math::Literals;
  using namespace definitions;

  typedef SceneGraph::Object<SceneGraph::MatrixTransformation3D> Object3D;
  typedef SceneGraph::Scene<SceneGraph::MatrixTransformation3D>  Scene3D;

  class Viewer : public Platform::Sdl2Application {
  public:
    explicit Viewer(const Arguments &arguments);

  private:
    void treatArguments(int argc, char **argv);
  
    void drawEvent() override;
    void drawGui();
    void update();

    void viewportEvent(const Vector2i &size) override;

    void keyPressEvent(KeyEvent &event) override;
    void keyReleaseEvent(KeyEvent &event) override;
    void mousePressEvent(MouseEvent &event) override;
    void mouseReleaseEvent(MouseEvent &event) override;
    void mouseMoveEvent(MouseMoveEvent &event) override;
    void mouseScrollEvent(MouseScrollEvent &event) override;
    void textInputEvent(TextInputEvent &event) override;

    void mouseTranslation(MouseMoveEvent const &event, Vector2 delta);
    void mouseRotation(MouseMoveEvent const &event, Vector2 delta);
    void mouseZoom(MouseMoveEvent const &event, Vector2 delta);
    FLOAT scriptedZoom(float delta);
    void mouseZoom(float delta);
    void mousePan(MouseMoveEvent const &event, Vector2 delta);

    Vector3 mousePlanePosition(Vector2i mouseScreenPos);
    void reset();
#ifdef PROJECTED_GRID
    void updateProjectedGrid();
#endif
    void exportMitsuba();
    void script();
    std::ofstream& scriptMitsuba(std::ofstream& file);
  
    Scene3D                     _scene;
    Object3D *                  _cameraObject;
    SceneGraph::Camera3D *      _camera;
    SceneGraph::DrawableGroup3D _drawables;

    Vector2i _previousMousePosition;
    Vector3  _mousePlanePosition[2];
    Vector3  _source;

    MagnumImGui _gui;

    DrawablePlane * plane;
    DrawablePlane * plane_obs;
    DrawableSphere *sphere;
    DrawableLine *  line;

    WaterSurface _surface;
    uint time = 0;
    uint stop_time = 1e4; // program stopped after <stop_time> frame

    FLOAT prev_solve_time;
    FLOAT prev_fps;
    bool up_cam;

    uint export_step_m;
    bool export_m;
    std::string export_file_m;
    std::string str_plot;
    std::ofstream stream_plot;
  
    bool zoom;
    Vector3 target_lookat;
  };

} // namespace Magnum

#endif
