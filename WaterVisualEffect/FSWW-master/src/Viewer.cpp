#include "Viewer.hpp"
#include "Times.hpp"
#include <SDL2/SDL.h>
#include <SDL2/SDL_image.h>
#include "ui_parameters.hpp"
#include "settings.hpp"
#include "error.hpp"
#include <fstream>

using namespace Magnum;
using namespace ui_parameters;
using namespace settings;


void help() {
  std::cout<<"\n     *** WAVE: Help ***\n"<<std::endl;
  std::cout<<"Synopsis: \n     .\\main <options>\n\nOptions:"<<std::endl;
  std::cout<<"     -l, -load <file>: load configuration file"<<std::endl;
  std::cout<<"     -e, -export <name>: export amplitude grid for each frequencies in the file <name>.obj"<<std::endl;
  std::cout<<"     -i, -import <name>: import amplitude grid for each frequencies in the file <name>.obj"<<std::endl;
  std::cout<<"     -stop <t>: stop animation and exit at time t"<<std::endl;
  std::cout<<"     -em <name>: export heightfields and render files in a set of files <name><frame number>.ong and <name><frame number>.xml"<<std::endl;
  std::cout<<"     -es, -export_step <n>: export every n frames"<<std::endl;
  std::cout<<"     -h, -help: print help\n"<<std::endl;
  exit(0);
}

void Viewer::treatArguments(int argc, char **argv) {
  for (int i = 1;  i < argc; ++i) {
    std::string s(argv[i]);
    if (s == "-l" || s == "-load") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help();
      }
      std::cout<<"Loading configuration file:"<<" "<<argv[i+1]<<std::endl;
      _surface.setImportConf(argv[i+1]);
      ++i;
    } else if (s == "-i" || s == "-import") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help();
      }
      std::cout<<"Importing"<<" "<<argv[i+1]<<std::endl;
      _surface.setImport(argv[i+1]);
      ++i;
    } else if (s == "-e" || s == "-export") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help();
      }
      std::cout<<"Exporting"<<" "<<argv[i+1]<<std::endl;
      _surface.setExport(argv[i+1]);
      ++i;
    } else if (s == "-em") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help();
      }
      std::cout<<"Exporting (mitsuba) "<<" "<<argv[i+1]<<std::endl;
      export_m = true;
      export_file_m = argv[i+1];
      std::stringstream ss_plot;
      ss_plot <<export_file_m<<"_plot.txt";
      str_plot = std::string(ss_plot.str());
      ++i;
    } else if (s == "-d" || s == "-data") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help();
      }
      std::cout<<"Saving amplitude data in"<<" "<<argv[i+1]<<std::endl;
      _surface.setData(argv[i+1]);
      ++i;
    } else if (s == "-stop") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help();
      }
      std::cout<<"Stop at t = "<<argv[i+1]<<std::endl;
      // _surface.setStopTime(atoi(argv[i+1]));
      stop_time = atoi(argv[i+1]);
      ++i;
    } else if (s == "-export_step" || s == "-es") {
      if (argc < i + 2) {
  	std::cerr<<"\nERROR: wrong number of arguments\n"<<std::endl;
  	help();
      }
      std::cout<<"Export every "<<argv[i+1]<<" steps"<<std::endl;
      _surface.setExportStep(atoi(argv[i+1]));
      export_step_m = atoi(argv[i+1]);
      ++i;
    } else if (s == "-r" || s == "-run") {
      running_ = true;
    } else if (s == "-h" || s == "-help") {
      help();
    } else {
      std::cerr<<"\nERROR: Unknown option\n"<<std::endl;
      help();
    }
  }
}

Viewer::Viewer(const Arguments &arguments)
  : Platform::Application{arguments,
    Configuration{}
  .setTitle("Water Wave")
     .setWindowFlags(Sdl2Application::Configuration::
		     WindowFlag::Maximized)},
		 _surface() {
		   export_step_m = 1;
		   export_m = false;
		   zoom = false;
		   treatArguments(arguments.argc, arguments.argv);
		   _surface.reset();
     
		   Renderer::enable(Renderer::Feature::DepthTest);
		   Renderer::enable(Renderer::Feature::FaceCulling);

		   Renderer::setPointSize(12.0);
		   Renderer::setLineWidth(6.0);

		   /* Configure camera */
		   _cameraObject = new Object3D{&_scene};
		   _cameraObject->translate(Vector3::zAxis(4.0f)).rotateX(Rad{M_PI / 4});
		   _camera = new SceneGraph::Camera3D{*_cameraObject};
		   viewportEvent(defaultFramebuffer.viewport().size()); // set up camera

		   plane = (new DrawablePlane(&_scene, &_drawables, n_rows_ - 1, n_cols_ - 1));
		   plane_obs =
		     (new DrawablePlane(&_scene, &_drawables, n_rows_obs - 1, n_cols_obs - 1));
		   line = (new DrawableLine(&_scene, &_drawables, 2 * n_rows_ * n_cols_));

		   /* draw sources and boundary points as spheres, */ 
		   /* do not use if too many sources and bp */
		   // std::list<VEC2> boundaries;
		   // _surface.getObstacleBoundary(boundaries);
		   // for (it = boundaries.begin(); it != boundaries.end(); ++it) {
		   //   sphere = (new DrawableSphere(&_scene, &_drawables, 2, 3));
		   //   sphere->setVertices([this](int i, DrawableSphere::VertexData &v) {
	
		   // 	 v.position *= 0.005;
		   // 	 v.position += Vector3((*it)[0], (*it)[1], 0.021);
		   // 	 v.color = Color4{1.0, 1.0, 1., 1.};
		   //     });
		   // }

    		   plane_obs->setVertices([this](int i, DrawableMesh::VertexData &v) {
    		       int ix        = i / n_cols_obs;
    		       int iy        = i % n_cols_obs;
    		       v.position[2] = _surface.height_obs(ix, n_cols_obs-iy-1);
    		       VEC2 p = gridObs2viewer(ix, n_cols_obs-iy-1);
    		       v.position[0] = p(0);
    		       v.position[1] = p(1);
		     });
		   
		   Times::TIMES->init();
		   time = 0;
		   prev_solve_time = 0;
		   prev_fps = 0;
		   setMinimalLoopPeriod(50);
		   up_cam = false;
#ifdef PROJECTED_GRID
		   Vector2i ori = {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
		   Vector3 mpos3 = mousePlanePosition(ori);
		   target_lookat = mpos3;
		   _surface.setTargetLookAt(VEC2(target_lookat[0],target_lookat[1]));
		   updateProjectedGrid();
		   
#endif
#ifdef PLOT_RESULT
		   if (export_m) {
		     if (stream_plot.is_open()) {
		       stream_plot.close();
		     }
		     INFO("Writing "<<str_plot);
		     stream_plot.open(str_plot);
		     stream_plot<<"set view map\n";
		     stream_plot<<"unset key\n";
		     stream_plot<<"unset tics\n";
		     stream_plot<<"unset border\n";
		     stream_plot<<"unset colorbox\n";
		     stream_plot<<"set terminal png size 600, 600\n";
		   }
#endif
		 }

void Viewer::update() {
  try {
    _surface.update();

  } catch (std::exception& e) {
    std::cerr << "Exception catched : " << e.what() << std::endl;
    _surface.clear();
    throw;
  }
}

void Viewer::drawEvent() {
  Times::TIMES->tick(Times::total_time_);
  defaultFramebuffer.clear(FramebufferClear::Color | FramebufferClear::Depth);
  Times::TIMES->tick(Times::simu_time_);
  //#ifndef INTERACTIVE_
  // if (zoom && running_) {
  //   FLOAT d = scriptedZoom(0.00125);
  //   if (d < 0.4) {
  //     zoom = false;
  //   }
  // } else if (time == 12) {
  //   zoom = true;
  // }
  //#endif
  if (running_) {
    script();
  }

  if (up_cam) {
    up_cam = false;
  
#ifdef PROJECTED_GRID
    Vector2i ori =
      {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
    Vector3 mpos3 = mousePlanePosition(ori);
    target_lookat = mpos3;
    _surface.setTargetLookAt(VEC2(target_lookat[0],target_lookat[1]));
    updateProjectedGrid();
#endif
  }

  if (running_) {
    update();
    ++time;
  }
  if (time >= stop_time) {
    if (stream_plot.is_open()) {
      stream_plot<<"set term pop; set out;";
      stream_plot.close();
    }
    exit();
  }
   
  if (running_) {
    if (export_m) {
      exportMitsuba();
    }
  }

  Times::TIMES->tock(Times::simu_time_);
  
  Times::TIMES->tick(Times::display_time_);

  // Set up wave grid visualization
  plane->setVertices([this](int i, DrawableMesh::VertexData &v) {
   
#ifdef PROJECTED_GRID
      VEC3 p = _surface.getPosProjGrid(i);
      if (isnan(p(2)) || isinf(p(2))) {
	WARNING(!(isnan(p(2)) || isinf(p(2))), "undefined height value", "[Viewer::drawEvent] "<<p(2));
	p(2) = 0;
      }
      v.position = Vector3{p(0), p(1), p(2)};
#else
      int ix = i / n_cols_;
      int iy = n_cols_ - i % n_cols_;
      VEC3 p = _surface.getPosGrid(ix, iy);
      v.position = Vector3{p(0), p(1), p(2)};
#endif
      // v.color = Color4::fromHsv(Rad{10.f * M_PI * v.position[2]}, 0.7, 1.0);
    });

  plane_obs->setVertices([this](int i, DrawableMesh::VertexData &v) {
      int ix        = i / n_cols_obs;
      int iy        = i % n_cols_obs;
      v.position[2] = _surface.height_obs(ix, n_cols_obs-iy-1);
      VEC2 p = gridObs2viewer(ix, n_cols_obs-iy-1);
      v.position[0] = p(0);
      v.position[1] = p(1);
      // v.color       = Color4::fromHsv(Rad{10.f * M_PI * v.position[2]}, 0.7, 1.0);
    });
  

  line->_mesh.setPrimitive(MeshPrimitive::Lines);

  line->setVertices([this](int i, DrawableLine::VertexData &v) {
      int   ix  = (i / 2) / n_cols_;
      int   iy  = (i / 2) % n_cols_;
#ifdef PROJECTED_GRID
      VEC3 p = _surface.getPosProjGrid(ix, iy);
      v.position = Vector3{p(0), p(1), p(2)};
#else
      VEC3 p = _surface.getPosGrid(ix, iy);
      v.position = Vector3{p(0), p(1), p(2)};
#endif
      if (ix == n_rows_ / 2 && iy == n_cols_ / 2)
	v.color = Color4{1.0, 0., 0., 1.};
    });
      
  _camera->draw(_drawables);

  Times::TIMES->tock(Times::display_time_);
  Times::TIMES->tock(Times::total_time_);

  // INFO("Display time: "<<Times::TIMES->getTime(Times::display_time_));
  // INFO("Simu time: "<<Times::TIMES->getTime(Times::simu_time_));
  // INFO("Solve time: "<<Times::TIMES->getTime(Times::solve_time_)<<" --   Sum time: "<<Times::TIMES->getTime(Times::sum_up_time_));
  
  drawGui();
  swapBuffers();

  Times::TIMES->next_loop();
}

void Viewer::drawGui() {
  _gui.newFrame(windowSize(), defaultFramebuffer.viewport().size());
  ImGui::SetWindowFontScale(2);
  if (time%ampli_steps[0]-1 == 0) {
    prev_solve_time = Times::TIMES_UP->getTime(Times::solve_time_);
    prev_fps = 1.0/Times::TIMES->getTime(Times::total_time_);
  }
  
  std::stringstream ss;
  ss<<"FPS: "<<1.0/Times::TIMES->getTime(Times::total_time_)<<" "<<prev_fps;
  std::string str(ss.str());
  ImGui::Text(str.c_str());

  ss = std::stringstream();
  ss<<"Display time: "<<Times::TIMES->getTime(Times::display_time_)<<"s";
  str = ss.str();
  ImGui::Text(str.c_str());

  ss = std::stringstream();
  ss<<"Simu time: "<<Times::TIMES->getTime(Times::simu_time_)<<"s";
  str = ss.str();
  ImGui::Text(str.c_str());

  ss = std::stringstream();
  
  ss<<"Solve time (every "<<ampli_steps[0]<<" frame): "<<Times::TIMES_UP->getTime(Times::solve_time_)<<"s";
  str = ss.str();
  ImGui::Text(str.c_str());

  ss = std::stringstream();
  ss<<"Sum time: "<<Times::TIMES->getTime(Times::sum_up_time_)<<"s";
  str = ss.str();
  ImGui::Text(str.c_str());
  
  if (ImGui::Button("Reset")) {
    reset();
  }
  ImGui::Checkbox("Show In Field", &show_in_field);
  ImGui::Checkbox("Show Scattered Field", &show_scattered_field);

  ImGui::SliderFloat("Time step", &dt_, 0.01, 0.1);
  ImGui::SliderFloat("Offset Surface", &offset_, 0.1, 0.9);
  ImGui::SliderFloat("Step Sampling", &step_sampling_, 0.1, 1);

  ImGui::Separator();
  
  ImGui::SliderFloat("Amplitude", &height_ampli_, 0.01, 5);

  ImGui::Separator();
  ImGui::SliderFloat("Damping", &damping_, 0.0, 0.1);
  ImGui::Separator();
  
  ImGui::Text("Obstacles:");
  ImGui::Checkbox("Circle", &circle_);
  ImGui::Checkbox("Square", &square_);
  ImGui::Checkbox("Line", &line_);
  ImGui::Checkbox("Harbour", &harbour_);
  ImGui::Checkbox("Island", &island_);
  ImGui::Checkbox("Test", &test_);
  ImGui::Checkbox("Wavy", &wavy_);

  ImGui::Separator();
  
  ImGui::Text("Wave:");
  ImGui::Checkbox("Point Source (1, 15)", &point_source1_);
  ImGui::Checkbox("Point Source (5, 5)", &point_source2_);
  ImGui::Checkbox("Point Source (15, 15)", &point_source3_);
  ImGui::Checkbox("Linear Wave (1, 0)", &linear_wave1_);
  ImGui::Checkbox("Linear Wave (1, 1)", &linear_wave2_);
  ImGui::Checkbox("Linear Wave (0, 1)", &linear_wave3_);
  ImGui::Checkbox("User Defined Source", &user_def_source_);

  ImGui::Separator();

  ImGui::Checkbox("Neumann condition", &neumann);
 
  ImGui::Separator();
  
  //  ImGui::Checkbox("Gersner", &use_gersner);
  
  ImGui::Separator();

  ImGui::Checkbox("Random", &random_);

  ImGui::Separator();
  ImGui::SliderFloat("Drop Size", &size_drop_, _surface.minWL(), _surface.maxWL());
    
  _gui.drawFrame();

  redraw();
}

void Viewer::viewportEvent(const Vector2i &size) {
  defaultFramebuffer.setViewport({{}, size});

  _camera->setProjectionMatrix
    (Matrix4::perspectiveProjection(60.0_degf, Vector2{size}.aspectRatio(), 0.001f, 10000.0f));
}

void Viewer::keyPressEvent(KeyEvent &event) {
  if (_gui.keyPressEvent(event)) {
    redraw();
    return;
  }

  if (event.key() == KeyEvent::Key::Esc) {
    if (stream_plot.is_open()) {
      stream_plot<<"set term pop; set out;";
      stream_plot.close();
    }
    exit();
  }

  if (event.key() == KeyEvent::Key::Space) {
    user_def_pos_ = VEC2(_mousePlanePosition[0][0], _mousePlanePosition[0][1]);
  }

  if (event.key() == KeyEvent::Key::G) {
    use_gersner = !use_gersner;
  }

  if (event.key() == KeyEvent::Key::B) {
    createTabs();
  }

  if (event.key() == KeyEvent::Key::A) {
    approx_inter = !approx_inter;
  }
  if (event.key() == KeyEvent::Key::E) {
    _surface.drawHeighField("heightfiled.txt");
    _surface.exportAmplitude("ampli.txt");
    _surface.exportDirectivityAmplitude("ampli_dir.txt", 5);
    _surface.exportDirectivityAnalytic("ampli_dir_analytic.txt", 1, 5);
    _surface.exportSurfaceObj("test.obj");
  }

  if (event.key() == KeyEvent::Key::I) {
    show_in_field = !show_in_field;
  }
  if (event.key() == KeyEvent::Key::S) {
    show_scattered_field = !show_scattered_field;
  }

  if (event.key() == KeyEvent::Key::R) {
    reset();
  }

  if (event.key() == KeyEvent::Key::N) {
    neumann = !neumann;
  }

  if (event.key() == KeyEvent::Key::Enter) {
    running_ = !running_;
  }
  
  if (event.key() == KeyEvent::Key::Space) {
    if (!running_) {
      update();
    }
  }

  if (event.key() == KeyEvent::Key::W) {

    if (plane->_mesh.primitive() == MeshPrimitive::Points) {
      plane->_mesh.setPrimitive(MeshPrimitive::Triangles);
    } else {
      plane->_mesh.setPrimitive(MeshPrimitive::Points);
    }
  }

  redraw();
}

void Viewer::keyReleaseEvent(KeyEvent &event) {
  if (_gui.keyReleaseEvent(event)) {
    redraw();
    return;
  }
}

void Viewer::mousePressEvent(MouseEvent &event) {
  if (_gui.mousePressEvent(event)) {
    redraw();
    return;
  }
  if (event.button() == MouseEvent::Button::Middle) {
    Vector3 camOrg =
      _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
    Vector2i ori = {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
    Vector3 mpos3 = mousePlanePosition(ori);
    Vector3 dir = (mpos3 - camOrg).normalized();
    
    FLOAT lambda = 0.5;
    _mousePlanePosition[1] =
      lambda * _mousePlanePosition[0] + (1 - lambda) * _mousePlanePosition[1];
    _mousePlanePosition[0] = mousePlanePosition(event.position());
    _surface.addDrop(_mousePlanePosition[0][0], _mousePlanePosition[0][1], size_drop_);
  }

  _previousMousePosition = event.position();
  event.setAccepted();
}

void Viewer::mouseReleaseEvent(MouseEvent &event) {
  if (_gui.mouseReleaseEvent(event)) {
    redraw();
    return;
  }

  event.setAccepted();
  redraw();
}

void Viewer::mouseMoveEvent(MouseMoveEvent &event) {
  if (_gui.mouseMoveEvent(event)) {
    redraw();
    return;
  }

  FLOAT lambda = 0.5;
  _mousePlanePosition[1] =
    lambda * _mousePlanePosition[0] + (1 - lambda) * _mousePlanePosition[1];
  _mousePlanePosition[0] = mousePlanePosition(event.position());

  const Vector2 delta = Vector2{event.position() - _previousMousePosition} /
    Vector2{defaultFramebuffer.viewport().size()};

    if (event.buttons() & MouseMoveEvent::Button::Left)
      mouseRotation(event, delta);

    if (event.buttons() & MouseMoveEvent::Button::Right)
      mouseTranslation(event, delta);

    if (event.buttons() & MouseMoveEvent::Button::Middle)
      _surface.applyForce(_mousePlanePosition[0][0], _mousePlanePosition[0][1], 0.02);
    mousePan(event, delta);

    _previousMousePosition = event.position();
    event.setAccepted();
    redraw();
}

void Viewer::mouseScrollEvent(MouseScrollEvent &event) {
  if (_gui.mouseScrollEvent(event)) {
    redraw();
    return;
  }

  Vector2 mpos = Vector2{0, 0};
  Vector3 mpos3 = {mpos[0], mpos[1], -1.f};
  auto    trans =
    _cameraObject->transformation() * _camera->projectionMatrix().inverted();
  mpos3 = trans.transformPoint(mpos3);
  Vector3 camOrg =
    _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  Vector3 dir    = mpos3 - camOrg;
  FLOAT   lambda = -camOrg.z() / dir.z();
  if (lambda < 0) {
    lambda = 0.5;
  }
  FLOAT zcoef = lambda*dir.length();
  mouseZoom(-0.025f*zcoef*event.offset().y());

  up_cam = true;
}

void Viewer::textInputEvent(TextInputEvent &event) {
  if (_gui.textInputEvent(event)) {
    redraw();
    return;
  }
}

void Viewer::mouseRotation(MouseMoveEvent const &event,
			   Vector2               delta) {

  Vector2 mpos = Vector2{0, 0};
  Vector2i ori = {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
  Vector3 mpos3 = mousePlanePosition(ori);
  Vector3 camOrg =
    _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  auto right =
    (_cameraObject->transformation().transformVector(Vector3{1.0, 0.0, 0.0})).normalized();
  
  auto camPos =
    _cameraObject->transformation().transformPoint(Vector3{0.0, 0.0, 0.0});

  auto axis = cross(Vector3{0.f, 0.f, 1.f}, camPos.normalized()).normalized();

  _cameraObject->rotate(Rad{-3.0f * delta.y()}, right);
  _cameraObject->rotateZ(Rad{-3.0f * delta.x()});

  up_cam = true;
}

void Viewer::mouseTranslation(MouseMoveEvent const &event,
			      Vector2               delta) {

  Vector2 mpos = Vector2{0, 0};
  Vector2i ori = {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
  Vector3 mpos3 = mousePlanePosition(ori);
  Vector3 camOrg =
    _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  Vector3 dir    = mpos3 - camOrg;
  FLOAT coef = dir.length();
  
  auto camPos =
    _cameraObject->transformation().transformPoint(Vector3{0.0, 0.0, 0.0});
  auto dir1 =
    _cameraObject->transformation().transformVector(Vector3{-1.0, 0.0, 0.0});
  auto dir2 =
    _cameraObject->transformation().transformVector(Vector3{0.0, 1.0, 0.0});
  _cameraObject->translate(coef*delta.x() * dir1);
  _cameraObject->translate(coef*delta.y() * dir2);

  up_cam = true;
}

FLOAT Viewer::scriptedZoom(float delta) {
  Vector2 mpos = Vector2{0, 0};
  Vector2i ori = {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
  Vector3 mpos3 = mousePlanePosition(ori);
  Vector3 camOrg =
    _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  Vector3 dir    = target_lookat - camOrg;
  _cameraObject->translate(delta * dir);
  up_cam = true;
  return dir.length();
}

void Viewer::mouseZoom(float delta) {
  auto dir =
    _cameraObject->transformation().transformVector(Vector3{0.0, 0.0, 1.0});

  _cameraObject->translate(10.0f * delta * dir);
}

void Viewer::mouseZoom(MouseMoveEvent const &event, Vector2 delta) {
  auto dir =
    _cameraObject->transformation().transformVector(Vector3{0.0, 0.0, 1.0});

  _cameraObject->translate(10.0f * delta.y() * dir);
}

void Viewer::mousePan(MouseMoveEvent const &event, Vector2 delta) {}

Vector3 Viewer::mousePlanePosition(Vector2i mouseScreenPos) {
  Vector2 mpos = 2.0f * (Vector2{mouseScreenPos} /
			 Vector2{defaultFramebuffer.viewport().size()} -
                         Vector2{.5f, 0.5f});
  mpos[1] *= -1.f;

  Vector3 mpos3 = {mpos[0], mpos[1], -1.f};
  auto    trans =
    _cameraObject->transformation() * _camera->projectionMatrix().inverted();
  mpos3 = trans.transformPoint(mpos3);
  Vector3 camOrg =
    _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  Vector3 dir    = mpos3 - camOrg;
  FLOAT   lambda = -camOrg.z() / dir.z();
  if (lambda < 0) {
    lambda = 0;
  }
  mpos3          = camOrg + lambda * dir;

  return mpos3;
}



void Viewer::reset() {
  _surface.reset();

  _drawables = SceneGraph::DrawableGroup3D();

  plane = (new DrawablePlane(&_scene, &_drawables, n_rows_ - 1, n_cols_ - 1));

  plane_obs = (new DrawablePlane(&_scene, &_drawables, n_rows_obs - 1, n_cols_obs - 1));
  line = (new DrawableLine(&_scene, &_drawables, 2 * n_rows_ * n_cols_));

  /* draw sources and boundary points as spheres, */ 
  /* do not use if too many sources and bp */
  // std::list<VEC2> boundaries;
  // _surface.getObstacleBoundary(boundaries);
  // for (it = boundaries.begin(); it != boundaries.end(); ++it) {
  //   sphere = (new DrawableSphere(&_scene, &_drawables, 2, 3));
  //   sphere->setVertices([this](int i, DrawableSphere::VertexData &v) {
	
  // 	v.position *= 0.005;
  // 	v.position += Vector3((*it)[0], (*it)[1], 0.02);
  // 	v.color = Color4{1.0, 1.0, 1., 1.};
  //     });
  //  }

  plane_obs->setVertices([this](int i, DrawableMesh::VertexData &v) {
      int ix        = i / n_cols_obs;
      int iy        = i % n_cols_obs;
      v.position[2] = _surface.height_obs(ix, n_cols_obs-iy-1);
      VEC2 p = gridObs2viewer(ix, n_cols_obs-iy-1);
      v.position[0] = p(0);
      v.position[1] = p(1);
    });

  time = 0;
  up_cam = false;
#ifdef PROJECTED_GRID
  Vector2i ori = {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
  Vector3 mpos3 = mousePlanePosition(ori);
  target_lookat = mpos3;
  _surface.setTargetLookAt(VEC2(target_lookat[0],target_lookat[1]));
  updateProjectedGrid();
#endif
#ifdef PLOT_RESULT
  if (export_m) {
    if (stream_plot.is_open()) {
      stream_plot.close();
    }
    INFO("Writing "<<str_plot);
    stream_plot.open(str_plot);
    stream_plot<<"set view map\n";
    stream_plot<<"unset key\n";
    stream_plot<<"unset tics\n";
    stream_plot<<"unset border\n";
    stream_plot<<"unset colorbox\n";
    stream_plot<<"set terminal png size 600, 600\n";
  }
#endif

}

#ifdef PROJECTED_GRID
void Viewer::updateProjectedGrid() {
#pragma omp parallel for collapse(2)
  for (uint i = 0; i < n_rows_; ++i) {
    for (uint j = 0; j < n_cols_; ++j) {
      uint ix = i*defaultFramebuffer.viewport().sizeX()/n_rows_;
      uint jx = j*defaultFramebuffer.viewport().sizeY()/n_cols_;
      Vector3 proj = mousePlanePosition(Vector2i(ix, jx));
      _surface.setProjGrid(i, j, proj[0], proj[1]);
    }
  }
  _surface.updatePosGridCuda();

}

#endif


void Viewer::exportMitsuba() {
  Vector2i ori =
    {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
  Vector3 mpos3 = mousePlanePosition(ori);
  target_lookat = mpos3;
  Vector3 camOrg =
    _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  Vector3 dir    = target_lookat - camOrg;
  Vector3 left = cross(Vector3{0.f, 0.f, 1.f}, dir.normalized()).normalized();
  Vector3 axis = Vector3{0.f, 1.f, 0.f};//cross(dir.normalized(), left).normalized();
  
  std::stringstream ss_xml;
  uint n = time/export_step_m;
  std::string s0 = "";
  if (n < 10) {
    s0 = "000";
  } else if (n < 100) {
    s0 = "00";
  } else if (n < 1000) {
    s0 = "0";
  }
  
  ss_xml <<export_file_m<<s0<<n<<".xml";
  std::string str_xml(ss_xml.str());
  std::stringstream ss_dat;
#ifdef PLOT_RESULT
  ss_dat <<export_file_m<<s0<<n<<".dat";
  std::string str_dat(ss_dat.str());
  _surface.drawHeighField(str_dat);
  stream_plot<<"set output \""<<export_file_m<<"_2d_"<<s0<<n<<".png\"\n";
  stream_plot<<"splot '"<<str_dat<<"' with pm3d\n";
#endif
  std::stringstream ss_obj; 
  ss_obj <<export_file_m<<n<<".obj";
  std::string str_obj(ss_obj.str());
  _surface.exportSurfaceObj(str_obj);
  std::ofstream file;
  file.open(str_xml);
  file<<"<?xml  version=\"1.0\"  encoding=\"utf-8\"?>\n";
  file<<"<scene version=\"0.5.0\">\n";
  file<<"<integrator type=\"path\">\n";
  file<<"<integer name=\"maxDepth\" value=\"2\"/>\n";
  file<<"</integrator>\n";
  file<<"<shape type=\"obj\">\n";
  file<<"<string name=\"filename\" value=\""<<str_obj<<"\"/>\n";
  file<<"<bsdf type=\"plastic\">\n";
  file<<"<srgb name=\"diffuseReflectance\" value=\"#003862\"/>\n";
  file<<"<float name=\"intIOR\" value=\"3\"/>\n";
  file<<"<string name=\"extIOR\" value=\"air\"/>\n";
  file<<"</bsdf>\n";
  file<<"</shape>\n";

  _surface.exportObstacleMitsuba(file);
  scriptMitsuba(file);
  
  file<<"<sensor type=\"perspective\">\n";
  file<<"<float name=\"focusDistance\" value=\"2.78088\"/>\n";
  file<<"<float name=\"fov\" value=\"55\"/>\n";
  file<<"<string name=\"fovAxis\" value=\"x\"/>\n";
  file<<"<transform name=\"toWorld\">\n";
  file<<"<lookat target=\""<<target_lookat[0]<<", "<<target_lookat[1]<<", "<<target_lookat[2]<<"\" origin=\""<<camOrg[0]<<", "<<camOrg[1]<<", "<<camOrg[2]<<"\" up=\""<<axis[0]<<", "<<axis[1]<<", "<<axis[2]<<"\"/>\n";
  file<<" </transform>\n";
  file<<"<sampler type=\"ldsampler\">\n";
  file<<"<integer name=\"sampleCount\" value=\"16\"/>\n";
  file<<"</sampler>\n";
  file<<"<film type=\"hdrfilm\">\n";
  file<<"<boolean name=\"banner\" value=\"false\"/>\n";
  file<<"<integer name=\"height\" value=\"1080\"/>\n";
  file<<"<string name=\"pixelFormat\" value=\"rgb\"/>\n";
  file<<"<integer name=\"width\" value=\"1940\"/>\n";
  file<<"<rfilter type=\"gaussian\"/>\n";
  file<<"</film>\n";
  file<<"</sensor>\n";
  file<<"<emitter type=\"sunsky\">\n";
  file<<"<transform name=\"toWorld\">\n";
  file<<"<rotate x=\"1\" angle=\"90\"/>\n";
  file<<"</transform>\n";
  file<<"<float name=\"scale\" value=\"10\"/>\n";
  file<<"<float name=\"sunRadiusScale\" value=\"1\"/>\n";
  file<<"<vector name=\"sunDirection\" x=\"-1\" y=\"3\" z=\"2\"/>\n"; //other
  //file<<"<vector name=\"sunDirection\" x=\"5\" y=\"5\" z=\"2\"/>\n"; // lake
  // file<<"<vector name=\"sunDirection\" x=\"5\" y=\"0\" z=\"5\"/>\n"; //coast (boat)
  file<<"</emitter>\n";
  file<<"</scene>\n";
}


void Viewer::script() {

  /**boats **/
  // Vector3 p0 = Vector3(-0.497683, 0.334127, 3.2959);//-0.487125, 0.239223, 2.63672);
  // Vector3 p1 = Vector3(-0.486465, 0.218313, 2.47193);
  // Vector3 p2 = Vector3(0.165719, -0.128639, 0.905468);// -0.310521, -0.10425, 1.39046);
  // if (time == 1) {
  //   Vector3 camOrg =
  //     _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  //   Vector3 dir    = (p0 - camOrg);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // } else if (time*dt_ >= 20 && time*dt_ < 25) {
  //   Vector3 dir = (p1 - p0)/(FLOAT)(5.0/dt_);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // } else if (time*dt_ >= 50 && time*dt_ < 55) {
  //   Vector3 dir = (p2 - p1)/(FLOAT)(5.0/dt_);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // }

  
  /** rain **/
  // Vector3 p0 = Vector3(1.1743, -0.685673, 1.21373);
  // Vector3 p1 = Vector3(1.70951, 0.406636, 0.42906);
  // if (time == 1) {
  //   Vector3 camOrg =
  //     _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  //   Vector3 dir    = (p0 - camOrg);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // } else if (time >= 300 && time < 400) {
  //   Vector3 dir = (p1 - p0)/(FLOAT)(100);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // }


  /** sig **/
  // Vector3 p1 = Vector3(-0.373547, -0.771941, 0.443292);
  // Vector3 p0 = Vector3(-0.373547, -9.7779, 9.44925);
  // if (time == 1) {
  //   Vector3 camOrg =
  //     _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  //   Vector3 dir    = (p1 - camOrg);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // } else  if (time == 10) {
  //   Vector3 camOrg =
  //     _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  //   Vector3 dir    = (p0 - camOrg);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // } else if (time >= 30 && time < 1300) {
  //   Vector3 dir = (p1 - p0)/(FLOAT)(1000);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // }

  /** lake **/
  // if (time == 1) {
  // Vector3 p0 = //Vector3(-0.219811, -0.276124, 0.254966);//
  // Vector3(-0.199574, -0.318405, 0.288495);
  //  Vector3 camOrg =
  //    _cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  //  Vector3 dir    = (p0 - camOrg);
  //  _cameraObject->translate(dir);
    
  //  Vector2i ori = {defaultFramebuffer.viewport().sizeX()/2, defaultFramebuffer.viewport().sizeY()/2};
  //  Vector3 mpos3 = mousePlanePosition(ori);
  //  Vector3 look0 = (mpos3 - camOrg).normalized();
  //  Vector3 look1 = //Vector3(0.376079, 0.740654, -0.556772).normalized();
  //    (Vector3(0.0148947, 0.106238, 0) - camOrg).normalized();
  //  Vector3 up = cross(look0, look1);
  //  FLOAT a = up.length();
  //  up = up.normalized();
  //  _cameraObject->rotate(Rad{a}, up);
  //  up_cam = true;
  // }
  
  /** drops **/
  //  if (time == 1) {
  //    _surface.addDropW(VEC2(16.6742, 23.6137 ), 6, 2);
  //     _surface.addDropW(VEC2(11.0971, 20.5134), 2, 1400);
  //     _surface.addDropW(VEC2(10.4419, 14.1067), 0.8, 2200);
  //     _surface.addDropW(VEC2(9.57121, 10.4347), 0.275, 2950);
  //     _surface.addDropW(VEC2(9.60529, 9.76842), 0.10, 3500);
  //     //     _surface.addDropW(VEC2(9.02707, 9.55576), 0.036, 4200);
  //     _surface.addDropW(VEC2( 9.07655, 9.45899)/*9.0451, 9.49047)*/, 0.036, 4080);
  //   }

  //  Vector3 p_1 = Vector3(-0.0220484, -1.14064, 1.50306),
  //    p0 = Vector3(-0.0220484, -0.840013, 1.20244),
  //    p1 = Vector3(-0.305633, -0.414691, 0.580752),
  //    p2 = Vector3( -0.298155, -0.456909, 0.329779),// -0.305633 -0.269497 0.435554
  //    p3 = Vector3(-0.35117, -0.426461, 0.133325), //-0.295647 -0.314363 0.184124
  //    p4 = Vector3(-0.360225, -0.375317, 0.0298571),
  //    // p5 = Vector3(-0.396073, -0.375344, 0.0185165),
  //    p6 = Vector3(-0.395358, -0.376365, 0.00948717),//-0.397149, -0.372401, 0.0070309),
  //    p7 = Vector3(0.441551, -2.18081, 2.88064);//(-0.908, -1.31755, 1.80383);
   
  //  FLOAT t_1 = 1000, t0 = 1200, t1 = 1350,
  //     t2 = 2000, t3 = 2150,
  //     t4 = 2750, t5 = 2900,
  //     t6 = 3300, t7 = 3450,
  //    t8 = 3900, t9 = 4050,
  //     // t10 = 4300, t11 = 4350,
  //    // t12 = 4550, t13 = 4750;
  //    t12 = 4250, t13 = 4450;
  //  if (time == 1) {
  //    Vector3 camOrg =
  // 	_cameraObject->transformation().transformPoint(Vector3{0.f, 0.f, 0.f});
  //    Vector3 dir    = (p0 - camOrg);
  //    _cameraObject->translate(dir);
  //  // } else if (time >= 20 && time < t_1) {
  //  //  Vector3 dir    = (p0 - p_1)/(FLOAT)(t_1 - 20);
  //  //  _cameraObject->translate(dir);
  //    up_cam = true;
  // } else if (time >= t0 && time < t1) {
  //   Vector3 dir    = (p1 - p0)/(FLOAT)(t1 - t0);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  //    } else if (time >= t2 && time < t3) {
  //   Vector3 dir    = (p2 - p1)/(FLOAT)(t3 - t2);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // } else if (time >= t4 && time < t5) {
  //   Vector3 dir    = (p3 - p2)/(FLOAT)(t5 - t4);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // } else if (time >= t6 && time < t7) {
  //   Vector3 dir    = (p4 - p3)/(FLOAT)(t7 - t6);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // } else if (time >= t8 && time < t9) {
  //   Vector3 dir    = (p6 - p4)/(FLOAT)(t9 - t8);
  //   _cameraObject->translate(dir);
  //   up_cam = true;
  // // } else if (time >= t10 && time < t11) {
  // //   Vector3 dir    = (p6 - p5)/(FLOAT)(t11 - t10);
  // //  _cameraObject->translate(dir);
  // } else if (time >= t12 && time < t13) {
  //   Vector3 dir    = (p7 - p6)/(FLOAT)(t13 - t12);
  //   _cameraObject->translate(dir);
  //      up_cam = true;
  // }

}


	

std::ofstream& Viewer::scriptMitsuba(std::ofstream& file) {
  /* boats */
  // VEC2 pos0(8, 2);
  // FLOAT speed = 2;
  // FLOAT t0 = 0, t1 = 30;
  //  FLOAT maink = gravity_/(4.0*pow(speed*cos(35.3*M_PI/180.0f), 2));
  // FLOAT maxk = gravity_/(4.0*pow(speed, 2));
  // FLOAT mainwl = 2.0*M_PI/maink;
  // FLOAT maxwl = 2.0*M_PI/maxk;

  // int t_start = t0/dt_;
  // int t_stop = t1/dt_;
  
  // VEC2 step2d = speed*dt_*VEC2(0, 1);
  // if (time >= t_start && time < t_stop) {
  //   VEC2 pos_cur = pos0 + (time - t_start)*step2d;
  //   VEC2 pv = world2viewer(pos_cur);
  //  file<<"<shape type=\"heightfield\">\n";
  //  file<<"<string name=\"filename\" value=\"Textures/boat.png\"/>\n";
  //  file<<"<float name=\"scale\" value=\""<<0.2<<"\"/>\n";
  //  file<<"<transform name=\"toWorld\">\n";
  //  file<<"<rotate x=\"0\" y=\"0\" z=\"1\" angle=\"0\"/>\n";
  //  file<<"<scale x=\""<<0.25<<"\" y=\""<<-0.25<<"\" z=\"1\"/>\n";
  //  file<<"<translate x=\""<<pv(0)<<"\" y=\""<<pv(1)-0.2<<"\" z=\""<<-0.03<<"\"/>\n";
  //  //file<<"<translate x=\""<<pv(0)<<"\" y=\""<<pv(1)<<"\" z=\""<<-0.05<<"\"/>\n";
  //  file<<" </transform>\n";
  //  file<<"<bsdf type=\"diffuse\">\n";
  //  file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //  file<<"</bsdf>\n";
  //  file<<"</shape>\n";
  // }

  // pos0(11, 5);
  // speed = 1;
  // t0 = 20; t1 = 80;
  // maink = gravity_/(4.0*pow(speed*cos(35.3*M_PI/180.0f), 2));
  // maxk = gravity_/(4.0*pow(speed, 2));
  // mainwl = 2.0*M_PI/maink;
  // maxwl = 2.0*M_PI/maxk;

  //  t_start = t0/dt_;
  //  t_stop = t1/dt_;
  
  //  step2d = speed*dt_*VEC2(0, 1);
  //  if (time >= t_start && time < t_stop) {
  //    VEC2 pos_cur = pos0 + (time - t_start)*step2d;
  //    VEC2 pv = world2viewer(pos_cur);
  //    file<<"<shape type=\"heightfield\">\n";
  //    file<<"<string name=\"filename\" value=\"Textures/boat.png\"/>\n";
  //    file<<"<float name=\"scale\" value=\""<<0.05<<"\"/>\n";
  //    file<<"<transform name=\"toWorld\">\n";
  //    file<<"<rotate x=\"0\" y=\"0\" z=\"1\" angle=\"0\"/>\n";
  //    file<<"<scale x=\""<<0.05<<"\" y=\""<<-0.05<<"\" z=\"1\"/>\n";
  //    file<<"<translate x=\""<<pv(0)<<"\" y=\""<<pv(1)-0.05<<"\" z=\""<<-0.006<<"\"/>\n";
  //    //file<<"<translate x=\""<<pv(0)<<"\" y=\""<<pv(1)<<"\" z=\""<<-0.05<<"\"/>\n";
  //    file<<" </transform>\n";
  //    file<<"<bsdf type=\"diffuse\">\n";
  //    file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //    file<<"</bsdf>\n";
  //    file<<"</shape>\n";
  //  }

  // pos0(15, 7);
  // speed = 0.5;
  // t0 = 50; t1 = 170;
  // maink = gravity_/(4.0*pow(speed*cos(35.3*M_PI/180.0f), 2));
  // maxk = gravity_/(4.0*pow(speed, 2));
  // mainwl = 2.0*M_PI/maink;
  // maxwl = 2.0*M_PI/maxk;

  //  t_start = t0/dt_;
  //  t_stop = t1/dt_;
  
  //  step2d = speed*dt_*VEC2(0, 1);
  //  if (time >= t_start && time < t_stop) {
  //    VEC2 pos_cur = pos0 + (time - t_start)*step2d;
  //    VEC2 pv = world2viewer(pos_cur);
  //    file<<"<shape type=\"heightfield\">\n";
  //    file<<"<string name=\"filename\" value=\"Textures/boat.png\"/>\n";
  //    file<<"<float name=\"scale\" value=\""<<0.025<<"\"/>\n";
  //    file<<"<transform name=\"toWorld\">\n";
  //    file<<"<rotate x=\"0\" y=\"0\" z=\"1\" angle=\"0\"/>\n";
  //    file<<"<scale x=\""<<0.025<<"\" y=\""<<-0.025<<"\" z=\"1\"/>\n";
  //    file<<"<translate x=\""<<pv(0)<<"\" y=\""<<pv(1)-0.025<<"\" z=\""<<-0.003<<"\"/>\n";
  //    //file<<"<translate x=\""<<pv(0)<<"\" y=\""<<pv(1)<<"\" z=\""<<-0.05<<"\"/>\n";
  //    file<<" </transform>\n";
  //    file<<"<bsdf type=\"diffuse\">\n";
  //    file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //    file<<"</bsdf>\n";
  //    file<<"</shape>\n";
  //  }


  /* dropping bunnies */
  //  int t1 = 20, t2 = 1400, t3 = 2200, t4 = 2950, t5 = 3500, t6 = 4200, t7 = 4400;
  // FLOAT coef1 = 5, coef2 = 2, coef3 = 0.8, coef4 = 0.275, coef5 = 0.1, coef6 = 0.036, coef7 = 0.01;
  // VEC2 p1(16.6742, 23.6137 );
  // VEC2 p2(11.0971, 20.5134);
  // VEC2 p3(10.4419, 14.1067);
  // VEC2 p4(9.57121, 10.4347);
  // VEC2 p5(9.60529, 9.76842);
  // VEC2 p6(9.02707, 9.55576);
  // VEC2 p7(9.03707, 9.46106);
  // int t_drop;// =  _surface.getTDrop(5, 50);
  // int fall_time = 20;
  // FLOAT coef;
  // VEC2 p;
  // int start; //= t_drop - fall_time;
  // // FLOAT coef = 0.1*5;

  // coef = coef1;
  // t_drop = _surface.getTDrop(coef1, t1);
  // start = t_drop - fall_time;
  // p = world2viewer(p1);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  //   coef = coef2;
  // t_drop = _surface.getTDrop(coef2, t2);
  // start = t_drop - fall_time;
  // p = world2viewer(p2);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  // coef = coef3;
  // t_drop = _surface.getTDrop(coef3, t3);
  // start = t_drop - fall_time;
  // p = world2viewer(p3);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  // coef = coef4;
  // t_drop = _surface.getTDrop(coef4, t4);
  // start = t_drop - fall_time;
  // p = world2viewer(p4);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  // coef = coef5;
  // t_drop = _surface.getTDrop(coef5, t5);
  // start = t_drop - fall_time;
  // p = world2viewer(p5);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  // coef = coef4;
  // t_drop = _surface.getTDrop(coef4, t4);
  // start = t_drop - fall_time;
  // p = world2viewer(p4);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  // coef = coef5;
  // t_drop = _surface.getTDrop(coef5, t5);
  // start = t_drop - fall_time;
  // p = world2viewer(p5);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  // coef = coef6;
  // t_drop = _surface.getTDrop(coef6, t6);
  // start = t_drop - fall_time;
  // p = world2viewer(p6);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  //   coef = coef7;
  // t_drop = _surface.getTDrop(coef7, t7);
  // start = t_drop - fall_time;
  // p = world2viewer(p7);
  // if (time > start && time < start + 2*fall_time) {
  //   FLOAT z = 1.0 - ((FLOAT)time-start)/fall_time;
  //   INFO("z  "<<z<<" "<<time);
  //   file<<"<shape type=\"obj\">\n";
  //   file<<"<string name=\"filename\" value=\"bunny.obj\"/>\n";
  //   //  file<<"<float name=\"scale\" value=\""<<coef<<"\"/>\n";
  //   file<<"<transform name=\"toWorld\">\n";
  //   file<<"<rotate x=\"1\" y=\"0\" z=\"0\" angle=\"90\"/>\n";
  //   file<<"<scale x=\""<<0.3*coef<<"\" y=\""<<0.3*coef<<"\" z=\""<<0.3*coef<<"\"/>\n";
  //   file<<"<translate x=\""<<p(0)<<"\" y=\""<<p(1)<<"\" z=\""<<z<<"\"/>\n";
  //   file<<" </transform>\n";
  //   file<<"<bsdf type=\"diffuse\">\n";
  //   file<<"<srgb name=\"reflectance\" value=\"#6d7185\"/>\n";
  //   file<<"</bsdf>\n";
  //   file<<"</shape>\n";
  // }

  return file;
}
