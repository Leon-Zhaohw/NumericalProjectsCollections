/*
 This file is part of SSFR (Zephyr).
 
 Zephyr is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 Zephyr is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with Zephyr.  If not, see <http://www.gnu.org/licenses/>.
 
 Copyright 2018 Qiaodong Cui (qiaodong@ucsb.edu)
 */

// Viewer for 3D laplacian fluid simulation.
#include "setting.h"
#include "util/QUICKTIME_MOVIE.h"

#include "Eigen"
#include <fftw3.h>
#include <iostream>
#include <sstream>
#define GLUT_DISABLE_ATEXIT_HACK

#include "glvu.h"
#ifdef APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
//#include <glog/logging.h>
//#include <gflags/gflags.h>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <vector>

#include "util/SIMPLE_PARSER.h"
#include "3D/fluid_3D_DCT.h"
#include "3D/obstacle_box_3d.h"
#include "util/glPrintString.h"
#include "util/get_current_time.h"
#include "util/trackball.h"
#include "util/timer.h"

FLUID_3D_DCT* fluid = NULL;

int light;
// Toggle light
int lp;

static int fluid_win_id;
static int spectra_win_id;

GLfloat zt =-1.5f; // depth into the screen.
GLfloat xt = 0, yt = 0;

Eigen::Vector3f lpos(0,5,0);
GLfloat LightPosition[] = { lpos(0), lpos(1), lpos(2), 1.0f };

int window_x_res, window_y_res;
static int spectra_window_x = 512;
static int spectra_window_y = 512;

int mouse_ax, mouse_ay;
float mouse_x_norm, mouse_y_norm;

float _quat[4];  // view rotation quaternion
float boxXLength, boxYLength, boxZLength;
float aplha_multipler = 0.02;

int low_cut = 5;

bool add_density = false;
bool use_disk_smoke = false;
bool write_dens_pbrt = false;

// animate it this frame?
bool animate = true;
void runEverytime();

// Quicktime movie to capture to
QUICKTIME_MOVIE movie;
bool capture_video = true;
std::string preview_path("./preview/");
float total_step_time = 0;
int simulationSnapshots = 20;

GLVU glvu;

namespace {

void InitializeObstacleParams(ObstacleParams3D* params_, SIMPLE_PARSER& parser) {
  
  bool use_obstacles = parser.getBool("use_obstacles", false);
  bool move_obstacle = parser.getBool("move_obstacle", false);
  std::string obstacle_list_file = parser.getString("obstacle_list_file", "");
  
  if (!use_obstacles) {
    std::cout <<  "Do not use obstacles." << std::endl;
    params_->use_obstacles = false;
  } else {
   // params_->img_file_name = FLAGS_obstacle_img_file;
    params_->use_obstacles = true;
    params_->obstacle_force_scale = 1.0;
    params_->handle_obstacle_implicit = false;
    params_->move_obstacle = move_obstacle;
    params_->obstacle_type = "";
    params_->obstacle_file = "";
    params_->obstacle_list = obstacle_list_file;
  }
}
///////////////////////////////////////////////////////////////////////
// GL and GLUT callbacks
///////////////////////////////////////////////////////////////////////
void glutDisplay()
{
  glvu.BeginFrame();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // glEnable(GL_DEPTH_TEST);
    glRotatef(180.0, 0.0, 1.0, 0.0);
    glTranslatef(-boxXLength,-boxYLength,-boxZLength);
    // 
    DrawCube(boxXLength, boxYLength, boxZLength);
    glPushMatrix();
    fluid->DrawObstacles();
    glPopMatrix();
  
    glPushMatrix();
    glTranslatef(boxXLength,boxYLength,boxZLength);
    fluid->DrawSmoke();
    glPopMatrix();
    fluid->DrawParticles(0.01);
    /*glPushMatrix();
      glTranslatef(0.5, 0.5, 0.5);
      // fluid->density().draw();
      fluid->DrawSmoke();
      fluid->density().drawBoundingBox();
    glPopMatrix();
    
    glPushMatrix();
    fluid->DrawObstacles();
    glPopMatrix();
    */
    
    //drawAxes();
  glvu.EndFrame();
}

///////////////////////////////////////////////////////////////////////
// animate and display new result
///////////////////////////////////////////////////////////////////////
void glutIdle()
{
  runEverytime();
  if (capture_video) {
      movie.addFrameGL(); // This function will eat lots of memory afterwards.
  }
  if (write_dens_pbrt) {
      fluid->OutputSmokeToFolder();
  }
  glutPostRedisplay();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutKeyboard(unsigned char key, int x, int y)
{
  switch (key)
  {
    case 'p':
      animate = !animate;
      break;
    case 'q':
      exit(0);
      break;
    default:
      break;
  }
  glvu.Keyboard(key,x,y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutSpecial(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:
      break;
    case GLUT_KEY_RIGHT:
      break;
    case GLUT_KEY_UP:
      break;
    case GLUT_KEY_DOWN:
      break;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseClick(int button, int state, int x, int y)
{
  glvu.Mouse(button,state,x,y);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void glutMouseMotion(int x, int y)
{
  glvu.Motion(x,y);
}

//////////////////////////////////////////////////////////////////////////////
// open the GLVU window
//////////////////////////////////////////////////////////////////////////////
int glvuWindow()
{
  char title[] = "3D Viewer";
  glvu.Init(title,
            GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH,
            0, 0, 800, 800);

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glShadeModel(GL_SMOOTH);
    
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
    
  gluPerspective(45.0f,(GLfloat)800/(GLfloat)800,0.1f,100.0f);
    
  glMatrixMode(GL_MODELVIEW);

  glLightfv(GL_LIGHT1, GL_POSITION,LightPosition);
  glEnable(GL_LIGHT1);      
  
  glutDisplayFunc(&glutDisplay);
  glutIdleFunc(&glutIdle);
  glutKeyboardFunc(&glutKeyboard);
  glutSpecialFunc(&glutSpecial);
  glutMouseFunc(&glutMouseClick);
  glutMotionFunc(&glutMouseMotion);

  ///glEnable (GL_BLEND);
  //glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
        Eye(0,0,-1.5), LookAtCntr(0,0,1),  Up(0,1,0);

  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov, Aspect, Near, Far);

  //glvuVec3f center(0.25, 0.25, 0.25);
  glvuVec3f center(0., 0., 0.);
  glvu.SetWorldCenter(center);

  glutMainLoop();

  // Control flow will never reach here
  return EXIT_SUCCESS;
}

void SavePreviewAndQuit() {
  if (capture_video) {
    std::string datetime = getCurrentDateTime(true);
    std::stringstream outname;
    outname << preview_path << datetime << ".mov";
    movie.writeMovie(outname.str().c_str());
    // reset the movie object
    movie = QUICKTIME_MOVIE();
  }
  std::cout <<  "Average DCT time: " << fluid->totalDCTtime_ / static_cast<float>(simulationSnapshots) << std::endl;
  
  std::cout <<  "Total stepping time: " << total_step_time << std::endl;
  std::cout <<  "Exitting..." << std::endl;
  delete fluid;
  exit(0);
}

}  // namespace

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
void runEverytime()
{

  if (animate)
  {
    static int steps = 0;

    // step the sim
    cout << " Simulation step " << steps << endl;
    
    if (add_density) {
      if ( use_disk_smoke) {
        fluid->addSmokeColumn();
      } else {
        fluid->AddSmokeTestCase();
      } 
    }
    
    Timer timer;
    timer.Reset();
    
    fluid->step();
    total_step_time += timer.ElapsedTimeInSeconds();
    std::cout <<  "Time for step: " << timer.ElapsedTimeInSeconds() << std::endl;
  
   
    // check if we're done
    if (steps == simulationSnapshots)
    {
      SavePreviewAndQuit();
     //  exit(0);
    }

    steps++;
  }
}


int main(int argc, char ** argv) {
  //google::ParseCommandLineFlags(&argc, &argv, true);
  //google::InitGoogleLogging(argv[0]);
 
  SIMPLE_PARSER parser(argv[1]);
    
  const int fftw_num_threads = parser.getInt("FFTW_threads", 4);  
  // For problems of smaller size, use less threads are better.
  fftw_init_threads();
  fftw_plan_with_nthreads(fftw_num_threads);
  
  int amplify = 4;
  int xRes = parser.getInt("xRes", 48);
  int yRes = parser.getInt("yRes", 64);
  int zRes = parser.getInt("zRes", 48);
  
  const int max_lenght = std::max(std::max(xRes, yRes), zRes);
  boxXLength = static_cast<float>(xRes) / max_lenght * 0.5;
  boxYLength = static_cast<float>(yRes) / max_lenght * 0.5;
  boxZLength = static_cast<float>(zRes) / max_lenght * 0.5;
  
  Real vorticity = parser.getFloat("vorticity", 0);

  cout << " Using resoluion: " << xRes << " " << yRes << " " << zRes << endl;
  cout << " Using vorticity: " << vorticity << endl;

  unsigned int boundaries[6];
  boundaries[0] = parser.getInt("front", 1);
  boundaries[1] = parser.getInt("back", 1);
  boundaries[2] = parser.getInt("left", 1);
  boundaries[3] = parser.getInt("right", 1);
  boundaries[4] = parser.getInt("top", 0);
  boundaries[5] = parser.getInt("bottom", 0);

  string names[] = {"front", "back", "left", "right", "top", "bottom"};
  for (int x = 0; x < 6; x++)
  {
    cout << " Boundary on " << names[x].c_str() << "\tis set to " << flush;
    if (boundaries[x] == 0)
      cout << "Neumann " << endl;
    else
      cout << "Dirichlet " << endl;
  }
  simulationSnapshots = parser.getInt("simulation snapshots", 20);
  
  bool use_two_phase_smoke = parser.getBool("use_two_phase_smoke", false);
  float added_smoke_density = parser.getFloat("added_smoke_density", 1.0);
  add_density = parser.getBool("add_density", false);
  use_disk_smoke = parser.getBool("use_disk_smoke", false);
  std::string buoyancy_direction = parser.getString("buoyancy_direction", "x");
  Real buoyancy = parser.getFloat("buoyancy", 0.2);
  
  const double density_attenuate_factor =  parser.getFloat("density_attenuate_factor",
                                                           0.995);
  const bool attenuate_density = parser.getBool("attenuate_density", false);
  const std::string density_folder = parser.getString("density_folder","");
  write_dens_pbrt = parser.getBool("write_density_to_PBRT", false);
  
  const Real dt = parser.getFloat("dt", 0.1);
  int buoyancy_step = parser.getInt("buoyancy_step", simulationSnapshots);
  // Initialize components
  fluid = new FLUID_3D_DCT(xRes, yRes, zRes, amplify, use_two_phase_smoke,
                           added_smoke_density, buoyancy_direction , buoyancy, density_attenuate_factor,
                           attenuate_density, density_folder, dt,buoyancy_step,
                           &boundaries[0]);
  
  fluid->vorticityEps() = vorticity;

  //
  VEC3F source_pos; VEC3F source_length;
  source_pos[0] = parser.getFloat("source_xpos", 0.1);
  source_pos[1] = parser.getFloat("source_ypos", 0.1);
  source_pos[2] = parser.getFloat("source_zpos", 0.1);
  
  source_length[0] = parser.getFloat("source_length", 0.1);
  source_length[1] = parser.getFloat("source_width", 0.1);
  source_length[2] = parser.getFloat("source_height", 0.1);
  
  std::string source_smoke_file = parser.getString("source_smoke_file", "");
  
  ObstacleParams3D obsparam;
  InitializeObstacleParams(&obsparam, parser);
  fluid->InitializeObstacles(obsparam);
  fluid->InitializeSourceSmoke(source_pos, source_length, source_smoke_file);
  
  glutInit(&argc, argv);
  glvuWindow();

  return 0;
}
