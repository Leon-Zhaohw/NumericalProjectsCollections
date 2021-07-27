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

#include "setting.h"

#include <iostream>
#define GLUT_DISABLE_ATEXIT_HACK

#ifdef APPLE
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
//#include <glog/logging.h>
//#include <gflags/gflags.h>
#include <sys/time.h>
#include <string>
#include <vector>

#include "2D/laplacian_fluid_2D.h"
#include "util/timer.h"
#include "2D/obstacle_2d.h"
#include "util/QUICKTIME_MOVIE.h"
#include "util/get_current_time.h"
#include "util/SIMPLE_PARSER.h"

LaplacianFluid2D* fluid = NULL;

static int fluid_window_x = 1024, spectra_window_x = 1700; 
static int fluid_window_y = 1024, spectra_window_y = 512;
static int fluid_win_id;
static int spectra_win_id;
static int mouse_down[3];
static int omx, omy, mx, my;
const double spectra_draw_factor = 2.0;
double particle_length_factor = 0.1;
static bool addedDensity = false;
static bool dvel = false;
/*
DEFINE_int32(xRes, 512, "The xReslution of the fluid.");
DEFINE_int32(yRes, 512, "The yReslution of the fluid.");
DEFINE_int32(basis_dim_root, 20, "The root square of the basis.");
DEFINE_double(dt, 0.1, "The time step of the simulation.");
DEFINE_double(buoyancy, 0.0001, "The buoyancy of the fluid.");
DEFINE_double(viscosity, 0.005, "The viscosity of the fluid.");

DEFINE_bool(add_density, true, "Whether to add density into the fluid.");
DEFINE_double(added_smoke_density, 2.0, "The density of added smoke.");
DEFINE_int32(source_smoke_xpos, 0, "The x coordinate of the added smoke source.");
DEFINE_int32(source_smoke_ypos, 0, "The y coordinate of the added smoke source.");
DEFINE_int32(source_smoke_width, 0, "The width of the added smoke source.");
DEFINE_int32(source_smoke_height, 0, "The height of the added smoke source.");
DEFINE_bool(add_density_once, false, "Only add the density to the first of the frame");
DEFINE_bool(add_density_mirror, false, "whether to add a mirror to current density");

DEFINE_string(basis_type, "all_dirichlet", "The basis type. Can be all_dirichlet,"
              "three_dirichlet_one_neumann, two_neumann_x");
DEFINE_string(integrator_type, "trapezoidal", "The integrator_type. Can be RK4, semi_implicit"
              "BiCGSTAB");
DEFINE_int32(total_frame, 2000, "The total frame want to simulate.");
DEFINE_int32(total_num_particles, 1000, "The total number of particles");
DEFINE_double(particle_length_factor, 0.1, "How long is the particles.");
DEFINE_double(force_scale, 0.5, "The force scale to apply the force by mouse");
DEFINE_string(Tensor_file, "", "The file which stores all the tensor.");
// Obstacles.
DEFINE_bool(use_obstacles, false, "Whether to use obstacles in simulation.");
DEFINE_string(obstacle_img_file, "", "The file that contains the image of obstacle.");
DEFINE_double(obstacle_force_scale, 20, "The scale of the penalty force of the obstacle.");
DEFINE_bool(obstacle_implicit_method, true, "Whether to use implicit method to handle obstacle.");

DEFINE_bool(add_bouyancy_once, false, "Only add the bouyancy for the first iteration");
DEFINE_bool(capture_video, false, "Whether to capture the video.");
DEFINE_string(preview_fpath, "", "The folder to put the preview video.");
*/

QUICKTIME_MOVIE movie;

int FLAGS_xRes = 512;
int FLAGS_yRes = 512;
int FLAGS_basis_dim_root = 20;
double FLAGS_dt = 0.1;
double FLAGS_buoyancy = 0.0001;
double FLAGS_viscosity = 0.005;
  
bool FLAGS_add_density = true;
double FLAGS_added_smoke_density = 2.0;
int FLAGS_source_smoke_xpos = 0;
int FLAGS_source_smoke_ypos = 0;
int FLAGS_source_smoke_width = 0;
int FLAGS_source_smoke_height = 0;
bool FLAGS_add_density_once = false;
bool FLAGS_add_density_mirror = false;

std::string FLAGS_basis_type = "all_dirichlet";
std::string FLAGS_integrator_type = "trapezoidal";
int FLAGS_total_frame = 2000;
int FLAGS_total_num_particles = 1000;
double FLAGS_particle_length_factor = 0.1;
double FLAGS_force_scale = 0.5;
std::string FLAGS_Tensor_file = "";

bool FLAGS_use_obstacles = false;
std::string FLAGS_obstacle_img_file = "";
double FLAGS_obstacle_force_scale = 20;
bool FLAGS_obstacle_implicit_method = true;
bool FLAGS_add_bouyancy_once = false;
bool FLAGS_capture_video = false;
std::string FLAGS_preview_fpath = "";

namespace {
void InitializeObstacleParams(ObstacleParams* params_) {
  if (!FLAGS_use_obstacles) {
    std::cout <<  "Do not use obstacles." << std::endl;
    params_->use_obstacles = false;
  } else {
    params_->img_file_name = FLAGS_obstacle_img_file;
    params_->use_obstacles = true;
    params_->obstacle_force_scale = FLAGS_obstacle_force_scale;
    params_->handle_obstacle_implicit = FLAGS_obstacle_implicit_method;
  }
}

void SavePreviewAndQuit() {
  if (FLAGS_capture_video) {
    std::string datetime = getCurrentDateTime(true);
    std::stringstream outname;
    outname << FLAGS_preview_fpath << datetime << ".mov";
    movie.writeMovie(outname.str().c_str());
    // reset the movie object
    movie = QUICKTIME_MOVIE();
  }
  delete fluid;
  
  std::cout <<  "Exitting..." << std::endl;
  exit(0);
}

}  // namespace
static void fluid_pre_display ( void ) {
  glViewport ( 0, 0, fluid_window_x, fluid_window_y );
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ();
  gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
  glClearColor ( 1.0f, 1.0f, 1.0f, 1.0f );
  glClear ( GL_COLOR_BUFFER_BIT );
}

static void spectra_pre_display(void) {
  glViewport ( 0, 0, spectra_window_x, spectra_window_y);
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ();
  gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
  glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
  glClear ( GL_COLOR_BUFFER_BIT );
}

static void key_func ( unsigned char key, int x, int y ) {
  switch ( key ) {
    case 'c':
    case 'C':
      //clear_data ();
      break;
    case 'q':
    case 'Q':
      //free_data ();
      fluid->Quit();
      exit (0);
      break;
    case 'v':
    case 'V':
      dvel = !dvel;
      break;
    case 'r':
    case 'R':
      fluid->ReSeedParticles();
      break;
  }
}

static void SpecialInput(int key, int x, int y) {
  switch(key) {
    case GLUT_KEY_UP:
      fluid->AddVelocityToObstacles(0,0.1);
      break;
    case GLUT_KEY_DOWN:
      fluid->AddVelocityToObstacles(0,-0.1);
      break;
    case GLUT_KEY_LEFT:
      fluid->AddVelocityToObstacles(-0.1,0);
      break;
    case GLUT_KEY_RIGHT:
      fluid->AddVelocityToObstacles(0.1,0);
      break;
  }
}

static void mouse_func ( int button, int state, int x, int y ) {
  omx = mx = x;
  omx = my = y;

  mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y ) {
  mx = x;
  my = y;
}

static void reshape_func_fluid ( int width, int height ) {
  glutSetWindow ( fluid_win_id );
  glutReshapeWindow ( width, height );
  fluid_window_x = width;
  fluid_window_y = height;
}

static void reshape_func_extra(int width, int height) {
  glutSetWindow ( spectra_win_id );
  glutReshapeWindow ( width, height );
  spectra_window_x = width;
  spectra_window_y = height;
}

void get_from_UI () {
  int i, j;
  bool addDens = FLAGS_add_density;
  if (FLAGS_add_density_once) {
    addDens = addDens && !addedDensity;
  }
  
  if (addDens) {
    fluid->AddSmokeTestCase(FLAGS_source_smoke_xpos, FLAGS_source_smoke_ypos,
                            FLAGS_source_smoke_width, FLAGS_source_smoke_height);
    if (FLAGS_add_density_once) {
      fluid->AddForceImpluse(FLAGS_source_smoke_xpos, FLAGS_source_smoke_ypos,
                            FLAGS_source_smoke_width+00, FLAGS_source_smoke_height, 0.5);
    }
    addedDensity = true;
  }
  if ( !mouse_down[0] && !mouse_down[2] ) {return;}
  
  i = (int)((       mx /(float)fluid_window_x)*(FLAGS_xRes - 1));
  j = (int)(((fluid_window_y-my)/(float)fluid_window_y)*(FLAGS_yRes - 1));
  int index = i + j*FLAGS_xRes;
  if ( i<1 || i>(FLAGS_xRes - 1) || j<1 || j>(FLAGS_yRes - 1) ) return;
  if (mouse_down[0]) {
    fluid->AddForce(i,j,FLAGS_force_scale *(mx - omx),
                    FLAGS_force_scale * (omy - my));
  }
  if ( mouse_down[2] ) {
    fluid->AddSmoke(i,j, 5);
  }
  omx = mx;
  omy = my;
  return;
}

static void idle_func ( void ) {
  get_from_UI();
  {  fluid->Step();
    //NunOfFrames++;
  }
//  calculateFPS();
     // calculateFPS();
  
  if (fluid->quit_) {
    SavePreviewAndQuit();
  }
  
  glutSetWindow ( fluid_win_id );
  glutPostRedisplay ();
  glutSetWindow(spectra_win_id);
  glutPostRedisplay();
}

static void post_display ( void ) {
  glutSwapBuffers ();
}

static void fluid_display_func ( void ) {
  fluid_pre_display ();
  if ( dvel ) fluid->DrawVelocity();
  else  fluid->DrawDensity();
  //fluid->DrawVort();
  fluid->DrawParticles(particle_length_factor);
  fluid->DrawObstacles();
  // The visualization of the velocity subject the aliasing of the screen.
  
  if (FLAGS_capture_video) {
      movie.addFrameGL(); // This function will eat lots of memory afterwards.
  }
  
  post_display ();
  
}

static void spectra_display_func(void) {
  spectra_pre_display();
  // Plot the spectra coefficients.
  fluid->DrawCoefficients(spectra_draw_factor);
  // if (FLAGS_capture_video) {
  //    movie.addFrameGL(); // This function will eat lots of memory afterwards.
  // }
  post_display();
}

static void open_glut_window ( void ) {
  glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

  glutInitWindowPosition ( 0, 0 );
  glutInitWindowSize ( fluid_window_x, fluid_window_y );
  fluid_win_id = glutCreateWindow ( "2d" );

  glClearColor ( 1.0f, 1.0f, 1.0f, 1.0f );
  glClear ( GL_COLOR_BUFFER_BIT );
  glutSwapBuffers ();
  glClear ( GL_COLOR_BUFFER_BIT );
  glutSwapBuffers ();
  fluid_pre_display ();
  glutKeyboardFunc ( key_func );
  glutSpecialFunc(SpecialInput);
  glutMouseFunc ( mouse_func );
  glutMotionFunc ( motion_func );
  glutReshapeFunc ( reshape_func_fluid );
  glutIdleFunc ( idle_func );
  glutDisplayFunc ( fluid_display_func );
  
  glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
  glutInitWindowPosition ( 580, 0 );
  glutInitWindowSize ( spectra_window_x, spectra_window_x );
  spectra_win_id = glutCreateWindow("energy");
  glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
  glClear ( GL_COLOR_BUFFER_BIT );
  glutSwapBuffers ();
  glClear ( GL_COLOR_BUFFER_BIT );
  glutSwapBuffers ();
  spectra_pre_display();
  glutReshapeFunc (reshape_func_extra);
  glutDisplayFunc(spectra_display_func);
}
/*

*/

int main(int argc, char ** argv) {
  //google::ParseCommandLineFlags(&argc, &argv, true);
  //google::InitGoogleLogging(argv[0]);
  
  if (argc != 2)
  {
    std::cout << " Usage: " << argv[0] << " *.cfg" << std::endl;
    return 0;
  }
  SIMPLE_PARSER parser(argv[1]);
  FLAGS_xRes = parser.getInt("xRes", 512);
  FLAGS_yRes = parser.getInt("yRes", 512);
  FLAGS_basis_dim_root = parser.getInt("basis_dim_root", 20);
  FLAGS_dt = parser.getFloat("dt", 0.1);
  FLAGS_buoyancy = parser.getFloat("buoyancy", 0.0001);
  FLAGS_viscosity = parser.getFloat("viscosity", 0.005);
  
  FLAGS_add_density = parser.getBool("add_density", true);
  FLAGS_added_smoke_density = parser.getFloat("added_smoke_density", 2.0);
  FLAGS_source_smoke_xpos = parser.getInt("source_smoke_xpos", 0);
  FLAGS_source_smoke_ypos = parser.getInt("source_smoke_ypos", 0);
  FLAGS_source_smoke_width = parser.getInt("source_smoke_width", 0);
  FLAGS_source_smoke_height = parser.getInt("source_smoke_height", 0);
  FLAGS_add_density_once = parser.getBool("add_density_once", false);
  FLAGS_add_density_mirror = parser.getBool("add_density_mirror", false);
  // std::cout << FLAGS_add_density << "\n" << FLAGS_add_density_once << std::endl;
  
  FLAGS_basis_type = parser.getString("basis_type", "all_dirichlet");
  FLAGS_integrator_type = parser.getString("integrator_type", "trapezoidal");
  FLAGS_total_frame = parser.getInt("total_frame", 2000);
  FLAGS_total_num_particles = parser.getInt("total_num_particles", 1000);
  FLAGS_particle_length_factor = parser.getFloat("particle_length_factor", 0.1);
  FLAGS_force_scale = parser.getFloat("force_scale", 0.5);
  FLAGS_Tensor_file = parser.getString("Tensor_file", "");
  FLAGS_use_obstacles = parser.getBool("use_obstacles", false);
  FLAGS_obstacle_img_file = parser.getString("obstacle_img_file", "");
  FLAGS_obstacle_force_scale = parser.getFloat("obstacle_force_scale", 20.0);
  FLAGS_obstacle_implicit_method = parser.getBool("obstacle_implicit_method", true);
  FLAGS_add_bouyancy_once = parser.getBool("add_bouyancy_once", false);
  FLAGS_capture_video = parser.getBool("capture_video", false);
  FLAGS_preview_fpath = parser.getString("preview_fpath", "");
  
  glutInit (&argc, argv);
  
  fluid = new LaplacianFluid2D(FLAGS_xRes, FLAGS_yRes, FLAGS_basis_dim_root,
      FLAGS_dt, FLAGS_basis_type, FLAGS_integrator_type, FLAGS_buoyancy, FLAGS_viscosity,
      FLAGS_total_num_particles, FLAGS_total_frame, FLAGS_added_smoke_density, FLAGS_Tensor_file, FLAGS_add_density_mirror);
  particle_length_factor = FLAGS_particle_length_factor;
  ObstacleParams obstacleparams;
  InitializeObstacleParams(&obstacleparams);
  fluid->InitializeObstacles(obstacleparams);
  // Whether to initialize velocity from Dedalus simulation.
  // if (FLAGS_dedalus_velo_folder.size() != 0) {
  //  std::string UInitfile = FLAGS_dedalus_velo_folder + "UInit";
  //  std::string VInitfile = FLAGS_dedalus_velo_folder + "VInit";
  //  fluid->InitVelocityFromDedalus(UInitfile, VInitfile, FLAGS_dedalus_xRes, FLAGS_dedalus_yRes);
  //}
  fluid->addBounyancyOnce = FLAGS_add_bouyancy_once;
  
  open_glut_window ();
  glutMainLoop ();

  return 0;
}
