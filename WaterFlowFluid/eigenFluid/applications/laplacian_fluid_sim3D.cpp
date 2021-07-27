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

#include "3D/laplacian_fluid_3d.h"
#include "3D/obstacle_box_3d.h"
#include "util/trackball.h"
#include "util/glPrintString.h"
#include "util/get_current_time.h"
#include "util/SIMPLE_PARSER.h"

const double spectra_draw_factor = 3.0;
/*
DEFINE_int32(window_x, 1024, "The x resolution of the window");
DEFINE_int32(window_y, 1024, "The y resolution of the window");
DEFINE_int32(xRes, 128, "The x resolution of fluid.");
DEFINE_int32(yRes, 128, "The y resolution of fluid.");
DEFINE_int32(zRes, 256, "The z resolution of fluid.");
DEFINE_int32(basis_dim, 512, "The total number of the basis.");

DEFINE_double(dt, 0.1, "The time step of the simulation.");
DEFINE_double(buoyancy, 0.0001, "The buoyancy of the fluid.");
DEFINE_string(buoyancy_direction, "x", "The direction of bouyancy, can be x, y, z");
DEFINE_double(viscosity, 0.005, "The viscosity of the fluid.");
DEFINE_double(added_smoke_density, 2.0, "The density of added smoke.");

DEFINE_int32(total_num_particles, 1000, "The total number of particles");
DEFINE_double(particle_length_factor, 0.1, "How long is the particles.");
DEFINE_int32(total_frame, 2000, "The total frame want to simulate.");

DEFINE_string(basis_type, "all_dirichlet", "The basis type. Can be all_dirichlet,"
              "one_neumann, two_neumann_x");

DEFINE_string(constant_init_strategy, "principle_x", "Choose which direction as principle"
              " propagation direction. Can be principle_x, principle_y, principle_z, random.");

DEFINE_string(basis_tensor_fname, "", "The file contains the basis defination and tensor.");
DEFINE_string(density_folder, "",  "The folder to store all the density");

DEFINE_double(source_xpos, 0, "x position of smoke.");
DEFINE_double(source_ypos, 0, "y position of smoke.");
DEFINE_double(source_zpos, 0, "z position of smoke.");
DEFINE_double(source_length, 0, "Length of source smoke block.");
DEFINE_double(source_width, 0, "width of source smoke block.");
DEFINE_double(source_height, 0, "Height of source smoke block.");
DEFINE_bool(attenuate_density, false, "Whether to attenuate_smoke over time");
DEFINE_bool(add_density, true, "Whether to add smoke into the domain.");
DEFINE_double(density_attenuate_factor, 1.0, "The const for attenuate_density.");
DEFINE_bool(write_density_to_PBRT, true, "Whether to write density to PBRT.");
DEFINE_bool(use_MacCormack, false, "Whether to use MacCormack advection.");
DEFINE_bool(use_disk_smoke, false, "Whether to use hard coded disk smoke.");
DEFINE_bool(use_sphere_smoke, false, "Whether to use sphere smoke, center is source_pos, and radius is"
            " source_length");
DEFINE_bool(use_two_phase_smoke, false, "Whether to use two phase smoke");
DEFINE_string(source_smoke_file, "", "The file contains source smoke pos, etc.");

// Obstacles.
DEFINE_bool(use_obstacles, false, "Whether to use obstacles in simulation.");
DEFINE_bool(move_obstacle, false, "Whether to move the obstacle around.");
DEFINE_string(obstacle_type, "box", "Which kind of obstacle want, can be box, cylinder, sphere");
DEFINE_string(obstacle_file, "", "The file to read obstacle from.");
DEFINE_string(obstacle_list_file, "", "The file contains list of obstacles.");

DEFINE_string(coefficients_file, "", "The file to write out all the coefficients.");
DEFINE_string(coefficients_file_in, "", "The file to read in coefficient files.");

DEFINE_double(obstacle_force_scale, 20, "The scale of the penalty force of the obstacle.");
DEFINE_bool(obstacle_implicit_method, true, "Whether to use implicit method to handle obstacle.");
DEFINE_bool(capture_video, false, "Whether to capture the preview video.");
DEFINE_string(preview_fpath, "", "The folder to put the preview video.");
DEFINE_int32(FFTW_threads,  1, "Number of threads to use for FFTW.");

DEFINE_int32(buoyancy_step, -1, "The step of applying buoyancy forces.");

DEFINE_bool(interactive, false, "Whether is interactive");
// Solver type.
DEFINE_string(solver_type, "BiCGSTAB", "Which type of solver to use, symmetric_cg, BiCGSTAB");
DEFINE_bool(add_force_with_source, false, "Whether add force at the position of the source smoke.");
*/

int FLAGS_window_x = 1024;
int FLAGS_window_y = 1024;
int FLAGS_xRes = 128;
int FLAGS_yRes = 128;
int FLAGS_zRes = 256;
int FLAGS_basis_dim = 512;

double FLAGS_dt = 0.1;
double FLAGS_buoyancy = 0.0001;
std::string  FLAGS_buoyancy_direction = "x";
double FLAGS_viscosity = 0.005;
double FLAGS_added_smoke_density = 2.0;

int FLAGS_total_num_particles = 1000;
double FLAGS_particle_length_factor = 0.1;
int FLAGS_total_frame = 2000;

std::string FLAGS_basis_type = "all_dirichlet";
std::string FLAGS_constant_init_strategy = "principle_x";
std::string FLAGS_basis_tensor_fname = "";
std::string FLAGS_density_folder = "";

double FLAGS_source_xpos = 0;
double FLAGS_source_ypos = 0;
double FLAGS_source_zpos = 0;
double FLAGS_source_length = 0;
double FLAGS_source_width = 0;
double FLAGS_source_height = 0;
bool FLAGS_attenuate_density = false;
bool FLAGS_add_density = true;
double FLAGS_density_attenuate_factor = 1.0;
bool FLAGS_write_density_to_PBRT = true;
bool FLAGS_use_MacCormack = false;
bool FLAGS_use_disk_smoke = false;
bool FLAGS_use_sphere_smoke = false;
bool FLAGS_use_two_phase_smoke = false;
std::string FLAGS_source_smoke_file = "";

// Obstacles.
bool FLAGS_use_obstacles = false;
bool FLAGS_move_obstacle = false;
std::string FLAGS_obstacle_type = "box";
std::string FLAGS_obstacle_file = "";
std::string FLAGS_obstacle_list_file = "";
std::string FLAGS_coefficients_file = "";
std::string FLAGS_coefficients_file_in = "";

double FLAGS_obstacle_force_scale = 20;
bool FLAGS_obstacle_implicit_method = true;
bool FLAGS_capture_video = false;
std::string FLAGS_preview_fpath = "";
int FLAGS_FFTW_threads = 1;

int FLAGS_buoyancy_step = -1;
bool FLAGS_interactive = false;
// Solver type.
std::string FLAGS_solver_type = "BiCGSTAB";
bool FLAGS_add_force_with_source = false;
double FLAGS_basisweightMultiplierC = 0.0;

/* ascii codes for various special keys */
#define ESCAPE 27
#define PAGE_UP 73
#define PAGE_DOWN 81
#define UP_ARROW 72
#define DOWN_ARROW 80
#define LEFT_ARROW 75
#define RIGHT_ARROW 77

LaplacianFluid3D* fluid = NULL;
// Quicktime movie to capture to
QUICKTIME_MOVIE movie;
QUICKTIME_MOVIE spectra1Movie;
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

static int spectra_window_x = 800;
static int spectra_window_y = 600;

int mouse_ax, mouse_ay;
float mouse_x_norm, mouse_y_norm;

float _quat[4];  // view rotation quaternion
float boxXLength, boxYLength, boxZLength;
float aplha_multipler = 0.05;

int low_cut = 0;

std::string printstr;
bool is_pause = false;
bool add_smoke = true;
bool write_dens_pbrt = false;
bool draw_particles = true;

namespace {
void InitializeObstacleParams(ObstacleParams3D* params_) {
  if (!FLAGS_use_obstacles) {
    std::cout <<  "Do not use obstacles." << std::endl;
    params_->use_obstacles = false;
  } else {
   // params_->img_file_name = FLAGS_obstacle_img_file;
    params_->use_obstacles = true;
    params_->obstacle_force_scale = FLAGS_obstacle_force_scale;
    params_->handle_obstacle_implicit = FLAGS_obstacle_implicit_method;
    params_->move_obstacle = FLAGS_move_obstacle;
    params_->obstacle_type = FLAGS_obstacle_type;
    params_->obstacle_file = FLAGS_obstacle_file;
    params_->obstacle_list = FLAGS_obstacle_list_file;
  }
}

void SavePreviewAndQuit() {
  if (FLAGS_capture_video) {
    std::string datetime = getCurrentDateTime(true);
    std::stringstream outname;
    outname << FLAGS_preview_fpath << datetime << ".mov";
    std::stringstream spename;
    spename << FLAGS_preview_fpath << datetime << "_spe.mov";
    movie.writeMovie(outname.str().c_str());
    spectra1Movie.writeMovie(spename.str().c_str());
    
    // reset the movie object
    movie = QUICKTIME_MOVIE();
    spectra1Movie = QUICKTIME_MOVIE();
  }
  delete fluid;
  
  std::cout <<  "Exitting..." << std::endl;
  exit(0);
}

}  // namespace

GLvoid InitGLFluid(GLsizei Width, GLsizei Height) {
 

  glClearColor(1, 1,1, 1);
  glShadeModel(GL_SMOOTH);
    
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
    
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
    
  glMatrixMode(GL_MODELVIEW);

  glLightfv(GL_LIGHT1, GL_POSITION,LightPosition);
  glEnable(GL_LIGHT1);                           
 
}

//  Draw a bounding box.
void DrawCube() {
  glPushMatrix(); 
  glTranslatef(boxXLength,boxYLength,boxZLength);
  glBegin(GL_LINE_STRIP);
  glColor3f(0.f,0.f,0.f);
  glVertex3f(-boxXLength,-boxYLength,-boxZLength);
  glVertex3f(boxXLength,-boxYLength,-boxZLength);
  glVertex3f(boxXLength,boxYLength,-boxZLength);
  glVertex3f(-boxXLength, boxYLength,-boxZLength);
  glVertex3f(-boxXLength,-boxYLength,-boxZLength);
  
  glVertex3f(-boxXLength,-boxYLength,boxZLength);
  glVertex3f(boxXLength,-boxYLength,boxZLength);
  glVertex3f(boxXLength,boxYLength,boxZLength);
  glVertex3f(-boxXLength,boxYLength,boxZLength);
  glVertex3f(-boxXLength,-boxYLength,boxZLength);
  glVertex3f(-boxXLength,boxYLength,boxZLength);
  glVertex3f(-boxXLength,boxYLength,-boxZLength);
  glVertex3f(boxXLength,boxYLength,-boxZLength);
  glVertex3f(boxXLength,boxYLength,boxZLength);
  glVertex3f(boxXLength,-boxYLength,boxZLength);
  glVertex3f(boxXLength,-boxYLength,-boxZLength);
  glEnd();
  if (!FLAGS_interactive) {
  glTranslatef(-0.1f, -0.1f, -0.1f);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
    // x axis is red
    glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
    glVertex3f(1.0f, 0.0f, 0.0f);

    // y axis is green 
    glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 10.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    
    // z axis is blue
    glColor4f(0.0f, 0.0f, 10.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glEnd();
  }
  glLineWidth(1.0f);

  glPopMatrix();
 
}

// Window to draw the coefficients of the fluid.
static void spectra_pre_display(void) {
  glViewport ( 0, 0, spectra_window_x, spectra_window_y);
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ();
  gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
  glClearColor ( 1.0f, 1.0f, 1.0f, 1.0f );
  glClear ( GL_COLOR_BUFFER_BIT );
}

static void reshape_func_extra(int width, int height) {
  glutSetWindow ( spectra_win_id );
  glutReshapeWindow ( width, height );
  spectra_window_x = width;
  spectra_window_y = height;
}

static void post_display ( void ) {
  glutSwapBuffers ();
}

void DrawCursor(void) {
  glBegin(GL_TRIANGLE_STRIP);
  glVertex2f(-0.5, 0);
  glVertex2f(0.5,0);
  glVertex2f(0, 0.866);
  glVertex2f(-0.5, 0);
  glEnd();
  glBegin(GL_QUADS);
  glVertex2f(-0.25,0);
  glVertex2f(-0.25,-0.5);
  glVertex2f(0.25,-0.5);
  glVertex2f(0.25,0);
  glEnd();
}

// Fluid window.
void fluid_display_func(void) {
  glClearColor(1.f, 1.f, 1.f,1.f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  glTranslatef(xt,yt,zt);
  float m[4][4];
  build_rotmatrix(m, _quat);
  glMultMatrixf(&m[0][0]);
  glTranslatef(-boxXLength,-boxYLength,-boxZLength);
  
  DrawCube();
  // Drawfluid
  if (draw_particles) {
    fluid->DrawParticles(FLAGS_particle_length_factor);
  }
  glPushMatrix();
  fluid->DrawObstacles();
  glPopMatrix();
  
  glPushMatrix();
  glTranslatef(boxXLength,boxYLength,boxZLength);
  fluid->DrawSmoke(aplha_multipler, low_cut);
  glPopMatrix();  
  if (FLAGS_interactive) {
    glPushMatrix();
    glColor3f(0.0,0.0,0.0);
    /*glBegin(GL_QUADS);
    glVertex2f(mouse_x_norm, mouse_y_norm);
    glVertex2f(mouse_x_norm+0.02, mouse_y_norm);
    glVertex2f(mouse_x_norm+0.02, mouse_y_norm + 0.02);
    glVertex2f(mouse_x_norm, mouse_y_norm+0.02);*/
    glTranslatef(mouse_x_norm, mouse_y_norm, 0);
    glRotatef(30, 0, 0, 1);
    glScalef(0.04, 0.04,0.04);
    DrawCursor();
    glEnd();
    glPopMatrix();
  }
  if (FLAGS_capture_video) {
      movie.addFrameGL(); // This function will eat lots of memory afterwards.
  }
  post_display();
}

// Window to draw the coefficients of the fluid.
static void spectra_display_func(void) {
  spectra_pre_display();
  // Draw the spectra coefficients.
  fluid->DrawCoefficients(spectra_draw_factor);

  // glPrintString1(0.00,0.95 ,printstr.c_str());
  if (FLAGS_capture_video) {
    spectra1Movie.addFrameGL();
  }

  post_display();
}

void keyPressed(unsigned char key, int x, int y) {
  /* avoid thrashing this procedure */
  usleep(100);

  switch (key) {    
  case ESCAPE: 
    SavePreviewAndQuit();            	
    break;
  case 76: 
  case 108: // switch the lighting.
    printf("L/l pressed; light is: %d\n", light);
    light = light ? 0 : 1;   
    printf("Light is now: %d\n", light);
    if (!light) {
      glDisable(GL_LIGHTING);
    } else {
      glEnable(GL_LIGHTING);
    }
    break;
  case 't':
    aplha_multipler *= 0.8;
    break;
  case 'T':
    aplha_multipler *= 1.2;
    break;
  case 'x': // cut more
    low_cut += 1;
    if (low_cut > 255) low_cut = 255; 
    break;
  case 'X': // cut less
    low_cut -= 1;
    if (low_cut < 0) low_cut = 0; 
    break;
  case 'c':
  case 'C':
    fluid->ClearDensity();
    break;
  case 'r':
  case 'R':
    fluid->ReSeedParticles();
    break;
  case 'v':
  case 'V':
    fluid->ResetCoeff();
    break;
  case 'p':
  case 'P':
    is_pause = !is_pause;
    break;
  case 'a':
  case 'A':
    add_smoke = !add_smoke;
    break;
  case 'd':
  case 'D':
    fluid->PrintDebugInfo();
    break;
  case 'w':
  case 'W':
    write_dens_pbrt = !write_dens_pbrt;
    break;
  case 'm':
  case 'M':
    draw_particles = !draw_particles;
    break;
  case '[':
    fluid->SetObstacleVelo(VEC3F(-0.1,0,0));
    break;
  case ']':
    fluid->SetObstacleVelo(VEC3F(0.1,0,0));
    break;
    
  default:
    break;
  }
}

GLvoid ReSizeGLScene(GLsizei Width, GLsizei Height) {
  if (Height==0)
    Height=1;

  glViewport(0, 0, Width, Height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,100.0f);
  glMatrixMode(GL_MODELVIEW);
  window_x_res = Width;
  window_y_res = Height;
}

void mouse_rotate(int x, int y) {
  if (!FLAGS_interactive) {
    float spin_quat[4];

    float sx2 = window_x_res*0.5f, sy2 = window_y_res*0.5f;

    trackball(spin_quat,
        (mouse_ax - sx2) / sx2,
        -(mouse_ay - sy2) / sy2,
        (x - sx2) / sx2,
        -(y - sy2) / sy2);
    add_quats(spin_quat, _quat, _quat);
    mouse_ax = x;
    mouse_ay = y;
  } else {
    float sx2 = window_x_res*0.5f, sy2 = window_y_res*0.5f;
    float cur_pos_x = (x - sx2) / sx2;
    float cur_pos_y = -(y - sy2) / sy2;
    mouse_x_norm = cur_pos_x*0.8 + boxXLength;
    mouse_y_norm = cur_pos_y*0.8 + boxYLength;
    
    float prev_pos_x =  (mouse_ax - sx2) / sx2;
    float prev_pos_y = -(mouse_ay - sy2) / sy2;
    fluid->SetOmega(cur_pos_x - prev_pos_x, cur_pos_y - prev_pos_y, cur_pos_x, cur_pos_y);
    mouse_ax = x;
    mouse_ay = y;
  }
}

void MouseFunc(int button, int state, int x, int y) {
  if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    mouse_ax = x;
    mouse_ay = y;
  }
}

void specialKeyPressed(int key, int x, int y)  {
  /* avoid thrashing this procedure */
  usleep(100);

  switch (key) {    
  case GLUT_KEY_PAGE_UP: // move the cube into the distance.
    zt-=0.1f;
    break;
  case GLUT_KEY_PAGE_DOWN: // move the cube closer.
    zt+=0.1f;
    break;
  case GLUT_KEY_UP:
    yt += 0.1f;
    break;
  case GLUT_KEY_DOWN:
    yt -= 0.1f;
    break;
  case GLUT_KEY_LEFT:
    xt -= 0.1f;
    break;
  case GLUT_KEY_RIGHT:
    xt += 0.1f;
    break;
  
  default:
  break;
  }	
}

static void idle_func ( void ) {
  if (!is_pause) {
    fluid->Step();
    if (add_smoke) {
      if (FLAGS_use_disk_smoke) {
        fluid->AddSmokeTestCaseCylinder();
      } else if (FLAGS_use_sphere_smoke) {
        fluid->AddSmokeTestCaseSphere();
      }  else {        
        fluid->AddSmokeTestCase();
      }
    }
    if (write_dens_pbrt) {
      fluid->OutputSmokeToFolder();
    }
    if (fluid->isFinished()) {
      SavePreviewAndQuit();
    }
  }
  glutSetWindow ( fluid_win_id );
  glutPostRedisplay ();
}

static void open_glut_window ( void ) {
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA);
  /* get a 640 x 480 window */
  glutInitWindowSize(FLAGS_window_x, FLAGS_window_y);  
  /* the window starts at the upper left corner of the screen */
  glutInitWindowPosition(0, 0);  

  /* Open a window */  
  fluid_win_id = glutCreateWindow("3D Viewer.");  

  glutDisplayFunc(&fluid_display_func);  
  glutIdleFunc(&idle_func); 
  glutReshapeFunc(&ReSizeGLScene);
  glutKeyboardFunc(&keyPressed);
  glutSpecialFunc(&specialKeyPressed);
  glutMotionFunc(mouse_rotate);
  glutMouseFunc(MouseFunc);
  
  glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );
  glutInitWindowPosition ( 1700, 0 );
  glutInitWindowSize ( spectra_window_x, spectra_window_x );
  spectra_win_id = glutCreateWindow("energy");
  // glClearColor ( 1.f, 1.f, 1.f, 1.0f );
  glClear ( GL_COLOR_BUFFER_BIT );
  glutSwapBuffers ();
  glClear ( GL_COLOR_BUFFER_BIT );
  glutSwapBuffers ();
  spectra_pre_display();
  glutReshapeFunc (reshape_func_extra);
  glutDisplayFunc(spectra_display_func);
}
/*

 **/
int main(int argc, char ** argv) {
  //google::ParseCommandLineFlags(&argc, &argv, true);
  //google::InitGoogleLogging(argv[0]);
  
  if (argc != 2)
  {
    std::cout << " Usage: " << argv[0] << " *.cfg" << std::endl;
    return 0;
  }
  SIMPLE_PARSER parser(argv[1]);
  
  FLAGS_use_obstacles = parser.getBool("use_obstacles", 0);
  FLAGS_move_obstacle = parser.getBool("move_obstacle", 0);
  FLAGS_obstacle_type = parser.getString("obstacle_type", "box");
  FLAGS_obstacle_file = parser.getString("obstacle_file", "");
  FLAGS_obstacle_list_file = parser.getString("obstacle_list_file", "");
  FLAGS_coefficients_file = parser.getString("coefficients_file", "");
  FLAGS_coefficients_file_in = parser.getString("coefficients_file_in", "");
  FLAGS_obstacle_force_scale =  parser.getFloat("obstacle_force_scale", 20.0);
  FLAGS_obstacle_implicit_method = parser.getBool("obstacle_implicit_method", 1.0);
  FLAGS_capture_video = parser.getBool("capture_video", 0);
  FLAGS_preview_fpath = parser.getString("preview_fpath", "");
  FLAGS_FFTW_threads = parser.getInt("FFTW_threads", 1);
  FLAGS_buoyancy_step = parser.getInt("buoyancy_step", -1);
  FLAGS_interactive = parser.getBool("interactive", 0);
  FLAGS_solver_type = parser.getString("solver_type", "BiCGSTAB");
  FLAGS_add_force_with_source = parser.getBool("add_force_with_source", 0);
  
  FLAGS_window_x = parser.getInt("window_x", 1024);
  FLAGS_window_y = parser.getInt("window_y", 1024);
  FLAGS_xRes = parser.getInt("xRes", 128);
  FLAGS_yRes = parser.getInt("yRes", 128);
  FLAGS_zRes = parser.getInt("zRes", 256);
  FLAGS_basis_dim = parser.getInt("basis_dim", 512);
  FLAGS_dt = parser.getFloat("dt", 0.1);
  FLAGS_buoyancy = parser.getFloat("buoyancy", 0.0001);
  FLAGS_buoyancy_direction = parser.getString("buoyancy_direction", "x");
  FLAGS_viscosity = parser.getFloat("viscosity", 0.005);
  FLAGS_added_smoke_density = parser.getFloat("added_smoke_density",2.0);
  FLAGS_total_num_particles = parser.getInt("total_num_particles", 1000);
  FLAGS_particle_length_factor = parser.getFloat("particle_length_factor", 0.1);
  FLAGS_total_frame = parser.getInt("total_frame", 2000);
  
  FLAGS_basis_type = parser.getString("basis_type", "all_dirichlet");
  FLAGS_constant_init_strategy = parser.getString("constant_init_strategy", "principle_x");
  FLAGS_basis_tensor_fname = parser.getString("basis_tensor_fname", "");
  FLAGS_density_folder = parser.getString("density_folder", "");
  FLAGS_source_xpos = parser.getFloat("source_xpos", 0);
  FLAGS_source_ypos = parser.getFloat("source_ypos", 0);
  FLAGS_source_zpos = parser.getFloat("source_zpos", 0);
  FLAGS_source_length = parser.getFloat("source_length", 0);
  FLAGS_source_width = parser.getFloat("source_width", 0);
  FLAGS_source_height = parser.getFloat("source_height", 0);
  FLAGS_attenuate_density = parser.getBool("attenuate_density", 0);
  FLAGS_add_density = parser.getBool("add_density", 1);
  FLAGS_density_attenuate_factor = parser.getFloat("density_attenuate_factor", 1.0);
  FLAGS_write_density_to_PBRT = parser.getBool("write_density_to_PBRT", 1);
  FLAGS_use_MacCormack = parser.getBool("use_MacCormack", 0);
  FLAGS_use_disk_smoke = parser.getBool("use_disk_smoke", 0);
  FLAGS_use_sphere_smoke = parser.getBool("use_sphere_smoke", 0);
  FLAGS_use_two_phase_smoke = parser.getBool("use_two_phase_smoke", 0);
  FLAGS_source_smoke_file = parser.getString("source_smoke_file", "");
  FLAGS_basisweightMultiplierC = parser.getFloat("basisweightMultiplierC", 0.0);
  
  glutInit (&argc, argv);
  
  // For problems of smaller size, use less threads are better.
   fftw_init_threads();
   fftw_plan_with_nthreads(FLAGS_FFTW_threads);
  
  const int max_lenght = std::max(std::max(FLAGS_xRes, FLAGS_yRes), FLAGS_zRes);
  boxXLength = static_cast<float>(FLAGS_xRes) / max_lenght * 0.5;
  boxYLength = static_cast<float>(FLAGS_yRes) / max_lenght * 0.5;
  boxZLength = static_cast<float>(FLAGS_zRes) / max_lenght * 0.5;
  int buoyancy_step = FLAGS_buoyancy_step;
  // Default to apply buoyancy through whole simulation.
  if (buoyancy_step < 0) {
    buoyancy_step = FLAGS_total_frame;
  }
  
  // Initialize components
  fluid = new LaplacianFluid3D(FLAGS_xRes, FLAGS_yRes, FLAGS_zRes,FLAGS_basis_dim, FLAGS_dt,
    FLAGS_viscosity, FLAGS_buoyancy, FLAGS_buoyancy_direction, FLAGS_added_smoke_density,
    FLAGS_total_num_particles, FLAGS_basis_type, FLAGS_constant_init_strategy ,
    lpos, FLAGS_basis_tensor_fname,
    FLAGS_total_frame, FLAGS_density_folder,
    FLAGS_solver_type,
    FLAGS_attenuate_density, FLAGS_density_attenuate_factor, FLAGS_use_MacCormack, FLAGS_coefficients_file,
    FLAGS_use_two_phase_smoke, buoyancy_step, FLAGS_coefficients_file_in, FLAGS_basisweightMultiplierC);
  
  ObstacleParams3D obsparam;
  InitializeObstacleParams(&obsparam);
  fluid->InitializeObstacles(obsparam);
  fluid->InitializeSourceSmoke(VEC3F(FLAGS_source_xpos, FLAGS_source_ypos, FLAGS_source_zpos),
                               VEC3F(FLAGS_source_length, FLAGS_source_width, FLAGS_source_height),
                               FLAGS_source_smoke_file, FLAGS_use_disk_smoke, FLAGS_add_force_with_source);
  
  write_dens_pbrt = FLAGS_write_density_to_PBRT;
  add_smoke = FLAGS_add_density;
  
  std::stringstream printstr_;
  printstr_ << "basis_dim: " << fluid->GetBasiDim(); 
  printstr = printstr_.str();
  
  open_glut_window();
  trackball(_quat, 0.0, 0.0, 0.0, 0.0);
  /* Initialize our window. */
  InitGLFluid(FLAGS_window_x, FLAGS_window_y);
  window_x_res = FLAGS_window_x;
  window_y_res = FLAGS_window_y;
  
  /* Start Event Processing Engine */  
  glutMainLoop();  
    
  return 0;
}
