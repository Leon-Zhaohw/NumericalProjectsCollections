#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cfloat>

#include "gluvi.h"
#include "wallclocktime.h"
#include "openglutils.h"
#include "trimesh2d.h"
#include "dualfluidsim.h"
#include "sampleseeder.h"
#include "framestepper.h"
#include "triangle_wrapper.h"
#include "dist_funcs.h"
#include "geom_routines.h"
#include "marching_triangles.h"
#include "interpolator.h"
#include "scenes.h"

using namespace std;

//Simulation parameters
//---------------------
//toggle whether to dump out frames
#define MAKE_MOVIE


char output_location[256];


//estimate of the sim resolution - all other quantities computed relative to this.
float approx_dx = 0.02;  

//the size of the simulation box. (actual mesh domain size will be larger than this)
Vec2f solid_bounds_centre(0.5f, 0.5f );
Vec2f solid_bounds_dims(0.6f, 0.6f);

//volume conservation - use El Topo to correct for gradually accumulated volume loss, through small vertex adjustments
bool conserve_volume = true;

//screen dimensions
int win_height_pixels = 940;
int win_width_pixels = 560;

//length of each frame, and number of frames to simulate
float frame_length = 0.02f;
int max_frame = 2000;

//choose which pre-set scene to run, from the list below:
//stillpool_and_circle, sinewaveandcircle, twocirclesbeside, twocircles, bigcircle, sinusoid
float (*active_scene)(const Vec2f& pt) = sinewaveandcircle; 



//Simulation data
//---------------------
FrameStepper* stepper;

std::vector<Vec2f> g_vertices;
std::vector<Vec3ui> g_tris;

DualFluidSim* g_dualsim = NULL;
eltopo2d::SurfTrack* g_explicit_surface = NULL;

bool g_playing = false;
Vec2d domain_min, domain_max;

//Gluvi (UI) stuff
//---------------------
//adjust the camera position
Gluvi::PanZoom2D cam( 0.15, -0.02, 0.7);

double oldmousetime;
Vec2f oldmouse;
int frame = 0;

void display();
void mouse(int button, int state, int x, int y);
void drag(int x, int y);
void timer(int junk);

void keyboard(unsigned char key, int, int )
{
   if ( key == 'n' )
   {
      int junk = 0;
      timer( junk );
   }
   
   if ( key == ' ' )
   {
      g_playing = !g_playing;
      if ( g_playing )
      {
         timer( ~0 );
      }
   }
   
   glutPostRedisplay();
   
}

void setup_simulation() 
{
   std::cout << "Creating surface mesh." << std::endl;
  
   //This is a bit obtuse, but we'll get a grid version of our input signed distance field (active_scene)
   //then using marching_triangles to mesh it, then feed that into our surface tracker as the input surface.
   
   //Build grid-based distance field
   int divisions = 48;
   Array2f phi(divisions,divisions);
   Vec2f origin = Vec2f(domain_min);
   
   float max_domain_dim = max(domain_max[0] - domain_min[0], domain_max[1] - domain_min[1]);
   float march_dx = max_domain_dim/(float)divisions;
   for(int j = 0; j < divisions; ++j) {
      for(int i = 0; i < divisions; ++i) {
         Vec2f pos = origin + Vec2f(i,j)*march_dx;
         phi(i,j) = active_scene(pos);
      }  
   }

   //Surface mesh it.
   MarchingTriangles marching_triangles( origin, march_dx );
   marching_triangles.phi = phi;
   marching_triangles.contour_grid( );
   
   eltopo2d::SurfTrackInitializationParameters params;   
   params.m_proximity_epsilon = 1e-5;
   params.m_verbose = false;
   params.m_collision_safety = true;
   params.m_max_area_change = 0.04 * approx_dx * approx_dx;
   params.m_min_edge_length = 0.5 * approx_dx;
   params.m_max_edge_length = 1.25 * approx_dx;
   params.m_merge_proximity_epsilon = 1e-4;
   params.m_allow_topology_changes = true;
   params.m_perform_improvement = true;
   
   //Create the surface mesh from the marching results
   std::vector<Vec2d> xs( marching_triangles.x.size() );
   for ( unsigned int i = 0; i < xs.size(); ++i ) { xs[i] = Vec2d(marching_triangles.x[i]); }
   std::vector<Vec2ui> es( marching_triangles.edge.size() );
   for ( unsigned int i = 0; i < es.size(); ++i ) { es[i] = Vec2ui(marching_triangles.edge[i]); }
   std::vector<double> masses( marching_triangles.x.size(), 1.0 );
  

   //Initialize surface tracker
   g_explicit_surface = new eltopo2d::SurfTrack( xs, es, masses, params );
   g_explicit_surface->rebuild_static_broad_phase();

   
   std::cout << "Creating initial simulation mesh" << std::endl;
   
   TriMesh2D mesh;
     
   std::vector<int> expected_signs;
   
   SampleSeeder::generate_regular_samples( domain_min, domain_max, 2*approx_dx, g_vertices, expected_signs);
   SampleSeeder::generate_samples( *g_explicit_surface, approx_dx, g_vertices, expected_signs);
   compute_delaunay(g_vertices, g_tris);
   mesh.setup_mesh( g_vertices, g_tris );

   g_dualsim = new DualFluidSim( *g_explicit_surface );
   
   g_dualsim->setup( mesh, approx_dx, solid_bounds_centre, solid_bounds_dims, Vec2f(domain_min), Vec2f(domain_max));
   g_dualsim->conserve_volume = conserve_volume;
   g_dualsim->compute_liquidphi(expected_signs);
   g_dualsim->extrapolate_liquid_phi_into_solid();
   
   std::cout << "Done simulation setup." << std::endl;
   
}


int main(int argc, char **argv)
{

   Gluvi::init("Voronoi Fluid Simulator 2D", &argc, argv);
   
   Gluvi::winwidth=win_width_pixels;
   Gluvi::winheight=win_height_pixels;
   Gluvi::camera=&cam;
   Gluvi::userDisplayFunc=display;
   Gluvi::userMouseFunc=mouse;
   Gluvi::userDragFunc=drag;

   stepper = new FrameStepper(frame_length, max_frame);
   
   domain_min = Vec2d(solid_bounds_centre - 0.5f*solid_bounds_dims - 10.0f*Vec2f(approx_dx,approx_dx));
   domain_max = Vec2d(solid_bounds_centre + 0.5f*solid_bounds_dims + 10.0f*Vec2f(approx_dx,approx_dx));
   
   setup_simulation();

   set_time_base();
   glClearColor(1,1,1,1);
   glutKeyboardFunc(keyboard);
   glutTimerFunc(1000, timer, 0);
   
   if ( argc > 1 )
   {
      strncpy( output_location, argv[1], 256 );
   }
   else
   {
#ifdef WIN32
      strncpy( output_location, "c:/output/", 256 );
#else
      strncpy( output_location, "/var/tmp/voronoi2d/", 256 );
#endif     
   }

   
   std::cout << "Ouput location: " << output_location << std::endl;
   
   std::cout << "Beginning sim.  Press space to pause/resume." << std::endl;
   
   g_playing = true; //by default, set the sim running. Stop with space bar.

   Gluvi::run();

   return 0;
}


void draw_dual_sim();


void display(void)
{
   
   //Solid boundary

   Vec2f minx = solid_bounds_centre - 0.5f * Vec2f( solid_bounds_dims[0], solid_bounds_dims[1] );
   Vec2f maxx = solid_bounds_centre + 0.5f * Vec2f( solid_bounds_dims[0], solid_bounds_dims[1] );
   
   float saved_line_width;
   glGetFloatv( GL_LINE_WIDTH, &(saved_line_width) );
   
   glLineWidth( 3.0 );
   glColor3f(0,0,0);
   glBegin( GL_LINE_STRIP );
   glVertex2fv( minx.v );
   glVertex2f( maxx[0], minx[1] );
   glVertex2f( maxx[0], maxx[1] );   
   glVertex2f( minx[0], maxx[1] );   
   glVertex2fv( minx.v );
   glEnd();
   glLineWidth( saved_line_width );

   //Paint over the exterior to hide the outer mesh
   glColor3f(1,1,1);
   draw_box2d(minx - Vec2f(2,2), 2, 10);
   draw_box2d(minx - Vec2f(2,2), 10, 2);
   draw_box2d(minx + Vec2f(0,solid_bounds_dims[1]), 10, 2);
   draw_box2d(minx + Vec2f(solid_bounds_dims[0],0), 2, 10);
   
   //Simulation
   draw_dual_sim();
  
   
}


void draw_dual_sim() 
{
   const unsigned int SURFACE_VERTEX_SIZE = 8;
   const unsigned int SURFACE_SEGMENT_WIDTH = 4;
   const unsigned int VORONOI_EDGE_WIDTH = 2;
   const unsigned int VORONOI_SITE_SIZE = 4;
   
   glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);	// Make round points, not square points (doesn't work on certain platforms)
   
   glEnable(GL_LINE_SMOOTH);
   glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);		// Antialias the lines
   glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


   const std::vector<Vec2d>& xs = g_explicit_surface->m_positions;
   
   glPointSize(SURFACE_VERTEX_SIZE);
   glColor3f(0,0,0);

   //glBegin( GL_POINTS );
   for ( unsigned int i = 0; i < xs.size(); ++i )
   {
      draw_circle2d((Vec2f)xs[i], 0.003f, 20);
      //glVertex2dv( xs[i].v );
   }
   //glEnd();
   
   
   const std::vector<Vec2ui>& es = g_explicit_surface->m_mesh.edges;
   glColor3f(0,0,1);
   glLineWidth(SURFACE_SEGMENT_WIDTH);
   glBegin( GL_LINES );
   for ( unsigned int i = 0; i < es.size(); ++i )
   {
      if( es[i][0] == es[i][1] )
         continue;
      glVertex2dv( xs[es[i][0]].v );
      glVertex2dv( xs[es[i][1]].v );
   }
   glEnd();
   
 
   
   glColor3f( 0.3f, 0.3f, 0.3f );
   glLineWidth(VORONOI_EDGE_WIDTH);
   
   glBegin(GL_LINES);
   for(unsigned int i = 0; i < g_dualsim->mesh.edges.size(); ++i) 
   {      
      Vec2i edge = g_dualsim->mesh.edge_to_tri_map[i];
      if(edge[0] == -1 || edge[1] == -1)
         continue;
      
      Vec2f v0 = g_dualsim->mesh.tri_circumcentres[edge[0]];
      Vec2f v1 = g_dualsim->mesh.tri_circumcentres[edge[1]];      
      glVertex2fv(v0.v);
      glVertex2fv(v1.v);
   }
   glEnd();

   glPointSize(VORONOI_SITE_SIZE);
   //glBegin(GL_POINTS);
   for(unsigned int i = 0; i < g_dualsim->mesh.vertices.size(); ++i) 
   {
      if ( g_dualsim->liquid_phi[i] < 0.0 )
      {
         glColor3f(0,0,1);
      }
      else
      {
         glColor3f(1,0,0);         
      }
      draw_circle2d((Vec2f)g_dualsim->mesh.vertices[i], 0.003f, 20);
      
      //glVertex2fv( g_dualsim->mesh.vertices[i].v );
   }
   //glEnd();
   
}



void mouse(int button, int state, int x, int y)
{
   Vec2f newmouse;
   cam.transform_mouse(x, y, newmouse.v);
   double newmousetime=get_time_in_seconds();

   oldmouse=newmouse;
   oldmousetime=newmousetime;
   
   Gluvi::camera->click( button, state, x, y );
   
   glutPostRedisplay();
   
}

void drag(int x, int y)
{
   Vec2f newmouse;
   cam.transform_mouse(x, y, newmouse.v);
   double newmousetime=get_time_in_seconds();

   Gluvi::camera->drag( x, y );
   
   oldmouse=newmouse;
   oldmousetime=newmousetime;
}


void timer(int junk)
{
   
#ifdef MAKE_MOVIE

      char sgifileformat[256];
      sprintf( sgifileformat, "%s/screenshot%%04d.sgi", output_location );
      Gluvi::sgi_screenshot(sgifileformat, frame);    
   
#endif

   std::cout << "--------------------------- Frame " << frame << " -----------------------------" << std::endl;
   
   if(!stepper->done_simulation()) 
   {
      
      while(!stepper->done_frame()) 
      {
         float cfl_step = g_dualsim->get_cfl_step();
         float t = stepper->get_step_length(cfl_step);
         
         std::cout << "Taking substep of length: " << t << std::endl;
         
         g_dualsim->advance( t );
         stepper->advance_step( t );
      }
      
      frame++;
      stepper->advance_frame();
   }
   
   if(stepper->done_simulation())
      exit(0);
   
   glutPostRedisplay();
   
   if ( g_playing )
   {
      if(!stepper->done_simulation())
         glutTimerFunc( 1, timer, 0);
   }
   
}



