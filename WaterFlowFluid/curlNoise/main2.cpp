#include <cstdio>
#include "gluvi.h"
#include "example_2dmouse.h"
#include "example_2dwake.h"
#include "example_2dfancy.h"

float display_timestep=0.04;
float simulation_timestep=0.04;
Gluvi::PanZoom2D cam(0, 0, 0.75);
const unsigned int tracelength=3;
std::vector<std::vector<Vec2f> > xtrace(tracelength);
CurlNoise2 *example=0;

void display(void)
{
   // first draw particle traces
   glDisable(GL_LIGHTING);
   glLineWidth(1);
   glColor3f(0, 0, 0);
   for(unsigned int i=0; i<xtrace[0].size(); ++i){
      glBegin(GL_LINE_STRIP);
      for(unsigned int j=0; j<tracelength; ++j){
         if(i>=xtrace[j].size()) break;
         float grey=j/(tracelength+0.5f);
         glColor3f(grey, grey, grey);
         glVertex2fv(xtrace[j][i].v);
      }
      glEnd();
   }
   glPointSize(2);
   glColor3f(0, 0, 0);
   glBegin(GL_POINTS);
   for(unsigned int i=0; i<xtrace[0].size(); ++i){
      glVertex2fv(xtrace[0][i].v);
   }
   glEnd();
   // next draw any solids
   glLineWidth(2);
   glColor3f(0.5,0,0);
   for(unsigned int r=0; r<example->rigid.size(); ++r){
      glBegin(GL_LINES);
      for(unsigned int e=0; e<example->rigid[r].edge.size(); ++e){
         unsigned int i, j; assign(example->rigid[r].edge[e], i, j);
         glVertex2fv(example->rigid[r].world_x(i).v);
         glVertex2fv(example->rigid[r].world_x(j).v);
      }
      glEnd();
      if(example->rigid[r].radius>0){
         glBegin(GL_LINE_LOOP);
         const int n=30;
         for(int i=0; i<n; ++i){
            glVertex2fv((example->rigid[r].centre+example->rigid[r].radius*Vec2f(std::cos(i*2*M_PI/n), std::sin(i*2*M_PI/n))).v);
         }
         glEnd();
      }
   }
}

void timer(int junk)
{
   Vec2f v, midx;
   for(unsigned int j=tracelength-1; j>0; --j){
      xtrace[j]=xtrace[j-1];
   }
   for(unsigned int i=0; i<xtrace[0].size(); ++i){
      // almost RK2 (not updating t halfway through time step, so it's technically only 1st order accurate)
      example->get_velocity(xtrace[0][i], v);
      midx=xtrace[0][i]+0.5f*simulation_timestep*v;
      example->get_velocity(midx, v);
      xtrace[0][i]+=simulation_timestep*v;
   }
   example->seed_particles(xtrace[0], simulation_timestep);
   example->advance_time(simulation_timestep);
   glutPostRedisplay();
   glutTimerFunc(int(1000*display_timestep), timer, 0);
}

void mouse(int button, int state, int x, int y)
{
   if(state==GLUT_DOWN){
      Vec2f mouse_x;
      cam.transform_mouse(x, y, mouse_x.v);
      example->click(mouse_x[0], mouse_x[1]);
   }
}

void drag(int x, int y)
{
   Vec2f mouse_x;
   cam.transform_mouse(x, y, mouse_x.v);
   example->click(mouse_x[0], mouse_x[1]);
}

int main(int argc, char **argv)
{
   Gluvi::init("curl noise in 2D", &argc, argv);
   glClearColor(1,1,1,0);

   Gluvi::camera=&cam;
   Gluvi::userDisplayFunc=display;
   Gluvi::userMouseFunc=mouse;
   Gluvi::userDragFunc=drag;
   glutTimerFunc(int(1000*display_timestep), timer, 0);

   //example=new Example2DMouse();
   //example=new Example2DWake();
   example=new Example2DFancy();

   example->seed_particles(xtrace[0], simulation_timestep);

   Gluvi::run();
   return 0;
}

