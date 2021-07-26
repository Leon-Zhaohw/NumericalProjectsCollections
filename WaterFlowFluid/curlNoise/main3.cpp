#include <cstdio>
#include "gluvi.h"
#include "example_3dplume.h"

float display_timestep=0.05;
float simulation_timestep=0.025;
Gluvi::Target3D cam;
const unsigned int tracelength=3;
std::vector<std::vector<Vec3f> > xtrace(tracelength);
Example3DPlume *example=0;

void set_view(void)
{
}

void display(void)
{
   // first draw particle traces
   glDisable(GL_LIGHTING);
   glLineWidth(1);
   for(unsigned int i=0; i<xtrace[0].size(); ++i){
      glBegin(GL_LINE_STRIP);
      for(unsigned int j=0; j<tracelength; ++j){
         if(i>=xtrace[j].size()) break;
         float grey=(tracelength-j)/(float)tracelength;
         glColor3f(grey, grey, grey);
         glVertex3fv(xtrace[j][i].v);
      }
      glEnd();
   }
   glPointSize(1);
   glColor3f(1, 1, 1);
   glBegin(GL_POINTS);
   for(unsigned int i=0; i<xtrace[0].size(); ++i){
      glVertex3fv(xtrace[0][i].v);
   }
   glEnd();

   // next draw any solids
   glEnable(GL_LIGHTING);
   Gluvi::set_generic_lights();

   Gluvi::set_matte_material(.5,.5,.3);
   glBegin(GL_QUADS);
   glNormal3f(0, 1, 0);
   glVertex3f(-2, -0.01, -2);
   glVertex3f(2, -0.01, -2);
   glVertex3f(2, -0.01, 2);
   glVertex3f(-2, -0.01, 2);
   glEnd();

   Gluvi::set_matte_material(.4,.4,.8);
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glTranslatef(example->sphere_centre[0], example->sphere_centre[1], example->sphere_centre[2]);
   glutSolidSphere(example->sphere_radius, 20, 10);
   glPopMatrix();
}

void timer(int junk)
{
   Vec3f v, midx;
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

int main(int argc, char **argv)
{
   Gluvi::init("curl noise in 3D", &argc, argv);

   cam.dist=5;
   cam.target[0]=0.5;
   cam.target[1]=0.5;
   cam.target[2]=0.5;
   Gluvi::camera=&cam;
   Gluvi::userDisplayFunc=display;
   glutTimerFunc(int(1000*display_timestep), timer, 0);

   example=new Example3DPlume();
   example->seed_particles(xtrace[0], simulation_timestep);

   Gluvi::run();
   return 0;
}

