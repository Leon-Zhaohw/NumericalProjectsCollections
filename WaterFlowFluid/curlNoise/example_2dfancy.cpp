#include "example_2dfancy.h"
#include "bluenoise.h"

Example2DFancy::
Example2DFancy(void)
   : initial_ymin(-0.5),
     initial_ymax(1.2),
     initial_xmin(0.75),
     initial_xmax(3.5),
     wind_velocity(-0.3,0),
     envelope(0.17),
     noise_lengthscale(0.1),
     noise_gain(0.03)
{
   rigid.resize(4);

   rigid[0].add_strip(4, -0.05,-0.1, 0.,-0.1, 0.,0.1, 0.05,0.1);
   rigid[0].set_velocity(Vec2f(0,0), 2);
   rigid[0].set_pose(Vec2f(0.5,0.48), 0.3);

   rigid[1].add_strip(4, -0.1,-0.1, 0.,0.1, 0.1,-0.1, -0.1,-0.1);
   rigid[1].set_velocity(Vec2f(0,0), 1.5);
   rigid[1].set_pose(Vec2f(0.5,0.22), 0);

   rigid[2].add_strip(2, -0.1,0., 0.1,0.);
   rigid[2].set_velocity(Vec2f(0.1,0), 4);
   rigid[2].set_pose(Vec2f(0.65,0.35), -0.2);

   rigid[3].add_strip(2, 0.,-0.1, 0.,0.1);
   rigid[3].set_velocity(Vec2f(0.1,0), 4);
   rigid[3].set_pose(Vec2f(0.65,0.35), -0.2);
}

bool Example2DFancy::
seed_particles(std::vector<Vec2f> &x, float dt) const
{
   if(t==0){
      bluenoise_sample(0.015f, Vec2f(initial_xmin, initial_ymin), Vec2f(initial_xmax, initial_ymax), x);
      return true;
   }else
      return false; // no particles added after initial seeding at t=0
}

float Example2DFancy::
potential(float x, float y) const
{
   Vec2f px(x,y);
   float numer=(1/sqr(envelope))*cross(px,wind_velocity);
   float denom=1/sqr(envelope);
   float minphi=1e36;
   for(unsigned int r=0; r<rigid.size(); ++r){
      Vec2f dx=px-rigid[r].centre;
      float psi=cross(px,rigid[r].velocity) + (sqr(envelope)-mag2(dx))/2*rigid[r].angular_velocity;
      float phi=rigid[r].distance(px);
      if(phi<minphi) minphi=phi;
      float m=1/(sqr(phi)+1e-18);
      numer+=m*psi;
      denom+=m;
   }
   float base_psi=numer/denom;

   float d=minphi/noise_lengthscale;
   float g=smooth_step(1.1-x);
   float turb=g*ramp(d)*noise_gain*noise((px-t*wind_velocity)/noise_lengthscale);

   return base_psi+turb;
}

void Example2DFancy::
advance_time(float dt)
{
   for(unsigned int r=0; r<rigid.size(); ++r){
      rigid[r].update_pose(dt);
   }
   t+=dt;
   noise.set_time(0.5f*noise_gain/noise_lengthscale*t);
}

