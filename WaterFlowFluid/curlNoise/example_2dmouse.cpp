#include "example_2dmouse.h"
#include "bluenoise.h"

Example2DMouse::
Example2DMouse(void)
   : ymin(0), ymax(0.75),
     xmin(0), xmax(1),
     mouse_x(0), mouse_y(0),
     click_time(-1e37),
     mouse_radius(0.5),
     duration(5),
     noise_lengthscale(0.14),
     noise_gain(0.09)
{
   rigid.resize(3);

   rigid[0].centre=Vec2f(0.7, 0.6);
   rigid[0].radius=0.17;

   rigid[1].centre=Vec2f(0.18, 0.12);
   rigid[1].radius=0.13;

   rigid[2].x.resize(4);
   rigid[2].x[0]=Vec2f(xmin, ymin);
   rigid[2].x[1]=Vec2f(xmax, ymin);
   rigid[2].x[2]=Vec2f(xmax, ymax);
   rigid[2].x[3]=Vec2f(xmin, ymax);
   rigid[2].edge.resize(4);
   rigid[2].edge[0]=Vec2ui(0,1);
   rigid[2].edge[1]=Vec2ui(1,2);
   rigid[2].edge[2]=Vec2ui(2,3);
   rigid[2].edge[3]=Vec2ui(3,0);
}

float Example2DMouse::
solid_distance(float x, float y) const
{
   float d=1e37;
   for(unsigned int i=0; i<rigid.size(); ++i){
      d=min(d, rigid[i].distance(Vec2f(x,y)));
   }
   if(d<0) d=0; // treat interiors of solids as solid, not fluid
   return d;
}

bool Example2DMouse::
seed_particles(std::vector<Vec2f> &x, float dt) const
{
   if(t==0){
      bluenoise_sample(0.012f, Vec2f(xmin, ymin), Vec2f(xmax, ymax), x);
      for(int i=0; i<(int)x.size(); ++i){
         if(solid_distance(x[i][0], x[i][1])<=0){
            erase_unordered(x, i);
            --i;
         }
      }
      return true;
   }else
      return false; // no particles added after initial seeding at t=0
}

float Example2DMouse::
potential(float x, float y) const
{
   float p=noise_gain*smooth_step(1-(t-click_time)/duration);
   if(p==0) return 0;
   float distance_to_mouse=std::sqrt(sqr(x-mouse_x)+sqr(y-mouse_y));
   p*=smooth_step(1-distance_to_mouse/mouse_radius);
   if(p==0) return 0;
   p*=noise(x/noise_lengthscale, y/noise_lengthscale);
   // modify for boundaries
   float d=solid_distance(x,y);
   if(d<noise_lengthscale){
      p*=ramp(d/noise_lengthscale);
   }
   return p;
}

void Example2DMouse::
advance_time(float dt)
{
   t+=dt;
   noise.set_time(0.5f*noise_gain/noise_lengthscale*t);
}

