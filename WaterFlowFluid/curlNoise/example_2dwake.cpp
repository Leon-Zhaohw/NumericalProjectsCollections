#include "example_2dwake.h"
#include "bluenoise.h"

Example2DWake::
Example2DWake(void)
   : disc_centre(0.7,0.5),
     disc_radius(0.12),
     disc_influence(0.25),
     ymin(0),
     ymax(0.75),
     initial_xmin(1.5),
     initial_xmax(3.5),
     background_flow_speed(-0.25),
     wake_expansion(0.3)
{
   rigid.resize(1);
   rigid[0].centre=disc_centre;
   rigid[0].radius=disc_radius;
   noise_lengthscale.push_back(0.1);
   noise_gain.push_back(0.03);
   noise_lengthscale.push_back(0.06);
   noise_gain.push_back(0.03*0.35);
   noise_lengthscale.push_back(0.026);
   noise_gain.push_back(0.03*0.1);
}

float Example2DWake::
solid_distance(float x, float y) const
{
   return dist(Vec2f(x,y),disc_centre)-disc_radius;
}

bool Example2DWake::
seed_particles(std::vector<Vec2f> &x, float dt) const
{
   if(t==0){
      bluenoise_sample(0.01f, Vec2f(initial_xmin, ymin), Vec2f(initial_xmax, ymax), x);
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

float Example2DWake::
potential(float x, float y) const
{
   // begin with background laminary flow (adjusted so zero level goes through centre of disc)
   float p=-background_flow_speed*(y-disc_centre[1]);
   // modify for disc
   float d=solid_distance(x,y);
   if(d<disc_influence){
      p*=ramp(d/disc_influence);
   }
   // add turbulent wake that respects disc
   float wake_x=smooth_step(1-(x-disc_centre[0])/disc_radius);
   if(wake_x>0){
      float wake_y=smooth_step(1-std::fabs(y-disc_centre[1])/(1.5f*disc_radius+wake_expansion*(disc_centre[0]-x)));
      if(wake_y>0){
         float wake=wake_x*wake_y;
         float s=0;
         for(unsigned int i=0; i<noise_lengthscale.size(); ++i){
            s+=ramp(d/noise_lengthscale[i])*noise_gain[i]*noise((x-background_flow_speed*t)/noise_lengthscale[i], y/noise_lengthscale[i]);
         }
         p+=wake*s;
      }
   }
   return p;
}

void Example2DWake::
advance_time(float dt)
{
   t+=dt;
   noise.set_time(0.5f*noise_gain[0]/noise_lengthscale[0]*t);
}

