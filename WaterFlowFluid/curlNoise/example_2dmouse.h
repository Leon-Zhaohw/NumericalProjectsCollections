#ifndef EXAMPLE_2DMOUSE_H
#define EXAMPLE_2DMOUSE_H

#include "curlnoise.h"
#include "noise.h"

struct Example2DMouse: public CurlNoise2
{
   std::vector<Vec2f> disc_centre;
   std::vector<float> disc_radius;
   float ymin, ymax;
   float xmin, xmax;
   float mouse_x, mouse_y;
   float click_time;
   float mouse_radius;
   float duration;
   float noise_lengthscale;
   float noise_gain;
   FlowNoise2 noise;

   Example2DMouse(void);
   float solid_distance(float x, float y) const;

   // virtual functions overriding base class
   bool seed_particles(std::vector<Vec2f> &x, float dt) const;
   float potential(float x, float y) const;
   void advance_time(float dt);

   void click(float x, float y)
   {
      mouse_x=x;
      mouse_y=y;
      click_time=t;
   }
};

#endif
