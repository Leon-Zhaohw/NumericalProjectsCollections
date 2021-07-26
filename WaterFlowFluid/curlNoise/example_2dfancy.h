#ifndef EXAMPLE_2DFANCY_H
#define EXAMPLE_2DFANCY_H

#include "curlnoise.h"
#include "noise.h"

struct Example2DFancy: public CurlNoise2
{
   float initial_ymin, initial_ymax;
   float initial_xmin, initial_xmax;
   Vec2f wind_velocity;
   float envelope;
   float noise_lengthscale;
   float noise_gain;
   FlowNoise2 noise;

   Example2DFancy(void);

   // virtual functions overriding base class
   bool seed_particles(std::vector<Vec2f> &x, float dt) const;
   float potential(float x, float y) const;
   void advance_time(float dt);
};

#endif
