#ifndef EXAMPLE_3DPLUME_H
#define EXAMPLE_3DPLUME_H

#include "curlnoise.h"
#include "noise.h"

struct Example3DPlume: public CurlNoise3
{
   Vec3f sphere_centre;
   float sphere_radius;
   float plume_height;
   std::vector<float> noise_lengthscale;
   std::vector<float> noise_gain;
   FlowNoise3 noise;
   float ring_radius;
   float ring_speed;
   float rings_per_second;
   float ring_magnitude;
   float ring_falloff;
   mutable unsigned int seed;
   float particles_per_second;
   float seed_radius;
   float initial_band;

   Example3DPlume(void);

   // components of vector noise
   float noise0(float x, float y, float z) const { return noise(x, y, z); }
   float noise1(float x, float y, float z) const { return noise(y+31.416f, z-47.853f, x+12.793f); }
   float noise2(float x, float y, float z) const { return noise(z-233.145f, x-113.408f, y-185.31f); }

   float distance_and_normal(float x, float y, float z, Vec3f &normal) const;
   void match_boundary(Vec3f &psi, float d, float lengthscale, const Vec3f &normal) const;

   // virtual functions overriding base class
   bool seed_particles(std::vector<Vec3f> &x, float dt) const;
   Vec3f potential(float x, float y, float z) const;
   void advance_time(float dt);
};

#endif
