#ifndef CURLNOISE_H
#define CURLNOISE_H

// base classes for 2d and 3d examples

#include "vec.h"
#include "rigidshape2.h"

struct CurlNoise2
{
   float t; // time
   float delta_x; // used for finite difference approximations of curl
   std::vector<RigidShape2> rigid; // list of rigid shapes in flow

   CurlNoise2(void)
      : t(0), delta_x(1e-4)
   {}

   virtual ~CurlNoise2(void)
   {}

   virtual bool seed_particles(std::vector<Vec2f> &x, float dt) const
   { return false; }

   virtual float potential(float x, float y) const
   { return 0; }

   virtual void advance_time(float dt)
   {
      t+=dt;
   }

   virtual void click(float x, float y)
   {}

   void get_velocity(const Vec2f &x, Vec2f &v) const
   {
      v[0]=-(potential(x[0], x[1]+delta_x) - potential(x[0], x[1]-delta_x))/(2*delta_x);
      v[1]= (potential(x[0]+delta_x, x[1]) - potential(x[0]-delta_x, x[1]))/(2*delta_x);
   }
};

struct CurlNoise3
{
   float t; // time
   float delta_x; // used for finite difference approximations of curl

   CurlNoise3(void)
      : t(0), delta_x(1e-4)
   {}

   virtual ~CurlNoise3(void)
   {}

   virtual bool seed_particles(std::vector<Vec3f> &x, float dt) const
   { return false; }

   virtual Vec3f potential(float x, float y, float z) const
   { return 0; }

   virtual void advance_time(float dt)
   {
      t+=dt;
   }

   void get_velocity(const Vec3f &x, Vec3f &v) const
   {
      v[0]=( (potential(x[0], x[1]+delta_x, x[2])[2] - potential(x[0], x[1]-delta_x, x[2])[2])
            -(potential(x[0], x[1], x[2]+delta_x)[1] - potential(x[0], x[1], x[2]-delta_x)[1]) ) / (2*delta_x);
      v[1]=( (potential(x[0], x[1], x[2]+delta_x)[0] - potential(x[0], x[1], x[2]-delta_x)[0])
            -(potential(x[0]+delta_x, x[1], x[2])[2] - potential(x[0]-delta_x, x[1], x[2])[2]) ) / (2*delta_x);
      v[2]=( (potential(x[0]+delta_x, x[1], x[2])[1] - potential(x[0]-delta_x, x[1], x[2])[1])
            -(potential(x[0], x[1]+delta_x, x[2])[0] - potential(x[0], x[1]-delta_x, x[2])[0]) ) / (2*delta_x);
   }
};

#endif
