#ifndef NOISE_H
#define NOISE_H

#include "vec.h"

struct Noise2
{
   Noise2(unsigned int seed=171717);
   void reinitialize(unsigned int seed);
   float operator()(float x, float y) const;
   float operator()(const Vec2f &x) const { return (*this)(x[0], x[1]); }

   protected:
   static const unsigned int n=128;
   Vec2f basis[n];
   int perm[n];

   unsigned int hash_index(int i, int j) const
   { return perm[(perm[i%n]+j)%n]; }
};

struct Noise3
{
   Noise3(unsigned int seed=171717);
   void reinitialize(unsigned int seed);
   float operator()(float x, float y, float z) const;
   float operator()(const Vec3f &x) const { return (*this)(x[0], x[1], x[2]); }

   protected:
   static const unsigned int n=128;
   Vec3f basis[n];
   int perm[n];

   unsigned int hash_index(int i, int j, int k) const
   { return perm[(perm[(perm[i%n]+j)%n]+k)%n]; }
};

struct Noise4
{
   Noise4(unsigned int seed=171717);
   void reinitialize(unsigned int seed);
   float operator()(float x, float y, float z, float t) const;
   float operator()(const Vec4f &x) const { return (*this)(x[0], x[1], x[2], x[3]); }

   protected:
   static const unsigned int n=128;
   Vec4f basis[n];
   int perm[n];

   unsigned int hash_index(int i, int j, int k, int l) const
   { return perm[(perm[(perm[(perm[i%n]+j)%n]+k)%n]+l)%n]; }
};

// FlowNoise classes - time varying versions of some of the above

struct FlowNoise2: public Noise2
{
   FlowNoise2(unsigned int seed=171717, float spin_variation=0.2);
   void set_time(float t); // period of repetition is approximately 1

   protected:
   Vec2f original_basis[n];
   float spin_rate[n];
};

struct FlowNoise3: public Noise3
{
   FlowNoise3(unsigned int seed=171717, float spin_variation=0.2);
   void set_time(float t); // period of repetition is approximately 1

   protected:
   Vec3f original_basis[n];
   float spin_rate[n];
   Vec3f spin_axis[n];
};

#endif
