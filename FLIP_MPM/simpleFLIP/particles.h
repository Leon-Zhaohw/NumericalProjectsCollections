#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include "grid.h"
#include "vec2.h"

struct Particles{
   Grid &grid;
   int np; // number of particles
   std::vector<Vec2f> x, u; // positions and velocities
   // transfer stuff
   Array2f sum;

   Particles(Grid &grid_)
      :grid(grid_), np(0),
       sum(grid_.pressure.nx+1, grid_.pressure.ny+1)
   {}

   void add_particle(const Vec2f &px, const Vec2f &pu);
   void transfer_to_grid(void);
   void update_from_grid(void);
   void move_particles_in_grid(float dt);
   void write_to_file(const char *filename_format, ...);

   private:
   template<class T> void accumulate(T &accum, float q, int i, int j, float fx, float fy);
};

#endif
