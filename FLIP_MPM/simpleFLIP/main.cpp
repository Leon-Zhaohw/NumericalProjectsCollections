#include <cstdio>
#include "particles.h"
#include "util.h"

using namespace std;

float fluidphi(Grid &grid, float x, float y)
{
   //return y-0.5*grid.ly;
   //return min(sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.625*grid.ly))-0.02*grid.ly, y-0.6*grid.ly);
   //return min(sqrt(sqr(x-0.3333*grid.lx)+sqr(y-0.71*grid.ly))-0.3*grid.ly, y-0.2*grid.ly);
   //return max(y-0.8*grid.ly, -sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.2*grid.ly))+0.1*grid.lx);
   //return sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.75*grid.ly))-0.15*grid.lx;
   return min(y-0.05*grid.ly, sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.5*grid.ly))-0.05*grid.lx);
   //return 0.75*grid.lx-x;
   //return max(x-0.75*grid.lx, 0.25*grid.lx-x, y-0.75*grid.ly, 0.25*grid.ly-y);
}

void project(Grid &grid, float &x, float &y, float current, float target)
{
   float dpdx=(fluidphi(grid, x+1e-4, y)-fluidphi(grid, x-1e-4, y))/2e-4;
   float dpdy=(fluidphi(grid, x, y+1e-4)-fluidphi(grid, x, y-1e-4))/2e-4;
   float scale=(target-current)/sqrt(dpdx*dpdx+dpdy*dpdy);
   x+=scale*dpdx;
   y+=scale*dpdy;
}

void init_water_drop(Grid &grid, Particles &particles, int na, int nb)
{
   int i, j, a, b;
   float x, y, phi;

   for(i=1; i<grid.marker.nx-1; ++i){
      for(j=1; j<grid.marker.ny-1; ++j){
         for(a=0; a<na; ++a){
            for(b=0; b<nb; ++b){
               x=(i+(a+0.1+0.8*rand()/(double)RAND_MAX)/na)*grid.h;
               y=(j+(b+0.1+0.8*rand()/(double)RAND_MAX)/nb)*grid.h;
               phi=fluidphi(grid, x, y);
               if(phi>-0.25*grid.h/na)
                  continue;
               else if(phi>-1.5*grid.h/na){
                  project(grid, x, y, phi, -0.75*grid.h/na);
                  phi=fluidphi(grid, x, y);
                  project(grid, x, y, phi, -0.75*grid.h/na);
                  phi=fluidphi(grid, x, y);
               }
               particles.add_particle(Vec2f(x,y), Vec2f(0,0));
            }
         }
      }
   }
}

void advance_one_step(Grid &grid, Particles &particles, double dt)
{
   for(int i=0; i<5; ++i)
      particles.move_particles_in_grid(0.2*dt);
   particles.transfer_to_grid();
   grid.save_velocities();
   grid.add_gravity(dt);
   grid.compute_distance_to_fluid();
   grid.extend_velocity();
   grid.apply_boundary_conditions();
   grid.make_incompressible();
   grid.extend_velocity();
   grid.get_velocity_update();
   particles.update_from_grid();
}

void advance_one_frame(Grid &grid, Particles &particles, double frametime)
{
   double t=0;
   double dt;
   bool finished=false;
   while(!finished){
      dt=2*grid.CFL();
      if(t+dt>=frametime){
         dt=frametime-t;
         finished=true;
      }else if(t+1.5*dt>=frametime)
         dt=0.5*(frametime-t);
      printf("advancing %g (to %f%% of frame)\n", dt, 100*(t+dt)/frametime);
      advance_one_step(grid, particles, dt);
      t+=dt;
   }
}

int main(int argc, char **argv)
{
   Grid grid(9.8, 50, 50, 1);
   Particles particles(grid);
   char *outputpath=".";

   if(argc>1) outputpath=argv[1];
   else printf("using default output path...\n");
   printf("Output sent to %s\n", outputpath);

   init_water_drop(grid, particles, 2, 2);
   particles.write_to_file("%s/frameparticles%04d", outputpath, 0);

   for(int i=1; i<101; ++i){
      printf("===================================================> step %d...\n", i);
      advance_one_frame(grid, particles, 1./30);
      particles.write_to_file("%s/frameparticles%04d", outputpath, i);
   }

   return 0;
}
