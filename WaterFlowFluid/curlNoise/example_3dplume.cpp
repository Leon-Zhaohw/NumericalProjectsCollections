#include "example_3dplume.h"

Example3DPlume::
Example3DPlume(void)
   : sphere_centre(0.3, 2, 0),
     sphere_radius(0.5),
     plume_height(8),
     noise_lengthscale(3),
     noise_gain(3),
     ring_radius(0.6),
     ring_speed(0.3),
     rings_per_second(1),
     ring_magnitude(7),
     ring_falloff(0.7),
     seed(0),
     particles_per_second(200),
     seed_radius(0.5),
     initial_band(0.1)
{
   noise_lengthscale[0]=0.4;
   noise_gain[0]=1;
   noise_lengthscale[1]=0.23;
   noise_gain[0]=0.5;
   noise_lengthscale[2]=0.11;
   noise_gain[0]=0.25;
}

float Example3DPlume::
distance_and_normal(float x, float y, float z, Vec3f &normal) const
{
   float phi=y;
   Vec3f u(x-sphere_centre[0], y-sphere_centre[1], z-sphere_centre[2]);
   float d=mag(u);
   if(d-sphere_radius<phi){
      phi=d-sphere_radius;
      normal=u/(d+1e-10);
   }else{
      normal=Vec3f(0,1,0);
   }
   return phi;
}

void Example3DPlume::
match_boundary(Vec3f &psi, float d, float lengthscale, const Vec3f &normal) const
{
   float alpha=ramp(std::fabs(d)/lengthscale);
   psi=alpha*psi+((1-alpha)*dot(psi,normal))*normal;
}

bool Example3DPlume::
seed_particles(std::vector<Vec3f> &x, float dt) const
{
   unsigned int num_new=(unsigned int)(dt*particles_per_second);
   x.reserve(x.size()+num_new);
   for(unsigned int i=0; i<num_new; ++i){
      float theta=randhashf(seed++, 0, 2*M_PI);
      float r=randhashf(seed++, 0, seed_radius);
      float y=randhashf(seed++, 0, initial_band);
      x.push_back(Vec3f(r*std::cos(theta), y, r*std::sin(theta)));
   }
   return true;
}

Vec3f Example3DPlume::
potential(float x, float y, float z) const
{
   Vec3f psi(0,0,0);
   Vec3f normal;
   float d=distance_and_normal(x, y, z, normal);
   // add turbulence octaves that respect boundaries, increasing upwards
   float height_factor=ramp(y/plume_height);
   for(unsigned int i=0; i<noise_lengthscale.size(); ++i){
      float sx=x/noise_lengthscale[i];
      float sy=y/noise_lengthscale[i];
      float sz=z/noise_lengthscale[i];
      Vec3f psi_i(noise0(sx,sy,sz), noise1(sx,sy,sz), noise2(sx,sy,sz));
      match_boundary(psi_i, d, noise_lengthscale[i], normal);
      psi+=height_factor*noise_gain[i]*psi_i;
   }
   // add rising vortex rings
   float ring_y=ring_speed*t;
   while(ring_y>0){
      float ry=y-ring_y;
      float rr=std::sqrt(x*x+z*z);
      float rmag=ring_magnitude/(sqr(rr-ring_radius)+sqr(rr+ring_radius)+sqr(ry)+ring_falloff);
      Vec3f rpsi=rmag*Vec3f(z, 0, -x);
      match_boundary(rpsi, d, ring_radius, normal);
      psi+=rpsi;
      ring_y-=ring_speed/rings_per_second;
   }

   return psi;
}

void Example3DPlume::
advance_time(float dt)
{
   t+=dt;
   noise.set_time(0.5f*noise_gain[0]/noise_lengthscale[0]*t);
}

