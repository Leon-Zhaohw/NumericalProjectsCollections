#ifndef RIGIDSHAPE2_H
#define RIGIDSHAPE2_H

#include <cmath>
#include <cstdarg>
#include "vec.h"

struct RigidShape2
{
   Vec2f centre;
   float rotation;
   float radius; // positive if we include a disc in the shape
   std::vector<Vec2f> x;
   std::vector<Vec2ui> edge;
   Vec2f velocity;
   float angular_velocity;
   // rigid motion with respect to the origin, not centre
   Vec2f global_velocity;
   float global_angular_velocity;

   float c, s; // cos and sin of rotation

   RigidShape2(void)
      : centre(0,0), rotation(0), radius(-1), velocity(0,0), angular_velocity(0), global_velocity(0,0), global_angular_velocity(0), c(1), s(0)
   {}

   void set_disc_radius(float radius_)
   {
      radius=radius_;
   }

   void add_strip(int num_points, ...)
   {
      unsigned int s=x.size();
      va_list ap;
      va_start(ap, num_points);
      for(int i=0; i<num_points; ++i){
         double px=va_arg(ap, double);
         double py=va_arg(ap, double);
         x.push_back(Vec2f((float)px,(float)py));
      }
      va_end(ap);
      for(int i=1; i<num_points; ++i){
         edge.push_back(Vec2ui(s+i-1,s+i));
      }
   }

   void set_velocity(const Vec2f &new_velocity, float new_angular_velocity)
   {
      velocity=new_velocity;
      angular_velocity=new_angular_velocity;
      global_velocity=velocity-angular_velocity*perp(centre);
      global_angular_velocity=angular_velocity;
   }

   void set_pose(const Vec2f &new_centre, float new_rotation)
   {
      centre=new_centre;
      rotation=new_rotation;
      c=std::cos(rotation);
      s=std::sin(rotation);
      global_velocity=velocity-angular_velocity*perp(centre);
      global_angular_velocity=angular_velocity;
   }

   void update_pose(float dt)
   {
      set_pose(centre+dt*velocity, rotation+dt*angular_velocity);
   }

   Vec2f world_x(unsigned int i) const
   {
      return centre+Vec2f(c*x[i][0]+s*x[i][1], -s*x[i][0]+c*x[i][1]);
   }

   float distance(const Vec2f &p) const
   {
      Vec2f q(c*(p[0]-centre[0])-s*(p[1]-centre[1]), s*(p[0]-centre[0])+c*(p[1]-centre[1]));
      float d=1e37;
      for(unsigned int e=0; e<edge.size(); ++e){
         unsigned int i, j; assign(edge[e], i, j);
         Vec2f dx(x[j]-x[i]);
         float m2=mag2(dx);
         float s=dot(x[j]-q, dx)/m2;
         if(s>1) s=1; else if(s<0) s=0;
         float thisd=mag(s*x[i]+(1-s)*x[j]-q);
         if(thisd<d) d=thisd;
      }
      if(radius>0){
         float thisd=mag(q)-radius;
         if(thisd<d) d=thisd;
      }
      return d;
   }
};

#endif
