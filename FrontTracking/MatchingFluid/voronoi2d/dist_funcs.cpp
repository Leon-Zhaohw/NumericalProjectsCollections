#include "dist_funcs.h"

float circle(const Vec2f& point, const Vec2f& centre, float radius) {
   return dist(point,centre) - radius;
}

static float point_segment_dist(const Vec2f &x0, const Vec2f &x1, const Vec2f &x2)
{
   Vec2f dx(x2-x1);
   float m2=mag2(dx);
   // find parameter value of closest point on segment
   float s=clamp(dot(x2-x0, dx)/m2, 0.0f, 1.0f);
   // and find the distance
   return dist(x0,s*x1+(1-s)*x2);
}

float capsule(const Vec2f& point, const Vec2f& start, const Vec2f& end, float radius) {
   return point_segment_dist(point, start, end) - radius;
}

float box(const Vec2f& point, const Vec2f& centre, float width, float height) {
   Vec2f offset = fabs(point - centre);
   offset -= Vec2f(width/2,height/2);
   if(offset[0] < 0 && offset[1] < 0) {
      return max(offset);
   }
   else {
      float sum = 0;
      if(offset[0] > 0)
         sum += sqr(offset[0]);
      if(offset[1] > 0)
         sum += sqr(offset[1]);
      return sqrt(sum);
   }
}

Vec2f box_gradient( const Vec2f& point, const Vec2f& centre, float width, float height ) 
{

   float ratio = height / width;
   Vec2f rel_point = point - centre;
   
   if ( ratio * fabs( rel_point[0] ) > fabs( rel_point[1] ) )
   {
      if ( rel_point[0] > 0.0 )
      {
         return Vec2f( -1, 0 );
      }
      else
      {
         return Vec2f( 1, 0 );
      }
   }
   else
   {
      if ( rel_point[1] > 0.0 )
      {
         return Vec2f( 0, -1 );
      }
      else
      {
         return Vec2f( 0, 1 );
      }
   }      
   
   
}


float rotated_box(const Vec2f& point, const Vec2f& centre, float width, float height, float angle) {
   
   Vec2f local = rotate(point-centre, angle);
   Vec2f offset = fabs(local);
   offset -= Vec2f(width/2,height/2);
   if(offset[0] < 0 && offset[1] < 0) {
      return max(offset);
   }
   else {
      float sum = 0;
      if(offset[0] > 0)
         sum += sqr(offset[0]);
      if(offset[1] > 0)
         sum += sqr(offset[1]);
      return sqrt(sum);
   }
}

void make_circle(Vec2f centre, float radius, int segments, std::vector<Vec2f>& verts, std::vector<Vec2ui>& edges ) {
   verts.clear();
   edges.clear();
   float d_angle = 2.0f*(float)M_PI / (float)segments;
   float angle = 0;
   for(int i = 0; i < segments; ++i) {
      verts.push_back(centre + radius*Vec2f(cos(angle), sin(angle)));
      edges.push_back(Vec2ui(i,i+1));
      angle+=d_angle;
   }
   edges[edges.size()-1][1] = 0; //fix the last one
}

void make_box(Vec2f bottom_left, float height, float width, std::vector<Vec2f>& verts, std::vector<Vec2ui>& edges) {
   verts.clear();
   edges.clear();

   verts.push_back(bottom_left);
   verts.push_back(bottom_left+Vec2f(width,0));
   verts.push_back(bottom_left+Vec2f(width,height));
   verts.push_back(bottom_left+Vec2f(0,height));

   edges.push_back(Vec2ui(0,1));
   edges.push_back(Vec2ui(1,2));
   edges.push_back(Vec2ui(2,3));
   edges.push_back(Vec2ui(3,0));
}
