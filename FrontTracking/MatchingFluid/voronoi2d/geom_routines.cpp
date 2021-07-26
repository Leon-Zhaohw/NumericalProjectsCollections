#include "geom_routines.h"

bool point_in_box(Vec2f& pt, Vec2f& boxMin, float width, float height) {
   return pt[0] >= boxMin[0] && pt[1] >= boxMin[1] && pt[0] <= boxMin[0]+width && pt[1] <= boxMin[1]+height;
}

bool point_in_tri(const Vec2f& point, const Vec2f& p0, const Vec2f& p1, const Vec2f& p2) {
   Vec2f v0 = p2 - p0;
   Vec2f v1 = p1 - p0;
   Vec2f v2 = point - p0;

   // Compute dot products
   double dot00 = dot(v0, v0);
   double dot01 = dot(v0, v1);
   double dot02 = dot(v0, v2);
   double dot11 = dot(v1, v1);
   double dot12 = dot(v1, v2);

   // Compute barycentric coordinates
   double invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
   double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
   double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

   // Check if point is in triangle (TODO Make this more consistent)
   return (u >= -1e-6) && (v >= -1e-6) && (u + v <= 1+1e-6);

}

bool seg_hits_box(const Vec2f&a, const Vec2f&b, const Vec2f&x00, float width, float height) {
   //check top
   float top = x00[1]+height;
   if(a[1] < top && b[1] > top || a[1] > top && b[1] < top) {
      float t = (top - a[1]) / (b[1] - a[1]);
      float xpos = (1-t)*a[0] + t*b[0];
      if(xpos > x00[0] && xpos < x00[0] + width)
         return true;
   }

   //check bottom
   float bottom = x00[1];
   if(a[1] < bottom && b[1] > bottom || a[1] > bottom && b[1] < bottom) {
      float t = (bottom - a[1]) / (b[1] - a[1]);
      float xpos = (1-t)*a[0] + t*b[0];
      if(xpos > x00[0] && xpos < x00[0] + width)
         return true;
   }

   //check right
   float right = x00[0]+width;
   if(a[0] < right && b[0] > right || a[0] > right && b[0] < right) {
      float t = (right - a[0]) / (b[0] - a[0]);
      float ypos = (1-t)*a[1] + t*b[1];
      if(ypos > x00[1] && ypos < x00[1] + height)
         return true;
   }

   //check left
   float left = x00[0];
   if(a[0] < left && b[0] > left || a[0] > left && b[0] < left) {
      float t = (left - a[0]) / (b[0] - a[0]);
      float ypos = (1-t)*a[1] + t*b[1];
      if(ypos > x00[1] && ypos < x00[1] + height)
         return true;
   }
   
   return false;

}

bool tri_intersects_square(Vec2f& a, Vec2f& b, Vec2f& c, Vec2f& x00, float width, float height) 
{
   //There are more efficient ways of doing this.
   return  (point_in_box(a, x00, width, height) ||  //check tri points in the box
            point_in_box(b, x00, width, height) || 
            point_in_box(c, x00, width, height) ||
            point_in_tri(x00, a,b,c) ||            //check box points in the tri
            point_in_tri(x00+Vec2f(width,0), a,b,c) || 
            point_in_tri(x00+Vec2f(width,height), a,b,c) || 
            point_in_tri(x00+Vec2f(0,height), a,b,c) ||
            seg_hits_box(a, b, x00, width, height) ||                      //check tri segs hit the box
            seg_hits_box(a, c, x00, width, height) ||
            seg_hits_box(b, c, x00, width, height));
}

void get_bound_box(const Vec2f& a, const Vec2f& b, const Vec2f& c,  Vec2f& minPt, Vec2f& maxPt) {
   minPt[0] = min(a[0], min(b[0], c[0]));
   minPt[1] = min(a[1], min(b[1], c[1]));
   maxPt[0] = max(a[0], max(b[0], c[0]));
   maxPt[1] = max(a[1], max(b[1], c[1]));
}