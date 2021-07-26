#include "scenes.h"

#include "dist_funcs.h"

extern Vec2f solid_bounds_centre;
extern Vec2f solid_bounds_dims;

float sinusoid(const Vec2f& pt) {
   float sineval = pt[1] - 0.5f + 0.04f*(float)sin(7*pt[0]);
   //float sineval = pt[1] + 0.1*pt[0] - 0.5;
   float boxval = box(pt, Vec2f(0.5f,0.5f), 0.85f,0.72f);
   return max(boxval,sineval);
}

float stillpool(const Vec2f& pt) {
   float flatsurface = -0.35f + pt[1];
   float boxcurve = box(pt, solid_bounds_centre, solid_bounds_dims[0], solid_bounds_dims[1]);

   return max(flatsurface, boxcurve);
}

float sinewave(const Vec2f& pt) {
   float sinecurve = -0.35f + pt[1] + 0.02f*(float)sin(30*pt[0]);
   float boxcurve = box(pt, solid_bounds_centre, solid_bounds_dims[0], solid_bounds_dims[1]);

   return max(sinecurve, boxcurve);
}

float stillpool_and_circle(const Vec2f& pt) {

   return min(stillpool(pt), circle(pt, solid_bounds_centre + Vec2f(0.05f,0.06f), 0.1f));
   //return min(stillpool(pt), circle(pt, solid_bounds_centre + Vec2f(0.05f,0.065f), 0.1f));
}

float sinewaveandcircle(const Vec2f& pt) {
   
   return min(sinewave(pt), circle(pt, solid_bounds_centre, 0.1f));
}


float bigcircle(const Vec2f& pt) 
{
   static const Vec2f CIRCLE_CENTRE( solid_bounds_centre - Vec2f( 0.0f, 0.25f*solid_bounds_dims[1] ) );
   static const float CIRCLE_RADIUS = 0.11f;

   //Vec2f centre(0.5, 0.4);
   //float circlePhi = circle(pt,centre,rad);   
   float circlePhi = circle( pt, CIRCLE_CENTRE, CIRCLE_RADIUS );
   
   return circlePhi;
}

float twocircles(const Vec2f& pt) 
{
   return min(circle(pt, Vec2f(0.5f, 0.52f), 0.07f), circle(pt, Vec2f(0.5f, 0.42f), 0.07f));
}

float twocirclesbeside(const Vec2f& pt) 
{
   return min(circle(pt, Vec2f(0.35f, 0.31f), 0.08f), circle(pt, Vec2f(0.69f, 0.31f), 0.1f));
}
