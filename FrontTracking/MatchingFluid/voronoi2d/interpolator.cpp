
// ---------------------------------------------------------
//
//  interpolator.cpp
//  
// ---------------------------------------------------------

#include "interpolator.h"
#include <mat.h>


bool VERBOSE_INTERPOLATION = false;


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// find distance x0 from segment x1-x2

static float point_segment_distance( const Vec2f& x0, const Vec2f& x1, const Vec2f& x2, float& alpha )
{
   Vec2f dx(x2-x1);
   float m2=mag2(dx);

   if ( m2 == 0.0 ) 
   {
      // degenerate edge
      alpha = 1.0;
      return dist( x0, x1 );
   }
   
   // find parameter value of closest point on segment
   alpha = (float)( dot(x2-x0, dx)/m2 );
   if(alpha<0){
      alpha=0;
   }else if(alpha>1){
      alpha=1;
   }
   // and find the distance
   return mag(alpha*x1+(1-alpha)*x2 - x0);
}

// ---------------------------------------------------------------------------

static void get_closest_point_on_polygon( const Vec2f& point, const std::vector<Vec2f>& poly_points, int& index_a, float& alpha )
{
   float min_distance = 1e30f;
   
   for ( unsigned int i = 0; i < poly_points.size(); ++i )
   {
      unsigned int next = (i+1) % poly_points.size();
      float curr_alpha;
      float curr_dist = point_segment_distance( point, poly_points[i], poly_points[next], curr_alpha );
      
      if ( curr_dist < min_distance )
      {
         min_distance = curr_dist;
         index_a = i;
         alpha = curr_alpha;
      }
   }
   
}

// ---------------------------------------------------------------------------

static bool point_is_in_triangle( const Vec2f& point, const Vec2f& a, const Vec2f& b, const Vec2f& c )
{
   Vec2f v0 = c - a;
   Vec2f v1 = b - a;
   Vec2f v2 = point - a;
   
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
   return (u >= -1e-5) && (v >= -1e-5) && (u + v <= 1+2e-5);
   
}


// ---------------------------------------------------------------------------

static bool point_is_in_convex_polygon( const Vec2f& point, const std::vector<Vec2f>& poly_points )
{
   static float EPSILON = 1e-5f;
   
   assert( poly_points.size() > 2 );
   
   float last_sign = cross( poly_points[0] - point, poly_points.back() - point );
   
   for ( unsigned int i = 1; i < poly_points.size(); ++i )
   {
      int prev = ( i - 1 );
      assert( prev >= 0 );
      
      float sign = cross( poly_points[i] - point, poly_points[prev] - point );

      if ( sign >  EPSILON && last_sign < -EPSILON ) { return false; }
      if ( sign < -EPSILON && last_sign >  EPSILON ) { return false; }
         
      last_sign = sign;
   }
   
   return true;
   
}


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Interpolation functions
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------


Vec2f triangulated_polygon_interpolation( const Vec2f& point, const std::vector<Vec2f>& poly_points, const std::vector<Vec2f>& data )
{

   if ( poly_points.size() < 3 )
   {
      //shouldn't really ever hit this.
      return Vec2f( 0, 0 );
   }

   if (!point_is_in_convex_polygon( point, poly_points ) )
   {
      // get nearest point on the polygon
      int index_a = 0;
      float alpha = 0.0f;
      get_closest_point_on_polygon( point, poly_points, index_a, alpha );
      int index_b = (index_a+1) % poly_points.size();      
      alpha = clamp(alpha, 0.0f, 1.0f);
      return alpha*data[index_a] + (1-alpha)*data[index_b];
   }
   
   Vec2f centroid(0.0f, 0.0f);
   Vec2f data_average(0.0f, 0.0f);
   for ( unsigned int i = 0; i < poly_points.size(); ++i )
   {
      centroid += poly_points[i];
      data_average += data[i];
   }
   centroid /= (float) poly_points.size();
   data_average /= (float) poly_points.size();

   float min_distance_to_edge = 1e30f;
   unsigned int min_edge_index = ~0;
   float min_alpha = 1e30f;
   bool outer = false;

   for ( unsigned int i = 0; i < poly_points.size(); ++i )
   {
      unsigned int next = (i+1) % poly_points.size();
      const Vec2f& p1 = poly_points[i];
      const Vec2f& p2 = poly_points[next];
      if ( point_is_in_triangle( point, centroid, p1, p2 ) )
      {
         // Compute twice area of triangle ABC
         
         float AreaABC = cross(p1-centroid,p2-centroid);
         float inv = 1/AreaABC;
         if(AreaABC < 1e-10)
            return 0.5f*(data[i]+data[next]);
         
         // Compute a
         float AreaPBC = cross(p1-point,p2-point);
         float a = AreaPBC * inv;
         
         // Compute b
         float AreaPCA = cross(p2-point,centroid-point);
         float b = AreaPCA * inv;
         
         a = clamp(a, 0.0f,1.0f);
         b = clamp(b, 0.0f,1.0f);
         // Compute c
         float c = clamp(1.0f - a - b, 0.0f, 1.0f);
         
         return a*data_average + b*data[i] + c*data[next];
      }
      
      //check outer edge
      float curr_alpha;
      float curr_dist = point_segment_distance( point, poly_points[i], poly_points[next], curr_alpha );
      if ( curr_dist < min_distance_to_edge )
      {
         outer = true;
         min_distance_to_edge = curr_dist;
         min_alpha = curr_alpha;
         min_edge_index = i;
      }

      //check central edge
      curr_dist = point_segment_distance( point, poly_points[i], centroid, curr_alpha );
      if ( curr_dist < min_distance_to_edge )
      {
         outer = false;
         min_distance_to_edge = curr_dist;
         min_alpha = curr_alpha;
         min_edge_index = i;
      }
      
   }
   
   //std::cout << "warning: point not found inside any triangles, ";
   //std::cout << "min_distance_to_edge: " << min_distance_to_edge << std::endl;
   min_alpha = clamp(min_alpha, 0.0f, 1.0f);
   if(!outer)
      return min_alpha * data[min_edge_index] + (1.0f-min_alpha) * data_average;
   else {
      unsigned int next = (min_edge_index + 1) % poly_points.size();
      return min_alpha * data[min_edge_index] + (1.0f-min_alpha) * data[next];
   }
   
}


//
//// ---------------------------------------------------------------------------


Vec2f reconstruct( const std::vector<Vec2f>& normals, const std::vector<float>& velocities, bool& success )
{
   assert( !normals.empty() );
   success = true;
   if ( normals.size() == 1 )
   {
      success = false;
      return Vec2f(0,0);
   }

   Mat22f ATA;
   Vec2f ATb;
   
   if ( normals.size() == 2 )
   {
      //if the normals are too close to collinear, information in the missing axis is unreliable
      //which can cause arbitrary velocities
      if( fabs(fabs( dot(normals[0],normals[1]) ) - 1.0f) < 0.05 )
      {
         //bit of a hack to ensure we don't overshoot wildly. just ignore the missing dimension.
         //We should instead grab more nearby velocity face samples, to get a better/smoother velocity guess
         success = false;
         return Vec2f(0,0);
      }
      Mat22f A( normals[0][0], normals[1][0], normals[0][1], normals[1][1] );
      ATA = A.transpose() * A;
      ATb = A.transpose() * Vec2f( velocities[0], velocities[1] );
   }
   else
   {
      assert( normals.size() == 3 );
      

      Mat32f A( normals[0][0], normals[1][0], normals[2][0], normals[0][1], normals[1][1], normals[2][1] );
      ATA = A.transpose() * A;
      ATb = A.transpose() * Vec3f( velocities[0], velocities[1], velocities[2] );
      if( fabs(fabs( dot(normals[0],normals[1]) ) - 1.0f) < 5e-2 && fabs(fabs( dot(normals[1],normals[2]) ) - 1.0f) < 5e-2 ) {
         success = false;
         return Vec2f(0,0);
      }
   }
   
   float a = ATA(0,0);
   float b = ATA(0,1);
   float c = ATA(1,0);
   float d = ATA(1,1);
   float e = ATb[0];
   float f = ATb[1];
   
   float denom = a*d - b*c;
   
   if( denom == 0.0 )
   {
      assert(0);
   }
   
   float x = (e*d - b*f) / denom;
   float y = (a*f - e*c) / denom;
   
   Vec2f reconstructed_velocity( x, y );
   
   return reconstructed_velocity;   
}


// ---------------------------------------------------------------------------
   
