
// ---------------------------------------------------------
//
//  interpolator.h
//
// ---------------------------------------------------------

#include <array2.h>
#include <vec.h>


// triangulate the polygon using its centroid, then use barycentric interpolation to get values from the polygon points
//
Vec2f triangulated_polygon_interpolation( const Vec2f& point, 
                                          const std::vector<Vec2f>& poly_points, 
                                          const std::vector<Vec2f>& data );

// use generalized barycentric coordinates as in Meyer et al. to get values from the polygon points
//
Vec2f generalized_barycentric_interpolation( const Vec2f& point, 
                                             const std::vector<Vec2f>& poly_points, 
                                             const std::vector<Vec2f>& data );

//same thing, but don't filter for short edges (assumes the data passed in has no short edges)
Vec2f generalized_barycentric_interpolation_no_filter( const Vec2f& point, const std::vector<Vec2f>& poly_points, const std::vector<Vec2f>& data );

// given a set of incident edge normals and their scalar normal velocities, reconstruct the 2d velocity at a point
//
Vec2f reconstruct( const std::vector<Vec2f>& normals, 
                   const std::vector<float>& velocities,
                   bool& success );

// given two incident edge normals and their scalar normal velocities, reconstruct the 2d velocity at a point
//
Vec2f reconstruct( const Vec2f& n1, 
                   float sample1, 
                   const Vec2f& n2, 
                   float sample2);



