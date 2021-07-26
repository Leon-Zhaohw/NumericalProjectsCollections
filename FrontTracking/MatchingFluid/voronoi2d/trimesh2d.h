#ifndef TRIMESH2D
#define TRIMESH2D

#include <vector>
#include "vec.h"
#include "geom_routines.h"

class TriMesh2D {

public:

   //Geometry
   std::vector<Vec2f> vertices;           //voronoi sites
   std::vector<Vec2f> tri_circumcentres;  //voronoi vertices

   //Neighbour relationships
   std::vector< std::vector<int> > vert_to_edge_map;
   std::vector<Vec2ui> edges; //edge_to_vert_map
   
   std::vector<Vec2i> edge_to_tri_map;    // == voronoi edges
   std::vector<Vec3ui> tri_to_edge_map;   // == voronoi vertex to voronoi edge map

   std::vector<Vec3ui> tris;  //tri_to_vert_map
   std::vector< std::vector<int> > vert_to_tri_map;
   
   //Other data
   std::vector<Vec3i> tri_edge_signs;
   std::vector<Vec2f> edge_normals;
   std::vector<float> edge_lengths;
   std::vector<Vec2f> edge_midpoints;
   
   std::vector< std::vector<int> > dual_edge_signs; //should be consistent with the tri_edge_signs
   std::vector<Vec2f> dual_edge_normals;
   std::vector<Vec2f> dual_edge_midpoints;
   std::vector<float> dual_edge_lengths; 
   
   std::vector<Vec2f> tri_barycentres;
   std::vector<float> tri_areas;
   std::vector<float> voronoi_areas;

   //Optional extra checks on mesh validity
   std::vector<bool> tri_delaunay;

   //Acceleration grid information
   std::vector< std::vector<int> > accel_grid;
   float accel_dx;
   float inv_accel_dx;
   Vec2f accel_origin;
   int accel_ni, accel_nj;
   
   std::vector< std::vector<int> > vert_accel_grid;
   
   TriMesh2D();
   void setup_mesh(const std::vector<Vec2f>& vertices, const std::vector<Vec3ui>& tris);
   
   int get_containing_tri(const Vec2f& point) const;
   
   void get_overlapping_cells( const Vec2f& centre, float radius, std::vector<int>& cells ) const;
   int get_nearest_vertex(const Vec2f& point) const;

   bool point_in_tri(const Vec2f& point, int tri_id)  const;
   
   void get_barycentric_weights_tri(int containing_tri, const Vec2f& pt, Vec3f& weights) const;
   void get_barycentric_weights_voronoi(int closest_vertex, const Vec2f& pt, std::vector<float>& weights) const;

   int get_edge_index( const Vec2ui& edge_query ) const;
   
   int get_shared_edge_index( unsigned int tri_a, unsigned int tri_b ) const;
   
private:
   void compute_edges();
   void compute_tris();
   void compute_nbrs();
   void compute_accelerator();
   

};

Vec2f circumcentre(Vec2f v0, Vec2f v1, Vec2f v2);
Vec2f barycentre(Vec2f v0, Vec2f v1, Vec2f v2);


#endif
