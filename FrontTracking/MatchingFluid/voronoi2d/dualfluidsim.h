#ifndef DUALFLUIDSIM_H
#define DUALFLUIDSIM_H

#include <vector>
#include "trimesh2d.h"

#include "surftrack.h"

typedef float (*dist_func)(const Vec2f &);

class VelocityFunctor {
public:
   
   virtual ~VelocityFunctor() {}
   
   virtual Vec2f operator()(const Vec2f& point) const {
      return Vec2f(0,0);
   };
};

class DualFluidSim {
public:

   eltopo2d::SurfTrack& explicit_surface;
   
   TriMesh2D mesh;
   float characteristic_distance;
   
   Vec2f domain_min, domain_max;
   Vec2f solid_bounds_centre, solid_bounds_dims;

   bool conserve_volume;

   //cache of clipped polygons to speed things up a bit.
   std::vector< std::vector<Vec2f> > per_vert_poly_verts;
   std::vector< std::vector<Vec2f> > per_vert_poly_vels;

   //Information for boundary conditions
   std::vector<float> liquid_phi; //located at voronoi sites (for GFM)
   std::vector<float> solid_phi;  //located at voronoi vertices (for FVM)
   std::vector<float> solid_weights; //FVM weights on faces

   //Velocities at various voronoi mesh locations
   std::vector<float> face_velocities; 
   std::vector<Vec2f> elt_velocities;
   
   std::vector<Vec2f> voronoi_vertex_velocities;
   std::vector<bool> voronoi_vertex_is_valid;
   std::vector<bool> voronoi_edge_is_valid;
   
   std::vector<double> pressures;
   
   std::vector<char> valid_elts;

   std::vector<float> temp_face;

   std::vector<Vec2f> markers;
   
   float target_volume;
   
   DualFluidSim( eltopo2d::SurfTrack& explicit_surface_ ) :
      explicit_surface( explicit_surface_ ), conserve_volume(false)
   {}
   
private:
   DualFluidSim();
   
public:
   
   void setup( TriMesh2D &mesh, float approx_dx, const Vec2f& solid_bounds_centre, const Vec2f& solid_bounds_dims,
                                                 const Vec2f& domain_min_, const Vec2f& domain_max_);
   void change_mesh(const TriMesh2D& newmesh, const std::vector<float>& newFaceVels);

  
   float solid_func(const Vec2f& pos);
   Vec2f solid_func_gradient(const Vec2f& pos);
   
   void compute_solids();
   
   float get_cfl_step();

   void compute_liquidphi(std::vector<int>& expected_signs);
   float interpolated_liquid_phi( const Vec2f& pt );
   void extrapolate_liquid_phi_into_solid();
   
   void advance_surface(float timestep);
   void advance( float timestep );

   void add_forces(float timestep);

   void volume_conservation( std::vector<Vec2d>& vertex_velocities, float dt );
   
   bool trace_rk2(const Vec2f& start, Vec2f& end, float dt, const VelocityFunctor& get_velocity );
   bool trace_rk3(const Vec2f& start, Vec2f& end, float dt, const VelocityFunctor& get_velocity );

   void advect_semilagrangian( float dt, 
                               const std::vector<Vec2f>& starting_points, 
                               std::vector<float>& new_velocities, 
                               const std::vector<Vec2f>& new_edge_normals, 
                               const VelocityFunctor& get_velocity );
   
   void solve_pressure();
               
   
   // ----------
   //
   // velocity functions
   //
   // ----------
   
 
   // get velocity using clipped polygon interpolation, reconstructing velocities 
   // at the polygon corners
   //
   Vec2f get_velocity_voronoi(const Vec2f& point );
   
   // get velocity using clipped polygon interpolation, assuming velocities 
   // have already been reconstructed on the voronoi vertices
   //
   Vec2f get_velocity_voronoi_consistent( const Vec2f& point );   
   

   // ----------
   //
   // helpers for polygon interpolation
   //
   // ----------
   
   bool triangles_are_adjacent_ordered( unsigned int tri1_index, 
                                        unsigned int tri2_index, 
                                        const unsigned int central_vertex, 
                                        int& shared_edge_index );
   
   // Clip a polygon against the solid, adding in zero normal velocities for new (solid) edges
   //   
   void clip_polygon_to_solid( const std::vector<int>& polygon_vertex_indices, 
                               const std::vector<float>& polygon_edge_velocities,
                               std::vector<Vec2f>& new_polygon_vertex_locations,
                               std::vector<float>& new_polygon_edge_velocities );
   
   // Clip a polygon against the solid, using existing vertex velocities if possible.
   // For new vertices (created by clipper), reconstruct using one existing edge velocity
   // and the u.n = 0 condition for the solid edge.
   //   
   void clip_polygon_to_solid_consistent( const std::vector<int>& polygon_vertex_indices, 
                                          const std::vector<Vec2f>& polygon_vertex_velocities, 
                                          const std::vector<float>& polygon_edge_velocities,
                                          std::vector<Vec2f>& new_polygon_vertex_locations,
                                          std::vector<Vec2f>& new_polygon_vertex_velocities );
   
   // Get the set of vertices surrounding the given voronoi site (in order)
   //
   void construct_voronoi_polygon( unsigned int voronoi_site_index, 
                                   std::vector<int>& polygon_vertices,
                                   std::vector<float>& polygon_edge_velocities );
   
   // For a given set of polygon vertices, and normal velocity components on the polygon edges, 
   // reconstruct 2D velocity vectors at each polygon vertex
   //
   void reconstruct_velocities_for_polygon( const std::vector<Vec2f>& polygon_vertex_locations, 
                                            const std::vector<float>& polygon_edge_velocities,
                                            std::vector<Vec2f>& polygon_vertex_velocities ); 
      
   // ----------
   //
   // extrapolation/velocity transfer functions
   //
   // ----------
   
 
   // reconstruct velocities on the voronoi vertices using the "consistent", triple-edge reconstruction
   // then extrapolate these velocities out
   //
   void extrapolate_consistent_vertex_velocities();
   
   // reconstruct velocities on the voronoi vertices using the "consistent", triple-edge reconstruction
   // then extrapolate these velocities out.  Use these extrapolated vertex velocities to get normal velocity
   // components on edges.
   //
   void extrapolate_voronoi_edge_velocities();

};


//
//// assuming we have 2d velocity vectors at all triangle vertices, use barycentric interpolation
////
//class DualFluidVelocityFunctor : public VelocityFunctor {
//   DualFluidSim& sim;
//public:
//   DualFluidVelocityFunctor(DualFluidSim& sim_):sim(sim_){}
//
//   Vec2f operator()(const Vec2f& pt) const {
//      int tri = sim.mesh.get_containing_tri(pt);
//      if(tri == -1)
//         return Vec2f(0,0);
//      else {
//         return sim.get_velocity_triangle(pt, sim.elt_velocities, tri);
//      }
//   }
//};


//
//// On a per-voronoi-region basis:
//// First clip the voronoi polygon against the solid, setting the normal velocity of new, solid edges to zero.
//// Then use the normal velocity samples on the polygon edges to reconstruct 2d velocity vectors at voronoi vertices.
//// Then use generalized barycentric interpolation to get velocity at a point in space.
////
//class DualFluidFancyVelocityFunctor : public VelocityFunctor {
//   DualFluidSim& sim;
//public:
//   DualFluidFancyVelocityFunctor(DualFluidSim& sim_):sim(sim_){}
//   
//   Vec2f operator()(const Vec2f& pt) const 
//   {
//      return sim.get_velocity_voronoi(pt);
//   }   
//};


// Assuming we have 2d velocity vectors on all voronoi vertices,
// First clip the voronoi polygon against the solid, constructing the velocity at the new, solid 
// points from the existing sample and the implied zero normal velocity on the solid boundary.
// Then use generalized barycentric interpolation to get velocity at a point in space.
//
class PolygonVelocityFunctor : public VelocityFunctor {
   DualFluidSim& sim;
public:
   
   virtual ~PolygonVelocityFunctor() {}
   
   PolygonVelocityFunctor(DualFluidSim& sim_):sim(sim_){}
   
   Vec2f operator()(const Vec2f& pt) const 
   {
      return sim.get_velocity_voronoi_consistent(pt);
   }   
};





#endif