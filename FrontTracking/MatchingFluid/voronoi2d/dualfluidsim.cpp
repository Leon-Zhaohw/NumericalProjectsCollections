#include <cfloat>
#include <queue>
#include <iostream>

#include "dualfluidsim.h"
#include "pressure_vor2d.h"
#include "dynamicsurface.h"
#include "sampleseeder.h"
#include "interpolator.h"
#include "surftrack.h"

#include "triangle_wrapper.h"
#include "wallclocktime.h"

#include "dist_funcs.h"
#include "makelevelset2.h"
#include "meshlevelset.h"


// -------------------------------------------------------------

// -------------------------------------------------------------

void DualFluidSim::setup(TriMesh2D& m, float approx_dx, const Vec2f& solid_bounds_centre_, const Vec2f& solid_bounds_dims_,
                         const Vec2f& domain_min_, const Vec2f& domain_max_) {

   characteristic_distance = approx_dx;
   solid_bounds_centre = solid_bounds_centre_;
   solid_bounds_dims = solid_bounds_dims_;
   domain_min = domain_min_;
   domain_max = domain_max_;
   
   mesh = m;

   elt_velocities.resize(mesh.vertices.size(), Vec2f(0,0));
   face_velocities.resize(mesh.edges.size(), 0);

   voronoi_vertex_velocities.resize( mesh.tris.size(), Vec2f(0,0) );
   
   liquid_phi.resize(mesh.vertices.size());
   solid_weights.resize(mesh.edges.size());
   solid_phi.resize(mesh.tris.size());
   
   compute_solids();
   
   target_volume = (float)explicit_surface.get_area();
   
   std::cout << "Initial volume " << target_volume << std::endl;
   
}

// -------------------------------------------------------------

void DualFluidSim::change_mesh(const TriMesh2D& newmesh, const std::vector<float>& newFaceVels) {
   //transfer the mesh and face velocities
   mesh = newmesh;
   face_velocities = newFaceVels;

   //resize the other dependent data
   elt_velocities.resize(mesh.vertices.size());
   
   liquid_phi.resize(mesh.vertices.size());
   solid_phi.resize(mesh.tris.size());
   solid_weights.resize(mesh.edges.size());
   
   pressures.resize(mesh.vertices.size());
   valid_elts.resize(mesh.vertices.size());

   temp_face.resize(mesh.edges.size());
}

// -------------------------------------------------------------

float DualFluidSim::solid_func(const Vec2f& pos) 
{
   float surround = box( pos, solid_bounds_centre, solid_bounds_dims[0], solid_bounds_dims[1]);
   return surround;
}

// -------------------------------------------------------------

Vec2f DualFluidSim::solid_func_gradient(const Vec2f& pos) 
{
   return box_gradient( pos, solid_bounds_centre, solid_bounds_dims[0], solid_bounds_dims[1] );
}


// -------------------------------------------------------------


//On a signed distance field, compute the fraction inside the (negative) isosurface
//along the 1D line segment joining phi_left and phi_right sample points.
float fraction_inside(float phi_left, float phi_right)
{
   //Note: should not generate divide by zero, because if
   //signs are different, and both values are != 0,
   //abs(left - right) is necessarily >= left, or right, alone, ie. not 0
   
   return 
      (phi_left >= 0 && phi_right >= 0)? //all empty
         0.0f : 
         (  (phi_left < 0 && phi_right < 0)? //all full
            1.0f:
            (
               (phi_left >= 0)?
               (1 - phi_left / (phi_left - phi_right)): //right inside 
               (phi_left / (phi_left - phi_right)) //left inside
            )
         );
}


void DualFluidSim::compute_solids() {
   
   
   //solid phi for a start, just a box
   for(unsigned int i = 0; i < mesh.tris.size(); ++i) {
      solid_phi[i] = solid_func(mesh.tri_circumcentres[i]);
   }
   
   //compute solid weights via edge fractions
   for(unsigned int i = 0; i < mesh.edges.size(); ++i) {
      Vec2i edge_tri_ids = mesh.edge_to_tri_map[i];
      
      //if one of the triangles doesn't exist, assume this face is closed.
      if(edge_tri_ids[0] == -1 || edge_tri_ids[1] == -1) {
         solid_weights[i] = 0;
         continue;
      }

      float phi0 = solid_phi[edge_tri_ids[0]];
      float phi1 = solid_phi[edge_tri_ids[1]];
      solid_weights[i] = fraction_inside(phi0,phi1);
      if(solid_weights[i] < 0.02) //not sure this is actually necessary
         solid_weights[i] = 0;

   }
   
}


// -------------------------------------------------------------


float DualFluidSim::get_cfl_step() {
   
   float max_vel = 0;
   for(unsigned int i = 0; i < face_velocities.size(); ++i)
      max_vel = max(max_vel, face_velocities[i]);

   std::cout << "Max velocity: " << max_vel << std::endl;
   
   return 0.5f*characteristic_distance / max_vel;
}


// -------------------------------------------------------------

void DualFluidSim::compute_liquidphi(std::vector<int>& expected_signs) 
{

   liquid_phi.clear();
   liquid_phi.resize( mesh.vertices.size(), ~0 );
   

   std::vector<Vec2d> double_vertices;
   double_vertices.reserve( mesh.vertices.size() );
   for ( unsigned int i = 0; i < mesh.vertices.size(); ++i )
   {
      double_vertices.push_back( Vec2d( mesh.vertices[i] ) );
   }

   std::vector<double> double_phi;
   make_mesh_level_set( explicit_surface,
                        mesh.tris,
                        double_vertices,
                        double_phi );

   for ( unsigned int i = 0; i < double_phi.size(); ++i )
   {
      liquid_phi[i] = (float) double_phi[i];
      if(expected_signs[i]!=0) {
         liquid_phi[i] = expected_signs[i] * fabs(liquid_phi[i]);
      }

   }
         
 
}

// -------------------------------------------------------------

float DualFluidSim::interpolated_liquid_phi( const Vec2f& pt )
{
   int tri_index = mesh.get_containing_tri( pt );
   if(tri_index == -1) 
      return 4*characteristic_distance;

   Vec3f weights;
   mesh.get_barycentric_weights_tri( tri_index, pt, weights );
   Vec3ui& tri = mesh.tris[tri_index];
   float interpolated_value =   weights[0] * liquid_phi[tri[0]]
                              + weights[1] * liquid_phi[tri[1]]
                              + weights[2] * liquid_phi[tri[2]];
   
   return interpolated_value;   
   
}

// -------------------------------------------------------------

void DualFluidSim::extrapolate_liquid_phi_into_solid() 
{   
   
   for ( unsigned int pass = 0; pass < 5; ++pass )
   {
      std::vector<float> new_liquid_phi = liquid_phi;

      for( unsigned int i = 0; i < mesh.vertices.size(); ++i ) 
      {
         const Vec2f& phi_location = mesh.vertices[i];
         float solid_phi = solid_func( phi_location );
         if ( solid_phi > 0.0 )
         {
            // voronoi cell is solid
            float sample_offset = solid_phi + 0.25f * characteristic_distance;
            Vec2f non_solid_sample = phi_location + sample_offset * solid_func_gradient( phi_location );

                        
            int tri_index = mesh.get_containing_tri( non_solid_sample );
            if ( tri_index >= 0 )
            {
               new_liquid_phi[i] = interpolated_liquid_phi( non_solid_sample ) - 0.1f * characteristic_distance;         
            }
         }
         
      }
      
      liquid_phi = new_liquid_phi;
   }   
   
}


// -------------------------------------------------------------



void DualFluidSim::advance_surface(float timestep) 
{

   PolygonVelocityFunctor get_velocity( *this );

   explicit_surface.improve_mesh();
   
   explicit_surface.topology_changes();
   
   float total_mag = 0.0;
   float max_mag = -1.0;
   int max_mag_vertex;
   
   std::vector<Vec2f> vertex_velocities( explicit_surface.m_positions.size() );
   for ( unsigned int i = 0; i < explicit_surface.m_positions.size(); ++i )
   {
      Vec2f x( explicit_surface.m_positions[i] );
      vertex_velocities[i] = get_velocity( x );
      
      if ( solid_func( x ) > -0.5f * characteristic_distance )
      {
         explicit_surface.m_masses[i] = 2.0;
      } 
      else
      {
         explicit_surface.m_masses[i] = 1.0;
      }
      
      total_mag += mag(vertex_velocities[i]);
      if ( mag(vertex_velocities[i]) > max_mag )
      {
         max_mag = mag(vertex_velocities[i]);
         max_mag_vertex = i;
      }
   }
   
   
   
   
   
   for ( unsigned int i = 0; i < explicit_surface.m_positions.size(); ++i )
   {   
      Vec2f end;
      trace_rk3( Vec2f(explicit_surface.m_positions[i]), end, timestep, get_velocity );
      //snap points back to the boundary   
      if ( solid_func( end ) > 0.0 )
      {
         end += 0.9f*solid_func( end ) * solid_func_gradient( end );
      }
      explicit_surface.m_newpositions[i] = Vec2d(end);
   }

   if(conserve_volume)
      volume_conservation( explicit_surface.m_newpositions, timestep );

   explicit_surface.integrate( timestep );
   
}

// ------------------------------------------------------------

void DualFluidSim::advance( float timestep ) 
{
   
   //Cache of per vertex interpolation polygons and data
   per_vert_poly_verts.clear();
   per_vert_poly_verts.resize(mesh.vertices.size());
   per_vert_poly_vels.clear();
   per_vert_poly_vels.resize(mesh.vertices.size());
   
   //set up interpolator
   PolygonVelocityFunctor get_velocity(*this);  
     
   std::cout << "Extrapolate velocities" << std::endl;
   extrapolate_consistent_vertex_velocities();
   

   std::cout << "Advecting surface." << std::endl;
   advance_surface( timestep );
  
   std::cout << "Remeshing domain. ";
   std::vector<int> expected_signs;
   std::vector<Vec2f> newPressurePoints;
   std::vector<Vec3ui> newTris;
   double dx = 2*characteristic_distance;
   
   // add regular points
   SampleSeeder::generate_regular_samples( Vec2d(domain_min), Vec2d(domain_max), dx, newPressurePoints, expected_signs );

   for ( unsigned int i = 0; i < newPressurePoints.size(); ++i )
   {
      for ( unsigned int j = i+1; j < newPressurePoints.size(); ++j )
      {
         assert( dist( newPressurePoints[i], newPressurePoints[j] ) > 1e-6 );
      }
   }

   // add surface-adaptive points
   explicit_surface.rebuild_static_broad_phase();
   SampleSeeder::generate_samples( explicit_surface, 0.5 * dx, newPressurePoints, expected_signs);

   for ( unsigned int i = 0; i < newPressurePoints.size(); ++i )
   {
      for ( unsigned int j = i+1; j < newPressurePoints.size(); ++j )
      {
         assert( dist( newPressurePoints[i], newPressurePoints[j] ) > 1e-6 );
      }
   }

   TriMesh2D newmesh;
   compute_delaunay(newPressurePoints, newTris);
   newmesh.setup_mesh(newPressurePoints, newTris);
   std::cout << "New tri count: " << newTris.size() << std::endl;
   
   std::cout << "Apply semi-Lagrangian advection." << std::endl;
   std::vector<float> new_face_velocities(newmesh.edges.size(), 0);
   advect_semilagrangian( timestep, newmesh.dual_edge_midpoints, new_face_velocities, newmesh.dual_edge_normals, get_velocity ); 

   change_mesh(newmesh, new_face_velocities);
   compute_solids();
 
   std::cout << "Adding gravity." << std::endl;
   add_forces(timestep);

   std::cout << "Compute signed distance field." << std::endl;;
   compute_liquidphi(expected_signs);
   extrapolate_liquid_phi_into_solid();

   std::cout << "Solve pressure projection." << std::endl;   
   solve_pressure();

}


// -------------------------------------------------------------


void DualFluidSim::add_forces(float dt) 
{
   
   //Apply some gravity
   Vec2f gravity(0.0f,-1.0f);
   for(unsigned int i = 0; i < mesh.edges.size(); ++i) 
   {
      float dv = dot(dt * gravity, mesh.dual_edge_normals[i]);
      face_velocities[i] += dv;
   }
  
   
}

// -------------------------------------------------------------

void DualFluidSim::volume_conservation( std::vector<Vec2d>& vertex_positions, float dt )
{
   // compute current volume
   float current_volume = (float)explicit_surface.get_area();
   float current_length = (float)explicit_surface.get_non_solid_length();
   
   //for some reason these are negative.
   target_volume = fabs(target_volume);
   current_volume = fabs(current_volume);

   //use a small fraction to gradually reintroduce volume.
   //could couple this to the timestep, but not too important.
   float dV = 0.04f*( target_volume - current_volume ) / current_length;
   
   for ( unsigned int i = 0; i < explicit_surface.m_positions.size(); ++i )
   {
      if( explicit_surface.m_mesh.vtxedge[i].size() == 0 ) { continue; }
      if( explicit_surface.m_masses[i] > 1.5 ) { continue; } //skip solid vertices
      
      Vec2d n = explicit_surface.get_vertex_normal( i );
      Vec2f average_normal( n );
      vertex_positions[i] -= Vec2d(dV * average_normal);
   }
   
}

   
// -------------------------------------------------------------


bool DualFluidSim::trace_rk2(const Vec2f& start, Vec2f& end, float dt, const VelocityFunctor& get_velocity ) {
   
   //Advect a single point, with an unknown starting triangle

   Vec2f vel = get_velocity(start); 

   //advance to midpoint
   end = start + 0.5f * dt * vel;

   vel = get_velocity(end);  

   //compute the final position
   end = start + dt * vel;
   
   return true;
}

// -------------------------------------------------------------

bool DualFluidSim::trace_rk3(const Vec2f& start, Vec2f& end, float dt, const VelocityFunctor& get_velocity ) 
{
   
   Vec2f x = start;
   
   Vec2f k1 = dt * get_velocity( start );
   Vec2f k2 = dt * get_velocity( start + 0.5f * k1 );
   Vec2f k3 = dt * get_velocity( start - k1 + 3.0f * k2 );
   end = start + 1.0f/6.0f * ( k1 + 4.0f * k2 + k3 );
   
   return true;
}

// -------------------------------------------------------------


void DualFluidSim::advect_semilagrangian( float dt, 
                                          const std::vector<Vec2f>& starting_points, 
                                          std::vector<float>& new_velocities, 
                                          const std::vector<Vec2f>& new_edge_normals,
                                          const VelocityFunctor& get_velocity ) 
{
   
   new_velocities.clear();
   new_velocities.assign(starting_points.size(), 0);
   
   for(unsigned int i = 0; i < starting_points.size(); ++i) 
   {
      Vec2f end;
      
      trace_rk3( starting_points[i], end, -dt, get_velocity );

      // interpolate new velocity at the point
      Vec2f vel = get_velocity( end );

      //dot with the new face normal
      float component = dot(vel, new_edge_normals[i]);
      
      new_velocities[i] = component;
   }
   
   
}


// -------------------------------------------------------------


void DualFluidSim::solve_pressure() {

   std::vector<float> wall_velocities(face_velocities.size(), 0);
   for(unsigned int i = 0; i < wall_velocities.size(); ++i) {
      wall_velocities[i] = 0;
   }

   pressures = pressure_solve_voronoi(mesh, explicit_surface, face_velocities, solid_weights, liquid_phi, wall_velocities);
}


// -------------------------------------------------------------

Vec2f DualFluidSim::get_velocity_voronoi_consistent(const Vec2f& point ) 
{
   assert( voronoi_vertex_velocities.size() == mesh.tri_circumcentres.size() );
   
   int closest_vert = mesh.get_nearest_vertex(point);
   
   if ( closest_vert < 0 )
   {
      return Vec2f(0,0);
   }
   
   std::vector<int> polygon_vertex_indices;
   std::vector<float> polygon_edge_normal_velocities;   
   
   std::vector<Vec2f> polygon_vertex_locations;
   std::vector<Vec2f> polygon_vertex_velocities;( polygon_vertex_indices.size() );   

   if(per_vert_poly_verts[closest_vert].size() == 0 ) {

      construct_voronoi_polygon( closest_vert, polygon_vertex_indices, polygon_edge_normal_velocities );
      
      if ( polygon_vertex_indices.size() == 0 )
      {
         return Vec2f(0,0);
      }
      polygon_vertex_locations.resize( polygon_vertex_indices.size() );
      polygon_vertex_velocities.resize( polygon_vertex_indices.size() );
      
      for ( unsigned int i = 0; i < polygon_vertex_indices.size(); ++i )
      {
         polygon_vertex_locations[i] = mesh.tri_circumcentres[ polygon_vertex_indices[i] ];
         polygon_vertex_velocities[i] = voronoi_vertex_velocities[ polygon_vertex_indices[i] ];    // possibly extrapolated velocity
      }
    
      
      bool short_edge_found = true;

      std::vector<Vec2f> filtered_poly_points = polygon_vertex_locations;
      std::vector<Vec2f> filtered_data = polygon_vertex_velocities;
      
      while ( short_edge_found )
      {
         short_edge_found = false;
         std::vector<Vec2f> new_filtered_poly_points;
         std::vector<Vec2f> new_filtered_data;
         
         for ( unsigned int i = 0; i < filtered_poly_points.size(); ++i )
         {
            unsigned int next = ( i + 1 ) % filtered_poly_points.size();
            float edge_length = dist( filtered_poly_points[next], filtered_poly_points[i] );
            if ( edge_length > 1e-5 )
            {
               new_filtered_poly_points.push_back( filtered_poly_points[next] );
               new_filtered_data.push_back( filtered_data[next] );
            }
            else
            {
               short_edge_found = true;
            }
         }
         
         filtered_poly_points = new_filtered_poly_points;
         filtered_data = new_filtered_data;
      }
      
      for ( unsigned int i = 0; i < filtered_data.size(); ++i )
      {
         assert( filtered_data[i][0] == filtered_data[i][0] );
         assert( filtered_data[i][1] == filtered_data[i][1] );
      }
      
      if ( filtered_poly_points.size() < 3 )
      {
         return Vec2f( 0, 0 );
      }

      // check we have no short edges
      for ( unsigned int i = 0; i < filtered_poly_points.size(); ++i )
      {
         unsigned int next = ( i + 1 ) % filtered_poly_points.size();
         float edge_length = dist( filtered_poly_points[next], filtered_poly_points[i] );
         assert ( edge_length > 1e-5 );
      }

      polygon_vertex_locations = filtered_poly_points;
      polygon_vertex_velocities = filtered_data;

      per_vert_poly_verts[closest_vert] = polygon_vertex_locations;
      per_vert_poly_vels[closest_vert] = polygon_vertex_velocities;

   }
   else {
      polygon_vertex_locations = per_vert_poly_verts[closest_vert];
      polygon_vertex_velocities = per_vert_poly_vels[closest_vert];
   }
   
   Vec2f interpolated_velocity = triangulated_polygon_interpolation(point, polygon_vertex_locations, polygon_vertex_velocities);

   return interpolated_velocity;
   
}



// -------------------------------------------------------------

//
// Detect if triangles are adjacent and are oriented clockwise around the given central vertex.
//

bool DualFluidSim::triangles_are_adjacent_ordered( unsigned int tri1_index, 
                                                   unsigned int tri2_index, 
                                                   const unsigned int central_vertex, 
                                                   int& shared_edge_index )
{
    
   const Vec3ui& tri1_edges = mesh.tri_to_edge_map[tri1_index];
   const Vec3ui& tri2_edges = mesh.tri_to_edge_map[tri2_index];
   
   shared_edge_index = -1;
   for ( unsigned int i = 0; i < 3; ++i )
   {
      for ( unsigned int j = 0; j < 3; ++j )
      {
         if ( tri1_edges[i] == tri2_edges[j] )
         {
            shared_edge_index = tri1_edges[i];
            break;
         }
      }
   }

   if ( shared_edge_index < 0 )
   {
      return false;
   }
    
   // Okay, they are adjacent.  Now check the orientation of the 
   // triangle wrt the common edge and central vertex.
   
   const Vec2ui& common_edge = mesh.edges[shared_edge_index];
      
   unsigned int other_vert = ~0;
   if ( common_edge[0] == central_vertex )
   {
      other_vert = common_edge[1];
   }
   else
   {
      assert( common_edge[1] == central_vertex );
      other_vert = common_edge[0];      
   }
   
   const Vec3ui& tri1 = mesh.tris[ tri1_index ];
   
   if ( tri1[0] == central_vertex ) 
   { 
      if ( tri1[1] == other_vert ) { return false; }
      assert( tri1[2] == other_vert );
   }
   else if ( tri1[1] == central_vertex )
   {
      if ( tri1[2] == other_vert ) { return false; }
      assert( tri1[0] == other_vert );      
   }
   else
   {
      assert( tri1[2] == central_vertex );
      if ( tri1[0] == other_vert ) { return false; }
      assert( tri1[1] == other_vert );      
   }
   
   return true;
   
}

// -------------------------------------------------------------


void DualFluidSim::construct_voronoi_polygon( unsigned int voronoi_site_index, 
                                              std::vector<int>& polygon_vertices,
                                              std::vector<float>& polygon_edge_velocities )
{
   
   polygon_vertices.clear();
   
   std::vector<int> incident_tris = mesh.vert_to_tri_map[voronoi_site_index];    // take a copy of this list
   
   if ( incident_tris.size() < 2 ) 
   {
      return; 
   }
   
   std::vector<int> sorted_incident_tris;   
   std::vector<int> sorted_edge_indices;
   
   sorted_incident_tris.push_back( incident_tris[0] );
   incident_tris.erase( incident_tris.begin() + 0 );
   
   while ( !incident_tris.empty() ) 
   {
      unsigned int i = 0;
      int shared_edge_index;
      for ( ; i < incident_tris.size(); ++i )
      {
         if ( triangles_are_adjacent_ordered( sorted_incident_tris.back(), incident_tris[i], voronoi_site_index, shared_edge_index ) )
         {
            break;
         }
      }

      if( i == incident_tris.size() )
      {
         // the vertex's triangle neighbourhood is not connected and closed
         return;
      }
      
      sorted_incident_tris.push_back( incident_tris[i] );
      incident_tris.erase( incident_tris.begin() + i );      
      sorted_edge_indices.push_back( shared_edge_index );
   }
   
   // final edge index
   int shared_edge_index;
   bool check = triangles_are_adjacent_ordered( sorted_incident_tris.back(), sorted_incident_tris.front(), voronoi_site_index, shared_edge_index );

   if ( !check )
   {
      // the vertex's triangle neighbourhood is not connected and closed
      return;
   }
   
   sorted_edge_indices.push_back( shared_edge_index ); 
   
   // verify
   for ( unsigned int j = 0; j < sorted_incident_tris.size() - 1; ++j )
   {
      int shared_edge_index = 0;
      assert( triangles_are_adjacent_ordered( sorted_incident_tris[j], sorted_incident_tris[j+1], voronoi_site_index, shared_edge_index ) );
      assert( shared_edge_index == sorted_edge_indices[j] );
   }
   
   std::vector<int> sorted_edge_signs;
   
   for ( unsigned int j = 0; j < sorted_edge_indices.size(); ++j )
   {   
      for ( unsigned int i = 0; i < mesh.vert_to_edge_map[voronoi_site_index].size(); ++i )
      {
         if ( sorted_edge_indices[j] == mesh.vert_to_edge_map[voronoi_site_index][i] )
         {
            sorted_edge_signs.push_back( mesh.dual_edge_signs[voronoi_site_index][i] );
            break;
         }
      }
   }
   
   
   //
   // construct polygon
   //
      
   polygon_vertices = sorted_incident_tris;   

   polygon_edge_velocities.clear();
   
   for ( unsigned int j = 0; j < polygon_vertices.size(); ++j )
   {
      polygon_edge_velocities.push_back( (float)sorted_edge_signs[j] * face_velocities[sorted_edge_indices[j]] );      
   }

   
}


// -------------------------------------------------------------

//
// reconstruct velocities on polygon vertices
//
/*

void DualFluidSim::reconstruct_velocities_for_polygon( const std::vector<Vec2f>& polygon_vertex_locations, 
                                                       const std::vector<float>& polygon_edge_velocities,
                                                       std::vector<Vec2f>& polygon_vertex_velocities ) 
{
   
   assert( polygon_vertex_locations.size() == polygon_edge_velocities.size() );
   assert( polygon_vertex_locations.size() > 2 );
   
   polygon_vertex_velocities.clear();
   polygon_vertex_velocities.reserve( polygon_vertex_locations.size() );
   
   for ( unsigned int j = 0; j < polygon_vertex_locations.size(); ++j )
   {
      int jprev = (j-1);
      if ( j == 0 ) { jprev = polygon_vertex_locations.size() - 1; }
      int jnext = (j+1) % polygon_vertex_locations.size();
      
      Vec2f n0 = perp( polygon_vertex_locations[j] - polygon_vertex_locations[jprev] );
      float v0 = polygon_edge_velocities[jprev];
      Vec2f n1 = perp( polygon_vertex_locations[jnext] - polygon_vertex_locations[j] );
      float v1 = polygon_edge_velocities[j];
            
      Vec2f reconstructed_velocity;
      
      if ( mag( n0 ) == 0 && mag( n1 ) != 0 ) 
      {
         reconstructed_velocity = normalized( n1 ) * v1;
      }
      else if ( mag( n0 ) != 0 && mag( n1 ) == 0 ) 
      {
         reconstructed_velocity = normalized( n0 ) * v0;
      }
      else if ( mag( n0 ) == 0 && mag( n1 ) == 0 ) 
      {
         reconstructed_velocity = Vec2f(0,0);
      }
      else
      {         
         reconstructed_velocity = reconstruct( normalized(n0), v0, normalized(n1), v1 );
      }
      
      polygon_vertex_velocities.push_back( reconstructed_velocity );      
      
   }
   
}
*/

// -------------------------------------------------------------

void DualFluidSim::extrapolate_consistent_vertex_velocities()
{
   voronoi_vertex_velocities.clear();
   voronoi_vertex_velocities.resize( mesh.tris.size(), Vec2f(0,0) );
   
   voronoi_vertex_is_valid.clear();
   voronoi_vertex_is_valid.resize( mesh.tris.size(), false );
   
   voronoi_edge_is_valid.clear();
   voronoi_edge_is_valid.resize( mesh.edges.size(), false );

   // for each voronoi vertex
   for ( unsigned int i = 0; i < mesh.tris.size(); ++i )
   {
      
      std::vector<Vec2f> incident_edge_normals;
      std::vector<float> incident_edge_velocities;
      
      // get the vertex's three incident voronoi edges
      const Vec3ui& incident_edges = mesh.tri_to_edge_map[i];
      
      // for each of the voronoi edges
      for ( unsigned int j = 0; j < 3; ++j )
      {
         unsigned int edge_index = incident_edges[j];
         
         // check the cells on either side of the voronoi edge
         // if one of the cells is liquid, *and* if the edge itself has positive weight, 
         // we can use the voronoi edge

         int cell_a = mesh.edges[edge_index][0];
         int cell_b = mesh.edges[edge_index][1];
         
         if ( ( cell_a >= 0 && liquid_phi[ cell_a ] < 0.0f ) || 
              ( cell_b >= 0 && liquid_phi[ cell_b ] < 0.0f ) )
         {
            if ( solid_weights[edge_index] > 0.0 )
            {
               incident_edge_normals.push_back( mesh.dual_edge_normals[edge_index] );
               incident_edge_velocities.push_back( face_velocities[edge_index] );
               voronoi_edge_is_valid[edge_index] = true;
            }
         } 
      }
      
      if ( incident_edge_normals.size() < 2 ) { continue; }
      assert( incident_edge_normals.size() < 4 );
      
      bool success;
      voronoi_vertex_velocities[i] = reconstruct( incident_edge_normals, incident_edge_velocities, success );
      voronoi_vertex_is_valid[i] = success;
   }
   
   
   std::vector<bool> new_valid = voronoi_vertex_is_valid;
   
    
   //// now extrapolate across voronoi edges, by pulling valid neighbour velocities from invalid points, and taking the average
   for ( unsigned int pass = 0; pass < 10; ++pass )
   {
      for ( unsigned int i = 0; i < mesh.tris.size(); ++i )
      {
         //if ( solid_phi[i] > 0.0f ) { continue; }
         
         if ( !voronoi_vertex_is_valid[i] )
         {
            // find neighbouring voronoi vertices which are not marked as valid
            
            const Vec3ui& edge_ids = mesh.tri_to_edge_map[i];
            voronoi_vertex_velocities[i] = Vec2f(0,0);
            
            int valid_count = 0;
            for ( int neighbour = 0; neighbour < 3; ++neighbour ) 
            {
               // determine neighbour triangle
               int edge = edge_ids[neighbour];
               const Vec2i& tri_ids = mesh.edge_to_tri_map[edge];
               int nbr_tri = (tri_ids[0] == (int)i? tri_ids[1] : tri_ids[0]);
               if(nbr_tri == -1)
                  continue;
               
               if ( voronoi_vertex_is_valid[nbr_tri] )
               {
                  ++valid_count;
                  voronoi_vertex_velocities[i] += voronoi_vertex_velocities[nbr_tri];
                  new_valid[i] = true;
               }
            }
            if(valid_count > 0)
               voronoi_vertex_velocities[i] /= (float)valid_count;

         }
      }
      
      voronoi_vertex_is_valid = new_valid;
   }

   //get rid of non-tangential components for good measure
   for(unsigned int i = 0; i < voronoi_vertex_velocities.size(); ++i) {
      // project out u.n for faces inside the solid
      if ( solid_func(mesh.tri_circumcentres[i]) > 0.0 )
      {
         Vec2f solid_normal = solid_func_gradient( mesh.tri_circumcentres[i] );
         float udotn = dot( voronoi_vertex_velocities[i], solid_normal );
         voronoi_vertex_velocities[i] -= solid_normal * udotn; 
      }
   }
}


// -------------------------------------------------------------







