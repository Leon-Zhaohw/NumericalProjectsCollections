// ---------------------------------------------------------
//
//  
//
//  
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include <queue>

#include "surftrack.h"
#include "gluvi.h"
#include "broadphasegrid.h"
#include "ccd_wrapper.h"

// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

namespace eltopo2d 
{
   
// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------


SurfTrack::SurfTrack( unsigned int num_vertices,
           const double* vertices,
           unsigned int num_edges,
           const int* edges,
           const double* masses,
           const SurfTrackInitializationParameters& initial_parameters ) :
   DynamicSurface( num_vertices, vertices, num_edges, edges, masses, initial_parameters.m_proximity_epsilon, initial_parameters.m_collision_safety, initial_parameters.m_verbose ),
   m_max_area_change( initial_parameters.m_max_area_change ),
   m_min_edge_length( initial_parameters.m_min_edge_length ),
   m_max_edge_length( initial_parameters.m_max_edge_length ),
   m_merge_proximity_epsilon( initial_parameters.m_merge_proximity_epsilon ),
   m_allow_topology_changes( initial_parameters.m_allow_topology_changes ),
   m_perform_improvement( initial_parameters.m_perform_improvement )
{
}
   
// ---------------------------------------------------------
// 
// ---------------------------------------------------------
   
SurfTrack::SurfTrack( const std::vector<Vec2d>& vs, 
                      const std::vector<Vec2ui>& es, 
                      const std::vector<double>& masses,
                      const SurfTrackInitializationParameters& initial_parameters ) :
   DynamicSurface( vs, es, masses, initial_parameters.m_proximity_epsilon, initial_parameters.m_collision_safety, initial_parameters.m_verbose ),
   m_max_area_change( initial_parameters.m_max_area_change ),
   m_min_edge_length( initial_parameters.m_min_edge_length ),
   m_max_edge_length( initial_parameters.m_max_edge_length ),
   m_merge_proximity_epsilon( initial_parameters.m_merge_proximity_epsilon ),
   m_allow_topology_changes( initial_parameters.m_allow_topology_changes ),
   m_perform_improvement( initial_parameters.m_perform_improvement )
{
}


// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

void SurfTrack::render( unsigned int options )
{
   // vertices
   glBegin(GL_POINTS);
   for ( unsigned int i = 0; i < m_positions.size(); ++i)
   {
      glVertex2dv( m_positions[i].v );
   }
   glEnd();

   // edges
   glBegin(GL_LINES);
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i)
   {
      const Vec2ui& e = m_mesh.edges[i];
      glVertex2dv( m_positions[e[0]].v );
      glVertex2dv( m_positions[e[1]].v );
   }
   glEnd();
   
}
   

static void generate_new_midpoint( unsigned int edge_index, const SurfTrack& surf, Vec2d& new_vertex_smooth_position )
{
   const Vec2ui& e = surf.m_mesh.edges[edge_index];
   new_vertex_smooth_position = 0.5 * ( surf.m_positions[ e[0] ] + surf.m_positions[ e[1] ] );   
}
   
// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

bool SurfTrack::split_edge( unsigned int edge_index )
{
   if ( m_verbose ) { std::cout << "splitting edge: " << edge_index << std::endl; }
   
   Vec2ui original_edge = m_mesh.edges[edge_index];
   
   unsigned int vertex_a = m_mesh.edges[edge_index][0];
   unsigned int vertex_b = m_mesh.edges[edge_index][1];
      
   if ( m_verbose ) { std::cout << vertex_a << ", " << vertex_b << std::endl; }
   
   Vec2d new_vertex_position = 0.5 * ( m_positions[ vertex_a ] + m_positions[ vertex_b ] );
   Vec2d new_vertex_smooth_position;
   
   generate_new_midpoint( edge_index, *this, new_vertex_smooth_position );
   
   bool use_smooth_point = ! ( split_edge_pseudo_motion_introduces_collision( new_vertex_position, 
                                                                              new_vertex_smooth_position, 
                                                                              edge_index, 
                                                                              vertex_a, 
                                                                              vertex_b ) );
   
   // check normal inversion?
   if ( use_smooth_point )
   {
   }
   
   if ( use_smooth_point == false )
   {
      if ( m_verbose ) { std::cout << "not using smooth subdivision\n" << std::endl; } 
      
      new_vertex_smooth_position = new_vertex_position;
      
      if ( split_edge_pseudo_motion_introduces_collision( new_vertex_position, 
                                                          new_vertex_smooth_position, 
                                                          edge_index, 
                                                          vertex_a, 
                                                          vertex_b ) )
      {
         if ( m_verbose )  { std::cout << "Even mid-point subdivision introduces collision.  Backing out." << std::endl; }
         return false;
      }
   }
 	else
   {
      if ( m_verbose ) { std::cout << "using smooth subdivision\n" << std::endl; }
   }
   
   Vec2d new_vertex_velocity = 0.5 * ( m_velocities[ vertex_a ] + m_velocities[ vertex_b ] );
   
   assert( m_volume_ids[ vertex_a ] == m_volume_ids[ vertex_b ] );
   
   double new_vertex_mass = min( m_masses[ vertex_a ], m_masses[ vertex_b ] );
   
   unsigned int new_vertex_index = add_vertex( new_vertex_smooth_position, new_vertex_velocity, new_vertex_mass, m_volume_ids[ vertex_a ] );

   remove_edge( edge_index );
      
   add_edge( Vec2ui( vertex_a, new_vertex_index ) );
   add_edge( Vec2ui( new_vertex_index, vertex_b ) );
   
   return true;
   
}


   
// ---------------------------------------------------------
// topological merging
// ---------------------------------------------------------

// ---------------------------------------------------------
// 
//
//
// ---------------------------------------------------------
   
bool SurfTrack::get_zipper_edges( unsigned int vertex_index, unsigned int edge_index, std::vector<Vec2ui>& new_edges )
{
   const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[vertex_index];
   assert( incident_edges.size() == 2 );

   assert( vertex_index == m_mesh.edges[ incident_edges[0] ][0] || vertex_index == m_mesh.edges[ incident_edges[0] ][1] );
   assert( vertex_index == m_mesh.edges[ incident_edges[1] ][0] || vertex_index == m_mesh.edges[ incident_edges[1] ][1] );
   
   unsigned int v_neighbour_a = m_mesh.edges[ incident_edges[0] ][0];          
   unsigned int v_neighbour_b = m_mesh.edges[ incident_edges[1] ][1];;
   Vec2ui e = m_mesh.edges[edge_index];
   if ( v_neighbour_a == vertex_index )
   {
      v_neighbour_a = m_mesh.edges[ incident_edges[1] ][0];
      v_neighbour_b = m_mesh.edges[ incident_edges[0] ][1];      
      assert( vertex_index == m_mesh.edges[ incident_edges[1] ][1] ); 
      assert( vertex_index == m_mesh.edges[ incident_edges[0] ][0] ); 
   }
   else
   {
      assert( vertex_index == m_mesh.edges[ incident_edges[0] ][1] ); 
      assert( vertex_index == m_mesh.edges[ incident_edges[1] ][0] ); 
   }
   
   assert( v_neighbour_a != vertex_index );
   assert( v_neighbour_b != vertex_index );

   if ( e[0] == v_neighbour_a || e[0] == v_neighbour_b || e[1] == v_neighbour_a || e[1] == v_neighbour_b )
   {
      return false;
   }
   
   Vec2ui new_edge = Vec2ui( v_neighbour_a, e[1] );
   if ( m_mesh.get_edge( new_edge ) == m_mesh.edges.size() )
   {
      new_edges.push_back( new_edge );
   }

   new_edge = Vec2ui( e[0], v_neighbour_b );
   if ( m_mesh.get_edge( new_edge ) == m_mesh.edges.size() )
   {
      new_edges.push_back( new_edge );
   }
      
                            
   return true;
}
 
   
// ---------------------------------------------------------
// 
//
//
// ---------------------------------------------------------
   
bool SurfTrack::zippering_introduces_collision( const std::vector<Vec2ui>& new_edges, const std::vector<unsigned int>& deleted_edges )
{

   // check between new edges
   
   for ( unsigned int new_e = 0; new_e < new_edges.size(); ++new_e )
   {
      const Vec2ui& new_edge = new_edges[new_e];
      
      for ( unsigned int other_new_e = new_e+1; other_new_e < new_edges.size(); ++other_new_e )
      {
         const Vec2ui& other_new_edge = new_edges[other_new_e];
         
         if( new_edge[0] != other_new_edge[0] && new_edge[0] != other_new_edge[1] && new_edge[1] != other_new_edge[0] && new_edge[1] != other_new_edge[1] )
         {
            if ( segment_segment_intersection( m_positions[new_edge[0]], new_edge[0],
                                              m_positions[new_edge[1]], new_edge[1], 
                                              m_positions[other_new_edge[0]], other_new_edge[0],
                                              m_positions[other_new_edge[1]], other_new_edge[1] ) )
            {
               return true;
            }
         }
         
      }
      
   }
      
   // check vs. existing edges
   
   for ( unsigned int new_e = 0; new_e < new_edges.size(); ++new_e )
   {
      const Vec2ui& new_edge = new_edges[new_e];
      Vec2d emin, emax;
      minmax( m_positions[new_edge[0]], m_positions[new_edge[1]], emin, emax );
             
      emin -= m_merge_proximity_epsilon * Vec2d(1,1);
      emax += m_merge_proximity_epsilon * Vec2d(1,1);

      std::vector<unsigned int> edge_candidates;
      m_broad_phase->get_potential_edge_collisions( emin, emax, edge_candidates );
      
      for(unsigned int i = 0; i < edge_candidates.size(); ++i )
      {
         bool is_deleted = false;
         for ( unsigned int j = 0; j < deleted_edges.size(); ++j )
         {
            if ( edge_candidates[i] == deleted_edges[j] )
            {
               is_deleted = true;
               break;
            }
         }
         if ( is_deleted ) { continue; }
            
         const Vec2ui& existing_edge = m_mesh.edges[ edge_candidates[i] ];
         
         if( new_edge[0] != existing_edge[0] && new_edge[0] != existing_edge[1] && new_edge[1] != existing_edge[0] && new_edge[1] != existing_edge[1] )
         {
            if ( segment_segment_intersection( m_positions[new_edge[0]], new_edge[0],
                                               m_positions[new_edge[1]], new_edge[1], 
                                               m_positions[existing_edge[0]], existing_edge[0],
                                               m_positions[existing_edge[1]], existing_edge[1] ) )
            {
               return true;
            }
         }
      }
   }
   
   return false;
}


   
double triangle_area( const Vec2d& a, const Vec2d& b, const Vec2d& c )
{
   return 0.5 * fabs( cross( b - a, c - a ) );
}     
   
// ---------------------------------------------------------
/// 
/// Control volume change
///
// ---------------------------------------------------------

   
bool SurfTrack::zippering_introduces_area_change( unsigned int vertex_index, unsigned int edge_index )
{   
   
   const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[vertex_index];
   assert( incident_edges.size() == 2 );

   unsigned int v_neighbour_a = m_mesh.edges[ incident_edges[0] ][0];          
   unsigned int v_neighbour_b = m_mesh.edges[ incident_edges[1] ][1];;

   if ( v_neighbour_a == vertex_index )
   {
      v_neighbour_a = m_mesh.edges[ incident_edges[1] ][0];
      v_neighbour_b = m_mesh.edges[ incident_edges[0] ][1];      
      assert( vertex_index == m_mesh.edges[ incident_edges[1] ][1] ); 
      assert( vertex_index == m_mesh.edges[ incident_edges[0] ][0] ); 
   }
   else
   {
      assert( vertex_index == m_mesh.edges[ incident_edges[0] ][1] ); 
      assert( vertex_index == m_mesh.edges[ incident_edges[1] ][0] ); 
   }
   
   assert( v_neighbour_a != vertex_index );
   assert( v_neighbour_b != vertex_index );
   
   unsigned int ea = m_mesh.edges[ edge_index ][0];
   unsigned int eb = m_mesh.edges[ edge_index ][1];
   
   double area = 0.0;
   
   area += triangle_area( m_positions[ v_neighbour_a ], m_positions[ eb ], m_positions[ vertex_index ] ); 
   area += triangle_area( m_positions[ vertex_index ], m_positions[ eb ], m_positions[ ea ] );
   area += triangle_area( m_positions[ vertex_index ], m_positions[ eb ], m_positions[ v_neighbour_b ] );
   
   return area > m_max_area_change;
                      
}
   
// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

bool SurfTrack::merge( unsigned int vertex_index, unsigned int edge_index )
{
   if ( m_mesh.vtxedge[vertex_index].size() != 2  )
   {
      if ( m_verbose ) { std::cout << "ZIPPER: vertex non-manifold" << std::endl; }
      return false;
   }
   
   if (    ( m_masses[m_mesh.edges[edge_index][0]] != m_masses[vertex_index] )
        || ( m_masses[m_mesh.edges[edge_index][1]] != m_masses[vertex_index] ) )
   {
      return false;
   }
   
   //
   // Get the 2 new edges which will join the two holes in the mesh
   //
   
   std::vector<Vec2ui> new_edges;

   if ( false == get_zipper_edges( vertex_index, edge_index, new_edges ) )
   {
      if ( m_verbose ) { std::cout << "ZIPPER: couldn't get a set of edges" << std::endl;   }
      return false;
   }
   
   // Keep a list of edges to delete
   std::vector<unsigned int> deleted_edges;
   deleted_edges.push_back( edge_index );
   deleted_edges.push_back( m_mesh.vtxedge[vertex_index][0] );
   deleted_edges.push_back( m_mesh.vtxedge[vertex_index][1] );
   
   //
   // Check the new edges for collision safety, ignoring the edges which will be deleted
   //
   
   if ( zippering_introduces_collision( new_edges, deleted_edges ) )
   {
      if ( m_verbose ) { std::cout << "ZIPPER: collision check failed" << std::endl; }
      return false;
   }
      
   //
   // Add the new edges
   //
   
   for ( unsigned int i = 0; i < new_edges.size(); ++i )
   {
      if ( m_verbose ) { std::cout << "adding edge " << new_edges[i] << std::endl; }
      add_edge( new_edges[i] );
   }
   

   //
   // Remove the old edges
   //
      
   remove_edge( deleted_edges[0] );
   remove_edge( deleted_edges[1] );
   remove_edge( deleted_edges[2] );
      
   remove_vertex( vertex_index );
   
   return true;
   
}


// ---------------------------------------------------------
// mesh cleanup
// ---------------------------------------------------------

// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

void SurfTrack::trim_non_manifold( const std::vector<unsigned int>& edge_indices )
{
}

// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

void SurfTrack::trim_non_manifold()
{
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      // ignore already deleted
      
      if ( m_mesh.edges[i][0] == m_mesh.edges[i][1] ) { continue; }
      
      // single edges
      
      Vec2ui e = m_mesh.edges[i];
      if ( ( m_mesh.vtxedge[e[0]].size() == 1 ) && ( m_mesh.vtxedge[e[1]].size() == 1 ) )
      {
         if ( m_verbose ) { std::cout << "removing single edge" << std::endl; }
         remove_edge( i );
         remove_vertex( e[0] );
         remove_vertex( e[1] );
         continue;
      }
      
      // two edges with same vertices
      
      unsigned int vertex_a = m_mesh.edges[i][0];
      unsigned int vertex_b = m_mesh.edges[i][1];
      std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[vertex_a];
      
      for ( unsigned int j = 0; j < incident_edges.size(); ++j )
      {
         if ( incident_edges[j] == i ) { continue; }
         const Vec2ui& another_edge = m_mesh.edges[ incident_edges[j] ];
         if ( another_edge[0] == vertex_b || another_edge[1] == vertex_b )
         {
            assert( another_edge[0] == vertex_a || another_edge[1] == vertex_a );
            if ( m_verbose ) { std::cout << "removing double edge: " << i << ", " << incident_edges[j] << std::endl; }
            if ( m_verbose ) { std::cout << "vertices: " << vertex_a << ", " << vertex_b << std::endl; }
            remove_edge(incident_edges[j]);
            remove_edge(i);
            break;
         }
      }
   }
   
   if ( m_verbose ) 
   {
      for ( unsigned int i = 0; i < m_mesh.vtxedge.size(); ++i )
      {
         if( m_mesh.vtxedge[i].size() != 0 && m_mesh.vtxedge[i].size() != 2 )
         {
            std::cout << "vertex incident on: " << std::endl;
            for ( unsigned int j = 0; j < m_mesh.vtxedge[i].size(); ++j )
            {
               std::cout << m_mesh.edges[ m_mesh.vtxedge[i][j] ] << std::endl;
            }
         }
      }
   }
   
}

// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

void SurfTrack::separate_singular_vertices()
{
   
}

// ---------------------------------------------------------
// mesh maintenance helpers
// ---------------------------------------------------------

// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

bool SurfTrack::split_edge_pseudo_motion_introduces_collision( const Vec2d& new_vertex_position, 
                                                               const Vec2d& new_vertex_smooth_position, 
                                                               unsigned int edge, 
                                                               unsigned int vertex_a, 
                                                               unsigned int vertex_b )
{
   if ( !m_collision_safety)
   {
      return false;
   }
   
   unsigned int dummy_index = m_positions.size();
   
   // new point vs all edges
   {
      
      Vec2d aabb_low, aabb_high;
      minmax( new_vertex_position, new_vertex_smooth_position, aabb_low, aabb_high );
      
      aabb_low -= m_proximity_epsilon * Vec2d(1,1);
      aabb_high += m_proximity_epsilon * Vec2d(1,1);
      
      std::vector<unsigned int> overlapping_edges;
      m_broad_phase->get_potential_edge_collisions( aabb_low, aabb_high, overlapping_edges );
      
      for ( unsigned int i = 0; i < overlapping_edges.size(); ++i )
      {
         if ( overlapping_edges[i] == edge ) { continue; }
         
         unsigned int edge_vertex_0 = m_mesh.edges[ overlapping_edges[i] ][0];
         unsigned int edge_vertex_1 = m_mesh.edges[ overlapping_edges[i] ][1];
                  
         Vec2ui sorted_edge( edge_vertex_0, edge_vertex_1 );
         if ( sorted_edge[0] > sorted_edge[1] ) { swap( sorted_edge[0], sorted_edge[1] ); }
         
         if ( point_segment_collision(  new_vertex_position, new_vertex_smooth_position, dummy_index,
                                        m_positions[ sorted_edge[0] ], m_positions[ sorted_edge[0] ], sorted_edge[0],
                                        m_positions[ sorted_edge[1] ], m_positions[ sorted_edge[1] ], sorted_edge[1] ) )
            
         {
            return true;
         }
      }
   }
   
   // 2 new edges vs all points
   
   {
      
      Vec2d edge_aabb_low, edge_aabb_high;
      
      // do one big query into the broadphase for both new edges
      minmax( new_vertex_position, new_vertex_smooth_position, 
              m_positions[ vertex_a ], m_positions[ vertex_b ],
              edge_aabb_low, edge_aabb_high );
      
      edge_aabb_low -= m_proximity_epsilon * Vec2d(1,1);
      edge_aabb_high += m_proximity_epsilon * Vec2d(1,1);
      
      std::vector<unsigned int> overlapping_vertices;
      m_broad_phase->get_potential_vertex_collisions( edge_aabb_low, edge_aabb_high, overlapping_vertices );
      
      for ( unsigned int i = 0; i < overlapping_vertices.size(); ++i )
      {
         if ( vertex_a != overlapping_vertices[i] ) 
         {
            
            if ( point_segment_collision( m_positions[ overlapping_vertices[i] ], m_positions[ overlapping_vertices[i] ], overlapping_vertices[i],
                                          m_positions[ vertex_a ], m_positions[ vertex_a ], vertex_a,
                                          new_vertex_position, new_vertex_smooth_position, dummy_index ) )
               
            {      
               return true;
            }
         }
         
         if ( vertex_b != overlapping_vertices[i] ) 
         {
            
            if ( point_segment_collision( m_positions[ overlapping_vertices[i] ], m_positions[ overlapping_vertices[i] ], overlapping_vertices[i],
                                          m_positions[ vertex_b ], m_positions[ vertex_b ], vertex_b,
                                          new_vertex_position, new_vertex_smooth_position, dummy_index ) )
            {          
               return true;
            }
         }
      }
      
   }
   
   return false;
}



// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

bool SurfTrack::collapse_edge_introduces_collision( unsigned int vertex_to_delete, unsigned int vertex_to_keep, unsigned int edge, const Vec2d& vertex_new_position )
{
   
   assert( m_newpositions.size() > 0 );
   
   if ( !m_collision_safety )  { return false; }

   // Change source vertex predicted position to superimpose onto dest vertex
   m_newpositions[vertex_to_keep] = vertex_new_position;
   m_newpositions[vertex_to_delete] = vertex_new_position;
   
   update_continuous_broad_phase( vertex_to_keep );
   update_continuous_broad_phase( vertex_to_delete );
      
   // Get the set of edges which move because of this motion
   std::vector<unsigned int> moving_edges;
   for ( unsigned int i = 0; i < m_mesh.vtxedge[vertex_to_keep].size(); ++i )
   {
      moving_edges.push_back( m_mesh.vtxedge[vertex_to_keep][i] );
   }
   for ( unsigned int i = 0; i < m_mesh.vtxedge[vertex_to_delete].size(); ++i )
   {
      moving_edges.push_back( m_mesh.vtxedge[vertex_to_delete][i] );
   }
   
   // Check this set of edges for collisions, holding everything else static
   
   for ( unsigned int i = 0; i < moving_edges.size(); ++i )
   { 
            
      const Vec2ui& current_edge = m_mesh.edges[ moving_edges[i] ];
      
      // Test the edge vs all other edges
      Vec2d aabb_low, aabb_high;
      minmax( m_positions[ current_edge[0] ], m_positions[ current_edge[1] ], 
              m_newpositions[ current_edge[0] ], m_newpositions[ current_edge[1] ],
              aabb_low, aabb_high );
      
      std::vector<unsigned int> overlapping_edges;
      m_broad_phase->get_potential_edge_collisions( aabb_low, aabb_high, overlapping_edges );
      
      for ( unsigned j=0; j < overlapping_edges.size(); ++j )
      {
         // Don't check against edges which are incident to the dest vertex
         bool edge_incident_to_dest = false;
         for ( unsigned int k = 0; k < moving_edges.size(); ++k )
         {
            if ( moving_edges[k] == overlapping_edges[j] )
            {
               edge_incident_to_dest = true;
               break;
            }
         }
         if ( edge_incident_to_dest )    { continue; }
                  
         const Vec2ui& other_edge =  m_mesh.edges[ overlapping_edges[j] ];
         
         if (   current_edge[0] == other_edge[0] || current_edge[1] == other_edge[0] 
             || current_edge[0] == other_edge[1] || current_edge[1] == other_edge[1] )
         {
            continue;
         }
         
         if ( point_segment_collision( m_positions[current_edge[0]], m_newpositions[current_edge[0]], current_edge[0],
                                       m_positions[other_edge[0]],   m_newpositions[other_edge[0]], other_edge[0],
                                       m_positions[other_edge[1]],   m_newpositions[other_edge[1]], other_edge[1] ) )
         {
            if ( m_verbose ) {
               std::cout << "point_segment_collision A: " << current_edge[0] << ", " << other_edge[0] << ", " << other_edge[1] << std::endl;
            }
            return true;
         }

         if ( point_segment_collision( m_positions[current_edge[1]], m_newpositions[current_edge[1]], current_edge[1],
                                       m_positions[other_edge[0]],   m_newpositions[other_edge[0]], other_edge[0],
                                       m_positions[other_edge[1]],   m_newpositions[other_edge[1]], other_edge[1] ) )
         {
            if ( m_verbose ) {
                std::cout << "point_segment_collision B: " << current_edge[1] << ", " << other_edge[0] << ", " << other_edge[1] << std::endl;
            }
            return true;
         }

         if ( point_segment_collision( m_positions[other_edge[0]],   m_newpositions[other_edge[0]], other_edge[0],
                                       m_positions[current_edge[0]], m_newpositions[current_edge[0]], current_edge[0],
                                       m_positions[current_edge[1]], m_newpositions[current_edge[1]], current_edge[1] ) )
         {
            if ( m_verbose ) {
                std::cout << "point_segment_collision C: " << other_edge[0] << ", " << current_edge[0] << ", " << current_edge[1] << std::endl;            
            }
            return true;
         }

         if ( point_segment_collision( m_positions[other_edge[1]],   m_newpositions[other_edge[1]], other_edge[0],
                                       m_positions[current_edge[0]], m_newpositions[current_edge[0]], current_edge[0],
                                       m_positions[current_edge[1]], m_newpositions[current_edge[1]], current_edge[1] ) )
         {
            if ( m_verbose ) {
               std::cout << "point_segment_collision D: " << other_edge[1] << ", " << current_edge[0] << ", " << current_edge[1] << std::endl;            
            }
            return true;
         }
         
      }
      
   }
   
   return false;
}
   
   
bool SurfTrack::collapse_edge_introduces_normal_inversion( unsigned int vertex_to_delete, unsigned int vertex_to_keep, unsigned int edge, const Vec2d& vertex_new_position ) const
{
   
   // check for normal inversion
      
   const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[vertex_to_keep];
   
   for ( unsigned int i = 0; i < incident_edges.size(); ++i )
   {
      if ( incident_edges[i] == edge ) { continue; }
      
      Vec2ui current_edge = m_mesh.edges[ incident_edges[i] ];
      
      Vec2d old_normal = get_edge_normal( current_edge );
      
      Vec2d vertex_a = m_positions[ current_edge[0] ];
      Vec2d vertex_b = m_positions[ current_edge[1] ];
      
      if ( current_edge[0] == vertex_to_keep )
      {
         vertex_a = vertex_new_position;
      }
      else
      {
         assert( current_edge[1] == vertex_to_keep );
         vertex_b = vertex_new_position;
      }
         
      Vec2d ev = vertex_b - vertex_a;
      ev /= mag(ev);
      
      Vec2d new_normal( -ev[1], ev[0] );
      
      if ( dot( new_normal, old_normal ) < 0.0 )
      {
         if ( m_verbose ) { std::cout << "collapse edge introduces normal inversion\n" << std::endl; }
         return true;
      } 
   }
   
   return false;
}


// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

bool SurfTrack::collapse_edge_introduces_area_change( unsigned int vertex_to_delete, unsigned int vertex_to_keep, unsigned int edge, const Vec2d& vertex_new_position )
{
   
   double area = 0;
   
   {
      Vec2d a = m_positions[ vertex_to_keep ];
      Vec2d b(0,0);
      const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[vertex_to_keep];
      for ( unsigned int i = 0; i < incident_edges.size(); ++i )
      {
         if ( incident_edges[i] == edge ) { continue; }
         
         const Vec2ui& incident_edge = m_mesh.edges[ incident_edges[i] ];
         if ( incident_edge[0] == vertex_to_keep )
         {
            b = m_positions[ incident_edge[1] ];
         }
         else
         {
            assert( incident_edge[1] == vertex_to_keep );
            b = m_positions[ incident_edge[0] ];
         }
      }     
      
      area += cross( a - vertex_new_position, b - vertex_new_position );
   }
   
   {
      Vec2d c = m_positions[ vertex_to_delete ];
      Vec2d d(0,0);
      const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[vertex_to_keep];
      for ( unsigned int i = 0; i < incident_edges.size(); ++i )
      {
         if ( incident_edges[i] == edge ) { continue; }
         
         const Vec2ui& incident_edge = m_mesh.edges[ incident_edges[i] ];
         if ( incident_edge[0] == vertex_to_keep )
         {
            d = m_positions[ incident_edge[1] ];
         }
         else
         {
            assert( incident_edge[1] == vertex_to_keep );
            d = m_positions[ incident_edge[0] ];
         }
      }
   
      area += cross( d - vertex_new_position, c - vertex_new_position );
   }
   
   if ( fabs(area) > m_max_area_change )
   {
      return true;
   }
   
   return false;
}


// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

bool SurfTrack::collapse_edge(unsigned int edge )
{
   unsigned int vertex_to_keep = m_mesh.edges[edge][0];
   unsigned int vertex_to_delete = m_mesh.edges[edge][1];
   
   // If we're disallowing topology changes, don't let an edge collapse form a degeneracy
   if ( false == m_allow_topology_changes )
   {
      std::vector<unsigned int>& vtx_a_edges = m_mesh.vtxedge[vertex_to_keep];
      std::vector<unsigned int>& vtx_b_edges = m_mesh.vtxedge[vertex_to_delete];
      
      std::vector<unsigned int> a_star, b_star;
      
      for ( unsigned int i = 0; i < vtx_a_edges.size(); ++i )
      {
         if ( m_mesh.edges[vtx_a_edges[i]][0] == vertex_to_keep )
         {
            a_star.push_back( m_mesh.edges[vtx_a_edges[i]][1] );
         }
         else
         {
            assert( m_mesh.edges[vtx_a_edges[i]][1] == vertex_to_keep );
            a_star.push_back( m_mesh.edges[vtx_a_edges[i]][0] );
         }
      }
      
      for ( unsigned int i = 0; i < vtx_b_edges.size(); ++i )
      {
         unsigned int b_neighbour = ~0;
         if ( m_mesh.edges[vtx_b_edges[i]][0] == vertex_to_delete )
         {
            b_neighbour = m_mesh.edges[vtx_b_edges[i]][1];
         }
         else
         {
            assert( m_mesh.edges[vtx_b_edges[i]][1] == vertex_to_delete );
            b_neighbour = m_mesh.edges[vtx_b_edges[i]][0];
         }
         
         for ( unsigned int j = 0; j < a_star.size(); ++j )
         {
            if ( a_star[j] == b_neighbour )
            {
               return false;
            }
            
         }
      }      
   }
   
   Vec2d vertex_new_position;

   if ( m_masses[vertex_to_keep] == m_masses[vertex_to_delete] )
   {
      generate_new_midpoint( edge, *this, vertex_new_position );
   }
   else if ( m_masses[vertex_to_keep] > m_masses[vertex_to_delete] )
   {
      vertex_new_position = m_positions[vertex_to_keep];
   }
   else
   {
      vertex_new_position = m_positions[vertex_to_delete];
   }
   
   if ( m_verbose ) { std::cout << "Collapsing edge.  Doomed vertex: " << vertex_to_delete << " --- Vertex to keep: " << vertex_to_keep << std::endl; }
   
   // Check vertex pseudo motion for collisions
   
   if ( mag ( m_positions[m_mesh.edges[edge][1]] - m_positions[m_mesh.edges[edge][0]] ) > 0 )
   {
      
      bool collision = collapse_edge_introduces_collision( vertex_to_delete, vertex_to_keep, edge, vertex_new_position );
      
      // Restore saved positions which were changed by the function we just called.
      m_newpositions[vertex_to_keep] = m_positions[vertex_to_keep];
      m_newpositions[vertex_to_delete] = m_positions[vertex_to_delete];
      
      if ( collision )
      {
         if ( m_verbose ) { std::cout << "edge collapse would introduce collision" << std::endl; }
         return false;
      }
      
      bool normal_inversion = collapse_edge_introduces_normal_inversion( vertex_to_delete, vertex_to_keep, edge, vertex_new_position );
      
      if ( normal_inversion )
      {
         if ( m_verbose ) { std::cout << "edge collapse would introduce normal inversion" << std::endl; }
         return false;         
      }
      
      bool area_change = collapse_edge_introduces_area_change( vertex_to_delete, vertex_to_keep, edge, vertex_new_position );
      
      if ( area_change )
      {
         if ( m_verbose ) { std::cout << "edge collapse would introduce change in area" << std::endl; }
         return false;         
      }
      
   }
   
   m_newpositions[vertex_to_keep] = m_positions[vertex_to_keep] = vertex_new_position;
   
   update_static_broad_phase( vertex_to_keep );
   
   // Find anything pointing to the doomed vertex and change it
   
   // copy the list of edges, don't take a refence to it
   std::vector< unsigned int > edges_incident_to_vertex = m_mesh.vtxedge[vertex_to_delete];    
   
   for ( unsigned int i=0; i < edges_incident_to_vertex.size(); ++i )
   {
      if ( edges_incident_to_vertex[i] == edge ) { continue; }
      
      Vec2ui new_edge = m_mesh.edges[ edges_incident_to_vertex[i] ];
      
      if ( new_edge[0] == vertex_to_delete )   { new_edge[0] = vertex_to_keep; }
      if ( new_edge[1] == vertex_to_delete )   { new_edge[1] = vertex_to_keep; }
      
      if ( m_verbose ) { std::cout << "adding updated edge: " << new_edge << std::endl; }
      
      unsigned int new_edge_index = add_edge( new_edge );
      (void) new_edge_index;
      
   }
   
   for ( unsigned int i=0; i < edges_incident_to_vertex.size(); ++i )
   {  
      if ( m_verbose )
      {
         std::cout << "removing vertex-incident edge: " << m_mesh.edges[ edges_incident_to_vertex[i] ] << std::endl;
      }
      
      remove_edge( edges_incident_to_vertex[i] );
   }
   
   // Delete vertex
   assert( m_mesh.vtxedge[vertex_to_delete].size() == 0 );
   remove_vertex( vertex_to_delete );
   
   return true;
   
}

   
// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

bool SurfTrack::split_pass()
{
   std::cout << "---------------------- El Topo: Edge split pass ----------------------" << std::endl;
   
   m_mesh.update_connectivity( m_positions.size() );
   rebuild_static_broad_phase();
   
   // whether a split operation was successful in this pass
   bool split_occurred = false;
   
   std::vector<unsigned int> edges_to_try;
   edges_to_try.clear();
   
   for( unsigned int i = 0; i < m_mesh.edges.size(); i++ )
   {    
      if ( m_mesh.edges[i][0] == m_mesh.edges[i][1] )   { continue; }    // skip deleted edges
      
      if ( m_masses[ m_mesh.edges[i][0] ] > 100.0 && m_masses[ m_mesh.edges[i][1] ] > 100.0 )     { continue; }    // skip solids
      
      unsigned int vertex_a = m_mesh.edges[i][0];
      unsigned int vertex_b = m_mesh.edges[i][1];
      
      assert( vertex_a < m_positions.size() );
      assert( vertex_b < m_positions.size() );
      
      Vec2d dp = m_positions[ vertex_b ] - m_positions[ vertex_a ];
      
      double current_length = mag( dp );
      
      if ( current_length > m_max_edge_length )
      {
         edges_to_try.push_back( i );
      }
   }
   
   while ( !edges_to_try.empty() )
   {
      unsigned int longest_edge_index = ~0;
      double longest_edge_length = -1.0;
      
      // find max edge length
      for( unsigned int i = 0; i < edges_to_try.size(); i++ )
      {
         unsigned int vertex_a = m_mesh.edges[edges_to_try[i]][0];
         unsigned int vertex_b = m_mesh.edges[edges_to_try[i]][1];
         Vec2d dp = m_positions[ vertex_b ] - m_positions[ vertex_a ];
         
         double current_length = mag( dp );
         
         if ( current_length > longest_edge_length )
         {
            longest_edge_length = current_length;
            longest_edge_index = i;
         } 
      }
      
      unsigned int longest_edge = edges_to_try[longest_edge_index];
      
      if ( longest_edge_length > m_max_edge_length )
      {
         
         if ( m_verbose ) 
         {
            std::cout << "splitting edge " << longest_edge << " / " << m_mesh.edges.size() << ", length = " << longest_edge_length << std::endl;
            std::cout << "edge " << m_mesh.edges[longest_edge];
         }
         
         if ( m_mesh.edges[longest_edge][0] == m_mesh.edges[longest_edge][1] )
         { 
            edges_to_try.erase( edges_to_try.begin() + longest_edge_index );
            continue;      // skip non-manifold and deleted edges
         }    
         
         bool result = split_edge( longest_edge );
         
         if ( m_verbose ) { std::cout << " result:  " << (result ? "ok" : "failed") << std::endl; }
         
         if (result)
         {
            split_occurred = true;
         }
      }
      
      edges_to_try.erase( edges_to_try.begin() + longest_edge_index );
      
   }
   
   return split_occurred;
   
}

   
// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

bool SurfTrack::collapse_pass()
{
   std::cout << "---------------------- El Topo: Edge collapse pass ----------------------" << std::endl;
   
   m_mesh.update_connectivity( m_positions.size() );

   rebuild_static_broad_phase();
   
   bool collapse_occurred = false;
   
   std::vector<SortableEdge> sortable_edges_to_try;
   
   for( unsigned int i = 0; i < m_mesh.edges.size(); i++ )
   {    
      if ( m_mesh.edges[i][0] == m_mesh.edges[i][1] )   { continue; }    // skip deleted edges
      if ( m_masses[ m_mesh.edges[i][0] ] > 100.0 )     { continue; }    // skip solids
      
      unsigned int vertex_a = m_mesh.edges[i][0];
      unsigned int vertex_b = m_mesh.edges[i][1];
      
      double current_length = mag( m_positions[ vertex_b ] - m_positions[ vertex_a ] );
      
      if ( current_length < m_min_edge_length )
      {
         sortable_edges_to_try.push_back( SortableEdge( i, current_length ) );
      }
   }
   
   std::sort( sortable_edges_to_try.begin(), sortable_edges_to_try.end() );
   
   if ( m_verbose )
   {
      std::cout << sortable_edges_to_try.size() << " candidate edges sorted" << std::endl;
      std::cout << "total edges: " << m_mesh.edges.size() << std::endl;
   }
   
   for ( unsigned int i = 0; i < sortable_edges_to_try.size(); ++i )
   {
      unsigned int e = sortable_edges_to_try[i].edge_index;
      double edge_length = mag( m_positions[m_mesh.edges[e][0]] - m_positions[m_mesh.edges[e][1]] );
      if ( edge_length < m_min_edge_length )
      {        
         if ( m_verbose )
         {
            std::cout << "collapsing edge " << e << " / " << m_mesh.edges.size() << ", length = " << edge_length << std::endl;
         }
         
         if ( m_mesh.edges[e][0] == m_mesh.edges[e][1] )
         {
            continue;      // skip deleted edges
         }    
         
         Vec2ui saved_edge = m_mesh.edges[e];
         
         bool result = collapse_edge( e );
         
         if ( m_verbose )
         {
            std::cout << " result: " << (result ? "ok" : "failed") << std::endl;
         }
         
         if ( result )
         {
            if ( m_verbose )
            {
               std::cout << "collapsed edge " << e << ": " << saved_edge << std::endl;
            }
            
            collapse_occurred = true;
            trim_non_manifold();
            rebuild_static_broad_phase();
         }
         
      }
   }
   
   trim_non_manifold();
   
   return collapse_occurred;
   
}

   
// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

void SurfTrack::merge_pass()
{
   std::cout << "---------------------- El Topo: Merging / topology change ----------------------" << std::endl;
   
   m_mesh.update_connectivity( m_positions.size() );
   
   rebuild_static_broad_phase();
   
   std::queue<Vec2ui> merge_candidates;
   
   //
   // Check point-edge proximities for zippering candidates
   //
   
   bool merge_occured = false;
   
   // sorted by proximity so we merge closest pairs first
   std::vector<SortablePointEdgeProximity> proximities;
   
   for(unsigned int i = 0; i < m_mesh.edges.size(); i++)
   {
      Vec2ui e0 = m_mesh.edges[i];
      
      if ( e0[0] == e0[1] ) { continue; }
      if ( m_masses[e0[0]] > 100 ) { continue; }
      
      Vec2d emin, emax;
      edge_static_bounds(i, emin, emax);
      emin -= m_merge_proximity_epsilon * Vec2d(1,1);
      emax += m_merge_proximity_epsilon * Vec2d(1,1);
      
      std::vector<unsigned int> point_candidates;
      m_broad_phase->get_potential_vertex_collisions( emin, emax, point_candidates );
      
      for(unsigned int j = 0; j < point_candidates.size(); j++)
      {
         unsigned int proximal_vertex_index = point_candidates[j];
         
         if ( m_masses[e0[0]] != m_masses[proximal_vertex_index] )
         {
            continue;
         }
         
         if( e0[0] != proximal_vertex_index && e0[1] != proximal_vertex_index )
         {
            double distance;
            Vec3d normal;
            
            point_segment_distance( false, 
                                   m_positions[proximal_vertex_index], proximal_vertex_index,
                                   m_positions[e0[0]], e0[0],
                                   m_positions[e0[1]], e0[1],
                                   distance );
            
            if ( distance < m_merge_proximity_epsilon )
            {
               proximities.push_back( SortablePointEdgeProximity( proximal_vertex_index, i, distance) );
            }
         }
      }
   }
   
   sort( proximities.begin(), proximities.end() );
   
   for ( unsigned int i = 0; i < proximities.size(); ++i )
   {
      unsigned int point_index = proximities[i].point_index;
      unsigned int edge_index = proximities[i].edge_index;
      
      if ( ( m_mesh.edges[edge_index][0] == point_index ) || ( m_mesh.edges[edge_index][1] == point_index ) ) 
      {
         continue;
      }
      
      if ( m_mesh.edges[edge_index][0] == m_mesh.edges[edge_index][1] )
      {
         continue;
      }
      
      if ( m_verbose ) { std::cout << "proximity: " << proximities[i].distance << std::endl; }
      
      if ( merge( point_index, edge_index ) )
      {
         if ( m_verbose ) 
         { 
            std::cout << "zippered vertex " << point_index << " and edge " << edge_index << std::endl; 
         }
         merge_occured = true;
         trim_non_manifold();
      }
   }   
}

   
// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

void SurfTrack::improve_mesh( )
{
   bool operation_occurred = true;
   unsigned int passes = 0;
   static unsigned int MAX_PASSES = 5;
   while ( m_perform_improvement && operation_occurred && passes++ < MAX_PASSES )
   {
      operation_occurred = false;
      operation_occurred |= split_pass();
      m_mesh.check_orientation();
      operation_occurred |= collapse_pass();      
      m_mesh.check_orientation();
   }
   
   partition_volumes();   
}

// ---------------------------------------------------------
/// 
///
///
// ---------------------------------------------------------

void SurfTrack::topology_changes( )
{
   if ( m_allow_topology_changes )
   {
      merge_pass();
      m_mesh.check_orientation();
      partition_volumes();
   }
   
}

   
   
}  // namespace


