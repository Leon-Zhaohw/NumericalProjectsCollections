// ---------------------------------------------------------
//
//  dynamicsurface.cpp
//  
//  An edge mesh with associated vertex locations and 
//  velocities.  Collision detection and solving.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include "dynamicsurface.h"
#include "broadphasegrid.h"

#include <vector>
#include <deque>
#include <queue>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else

#ifdef WIN32
#include <windows.h>
#endif

#include <GL/gl.h>

#endif

#include <ctime>

#include "vec.h"
#include "mat.h"

#include "ccd_wrapper.h"
#include "gluvi.h"

#include "edgemesh.h"

#include "wallclocktime.h"

#include <cassert>

#include "sparse_matrix.h"
#include "krylov_solvers.h"

#include <fstream>

//#include "lapack_wrapper.h"

#include "array2.h"


// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Extern globals
// ---------------------------------------------------------

namespace eltopo2d
{
   
// ---------------------------------------------------------
// Static function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Add a collision to the list as long as it doesn't have the same vertices as any other collisions in the list.
///
// ---------------------------------------------------------

static void add_unique_collision( std::vector<Collision>& collisions, const Collision& new_collision )
{
   for ( std::vector<Collision>::iterator iter = collisions.begin(); iter != collisions.end(); ++iter )
   {
      if ( iter->same_vertices( new_collision ) )
      {
         return;
      }
   }
   
   collisions.push_back( new_collision );
}

// ---------------------------------------------------------
///
/// Helper function: multiply transpose(A) * D * B
///
// ---------------------------------------------------------

static void AtDB(const SparseMatrixDynamicCSR &A, const double* diagD, const SparseMatrixDynamicCSR &B, SparseMatrixDynamicCSR &C)
{
   assert(A.m==B.m);
   C.resize(A.n, B.n);
   C.set_zero();
   for(int k=0; k<A.m; ++k)
   {
      const DynamicSparseVector& r = A.row[k];
      
      for( DynamicSparseVector::const_iterator p=r.begin(); p != r.end(); ++p )
      {
         int i = p->index;
         double multiplier = p->value * diagD[k];
         C.add_sparse_row( i, B.row[k], multiplier );
      }
   }
}

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

DynamicSurface::DynamicSurface( unsigned int num_vertices,
                               const double* vertices,
                               unsigned int num_edges,
                               const int* edges,
                               const double* masses,                  
                               double in_proximity_epsilon,
                               bool in_collision_safety,
                               bool in_verbose ) :
   
   m_proximity_epsilon( in_proximity_epsilon ),
   m_verbose( in_verbose ),
   m_collision_safety( in_collision_safety ),
   m_positions(0), m_newpositions(0), m_velocities(0), m_masses(0),
   m_volume_ids(0), m_volumes(0),
   m_mesh()
   
{
   
   if ( m_verbose ) { std::cout << "constructing dynamic surface" << std::endl; }
   
   for(unsigned int i = 0; i < num_vertices; i++)
   {
      m_positions.push_back( Vec2d( vertices[2*i], vertices[2*i+1] ) );
		m_newpositions.push_back( Vec2d( vertices[2*i], vertices[2*i+1] ) );
      m_velocities.push_back( Vec2d(0,0) );
      m_masses.push_back( masses[i] );
      m_volume_ids.push_back(0);
   }
   
   for(unsigned int i = 0; i < num_edges; i++)
   {
      m_mesh.edges.push_back( Vec2ui( edges[2*i], edges[2*i+1] ) );
   }
   
   m_mesh.update_connectivity( m_positions.size() );
   
   m_broad_phase = new BroadPhaseGrid();
   
   if ( m_verbose ) { std::cout << "constructed dynamic surface" << std::endl; }
   
   
}
   
// ---------------------------------------------------------
///
/// DynamicSurface constructor.  Copy edges and vertex locations.
///
// ---------------------------------------------------------

DynamicSurface::DynamicSurface( const std::vector<Vec2d>& vertex_positions, 
                                const std::vector<Vec2ui>& edges,
                                const std::vector<double>& masses,
                                double in_proximity_epsilon,
                                bool in_collision_safety,
                                bool in_verbose ) :
   m_proximity_epsilon( in_proximity_epsilon ),
   m_verbose( in_verbose ),
   m_collision_safety( in_collision_safety ),
   m_positions(0), m_newpositions(0), m_velocities(0), m_masses(0),
   m_volume_ids(0), m_volumes(0),
   m_mesh()
{
   
   if ( m_verbose ) { std::cout << "constructing dynamic surface" << std::endl; }
  
   for(unsigned int i = 0; i < vertex_positions.size(); i++)
   {
      m_positions.push_back( vertex_positions[i] );
		m_newpositions.push_back( vertex_positions[i] );
      m_velocities.push_back( Vec2d(0,0) );
      m_masses.push_back( masses[i] );
      m_volume_ids.push_back(0);
   }
   
   for(unsigned int i = 0; i < edges.size(); i++)
   {
      m_mesh.edges.push_back( edges[i] );
   }
   
   if ( m_verbose ) 
   {
      std::cout << "vs: " << m_positions.size() << std::endl;
      std::cout << "edges: " << m_mesh.edges.size() << std::endl;
   }
   
   m_mesh.update_connectivity( m_positions.size() );
   
   m_broad_phase = new BroadPhaseGrid();
   
   if ( m_verbose ) { std::cout << "constructed dynamic surface" << std::endl; }
   
}



// ---------------------------------------------------------
///
/// Add a triangle to the surface.  Update the underlying EdgeMesh and acceleration grid. 
///
// ---------------------------------------------------------

unsigned int DynamicSurface::add_edge( const Vec2ui& e )
{
   unsigned int new_edge_index = m_mesh.edges.size();
   m_mesh.nondestructive_add_edge( e[0], e[1] );
   
   return new_edge_index;
}


// ---------------------------------------------------------
///
/// Remove a triangle from the surface.  Update the underlying EdgeMesh and acceleration grid. 
///
// ---------------------------------------------------------

void DynamicSurface::remove_edge(unsigned int e)
{
   m_mesh.nondestructive_remove_edge( e );
}


// ---------------------------------------------------------
///
/// Add a vertex to the surface.  Update the acceleration grid. 
///
// ---------------------------------------------------------

unsigned int DynamicSurface::add_vertex( const Vec2d& new_vertex_position, 
                                         const Vec2d& new_vertex_velocity, 
                                         double new_vertex_mass, 
                                         unsigned int volume_id )
{
   m_positions.push_back( new_vertex_position );
   m_newpositions.push_back( new_vertex_position );
   m_velocities.push_back( new_vertex_velocity );
   m_masses.push_back( new_vertex_mass );
   m_volume_ids.push_back( volume_id );
   
   unsigned int new_vertex_index = m_mesh.nondestructive_add_vertex( );

   assert( new_vertex_index == m_positions.size() - 1 );
      
   return new_vertex_index;
}


// ---------------------------------------------------------
///
/// Remove a vertex from the surface.  Update the acceleration grid. 
///
// ---------------------------------------------------------

void DynamicSurface::remove_vertex( unsigned int vertex_index )
{
   m_mesh.nondestructive_remove_vertex( vertex_index );
   
   m_positions[ vertex_index ] = Vec2d( 0.0, 0.0 );
   m_newpositions[ vertex_index ] = Vec2d( 0.0, 0.0 );
   m_masses[vertex_index] = 0.0;
   m_volume_ids[ vertex_index ] = (unsigned int) ~0;
}

// --------------------------------------------------------
///
/// 
///
// --------------------------------------------------------
   
double DynamicSurface::get_curvature( unsigned int vertex_id ) const
{
   const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[vertex_id];
   if ( incident_edges.size() != 2 ) { return 0.0; }
   const Vec2ui& edge_a = m_mesh.edges[incident_edges[0]];
   const Vec2ui& edge_b = m_mesh.edges[incident_edges[1]];
   
   unsigned int prev, next;
   
   if ( edge_a[0] == vertex_id )
   {
      next = edge_a[1];
      assert( edge_b[1] == vertex_id );
      prev = edge_b[0];
   }  
   else
   {
      assert( edge_a[1] == vertex_id );
      prev = edge_a[0];
      assert( edge_b[0] == vertex_id );
      next = edge_b[1];
   }
   
   double a1 = 0.5 * (m_positions[next][0] - m_positions[prev][0]);
   double a2 = 0.5 * (m_positions[next][0] + m_positions[prev][0]) - m_positions[vertex_id][0];
   double b1 = 0.5 * (m_positions[next][1] - m_positions[prev][1]);
   double b2 = 0.5 * (m_positions[next][1] + m_positions[prev][1]) - m_positions[vertex_id][1];
   double ss = a1*a1 + b2*b2;
   
   return 2*(a1*b2 - a2*b1) / sqrt( ss*ss*ss );
   
}
   
   
// --------------------------------------------------------
///
/// 
///
// --------------------------------------------------------
   
double DynamicSurface::get_distance_to_surface( const Vec2d& x ) const
{

   Vec2d query_low = x - Vec2d(m_proximity_epsilon);
   Vec2d query_high = x + Vec2d(m_proximity_epsilon);
      
   std::vector<unsigned int> overlapping_edges;
   m_broad_phase->get_potential_edge_collisions( query_low, query_high, overlapping_edges );

   while ( overlapping_edges.empty() )
   {
      query_low -= Vec2d(0.1);
      query_high += Vec2d(0.1);
      m_broad_phase->get_potential_edge_collisions( query_low, query_high, overlapping_edges );
   }
   
   double min_distance = 1e30;
   
   for ( unsigned int j = 0; j < overlapping_edges.size(); ++j )
   {
      const Vec2ui& edge = m_mesh.edges[overlapping_edges[j]];
      point_segment_distance( true,
                              x, m_positions.size(),
                              m_positions[edge[0]], edge[0],
                              m_positions[edge[1]], edge[1],
                              min_distance );                             
   }
   
   return min_distance;
      
}
   

// --------------------------------------------------------
///
/// Determine if the specified point is inside the volume enclosed by this surface.
///
// --------------------------------------------------------
   
bool DynamicSurface::point_is_inside_volume( const Vec2d& x ) const
{
   std::vector<double> ss;
   std::vector<unsigned int> edges;

   //
   // The point is inside if there are an odd number of hits on a ray cast from the point. 
   // We'll cast four rays for numerical robustness.
   //
   
   unsigned int inside_votes = 0;
      
   // shoot a ray in the positive-x direction
   Vec2d ray_end = x + Vec2d( 1e+3, 0 );
   get_segment_collisions( x, ray_end, ss, edges );
   if ( ss.size() % 2 == 1 ) { ++inside_votes; }
   
   // negative x
   ray_end = x - Vec2d( 1e+3, 0 );
   ss.clear();
   get_segment_collisions( x, ray_end, ss, edges );
   if ( ss.size() % 2 == 1 ) { ++inside_votes; }

   // positive y
   ray_end = x + Vec2d( 0, 1e+3 );
   ss.clear();
   get_segment_collisions( x, ray_end, ss, edges );
   if ( ss.size() % 2 == 1 ) { ++inside_votes; }

   // negative y
   ray_end = x - Vec2d( 0, 1e+3 );
   ss.clear();   
   get_segment_collisions( x, ray_end, ss, edges );
   if ( ss.size() % 2 == 1 ) { ++inside_votes; }
   
   return ( inside_votes > 2 );
}

   
// --------------------------------------------------------
///
/// Determine volume IDs for all vertices
///
// --------------------------------------------------------

void DynamicSurface::partition_volumes()
{
   
   m_mesh.nondestructive_clear_unused();
   m_mesh.update_connectivity( m_positions.size() );

   clear_deleted_vertices();
   m_mesh.update_connectivity( m_positions.size() );
   
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      assert( m_mesh.edges[i][0] != m_mesh.edges[i][1] );
   }
   
   static const unsigned int UNASSIGNED = (unsigned int) ~0;
   
   m_volumes.clear();
   
   m_volume_ids.clear();
   m_volume_ids.resize( m_positions.size(), UNASSIGNED );
   
   unsigned int curr_volume = 0;
   
   while ( true )
   { 
      unsigned int next_unassigned_vertex;
      for ( next_unassigned_vertex = 0; next_unassigned_vertex < m_volume_ids.size(); ++next_unassigned_vertex )
      {
         if ( m_mesh.vtxedge[next_unassigned_vertex].size() == 0 ) { continue; }
         
         if ( m_volume_ids[next_unassigned_vertex] == UNASSIGNED )
         {
            break;
         }
      }
      
      if ( next_unassigned_vertex == m_volume_ids.size() )
      {
         break;
      }
      
      std::queue<unsigned int> open;
      open.push( next_unassigned_vertex );
      
      std::vector<unsigned int> volume_vertices;
      
      while ( false == open.empty() )
      {
         unsigned int vertex_index = open.front();
         open.pop();
         
         if ( m_mesh.vtxedge[vertex_index].size() == 0 ) { continue; }
         
         if ( m_volume_ids[vertex_index] != UNASSIGNED )
         {
            assert( m_volume_ids[vertex_index] == curr_volume );
            continue;
         }
         
         m_volume_ids[vertex_index] = curr_volume;
         volume_vertices.push_back( vertex_index );
         
         const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[vertex_index];
         
         for( unsigned int i = 0; i < incident_edges.size(); ++i )
         {
            unsigned int adjacent_vertex = m_mesh.edges[ incident_edges[i] ][0];
            if ( adjacent_vertex == vertex_index ) { adjacent_vertex = m_mesh.edges[ incident_edges[i] ][1]; }
            
            if ( m_volume_ids[adjacent_vertex] == UNASSIGNED )
            {
               open.push( adjacent_vertex );
            }
            else
            {
               assert( m_volume_ids[adjacent_vertex] == curr_volume );
            }
            
         } 
      }
      
      m_volumes.push_back( volume_vertices );
      
      ++curr_volume;
      
   }
   
   //
   // assert all vertices are assigned and share volume IDs with their neighbours
   //
   
   for ( unsigned int i = 0; i < m_volume_ids.size(); ++i )
   {
      if ( m_mesh.vtxedge[i].size() == 0 ) { continue; }
      
      assert( m_volume_ids[i] != UNASSIGNED );
      
      const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[i];    
      for( unsigned int j = 0; j < incident_edges.size(); ++j )
      {
         unsigned int adjacent_vertex = m_mesh.edges[ incident_edges[j] ][0];
         if ( adjacent_vertex == i ) { adjacent_vertex = m_mesh.edges[ incident_edges[j] ][1]; }
         assert( m_volume_ids[adjacent_vertex] == m_volume_ids[i] );         
      } 
      
   }
   
   
   //
   // Now we have partitioned vertices into connected sets (surfaces).  We need to find out if any of these surfaces are part of
   // the same volumes (e.g. inner and outer surface or a "bubble").
   //
   
   for ( unsigned int i = 0; i < m_volumes.size(); ++i )
   {
      
   }
      
}


// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------
   
void DynamicSurface::compute_all_edge_normals( std::vector<Vec2d>& midpoints, std::vector<Vec2d>& normals )
{
   midpoints.resize( m_mesh.edges.size() );
   normals.resize( m_mesh.edges.size() );
   
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      unsigned int ea = m_mesh.edges[i][0];
      unsigned int eb = m_mesh.edges[i][1];
      midpoints[i] = 0.5 * (m_positions[ea] + m_positions[eb]);
      normals[i] = normalized( -perp( m_positions[eb] - m_positions[ea] ) );
   }
      
}
   
   
// ---------------------------------------------------------
///
/// Remove all vertices not incident on any triangles.
///
// ---------------------------------------------------------

void DynamicSurface::clear_deleted_vertices( )
{

   unsigned int j = 0;
   
   for ( unsigned int i = 0; i < m_positions.size(); ++i )
   {
      std::vector<unsigned int>& inc_edges = m_mesh.vtxedge[i];

      if ( inc_edges.size() != 0 )
      {
         m_positions[j] = m_positions[i];
         m_newpositions[j] = m_newpositions[i];
         m_velocities[j] = m_velocities[i];
         m_masses[j] = m_masses[i];
         m_volume_ids[j] = m_volume_ids[i];
         
         for ( unsigned int t = 0; t < inc_edges.size(); ++t )
         {
            Vec2ui& edge = m_mesh.edges[ inc_edges[t] ];
            
            if ( edge[0] == i ) { edge[0] = j; }
            if ( edge[1] == i ) { edge[1] = j; }
         }
         
         ++j;
      }

   }
      
   m_positions.resize(j);
   m_newpositions.resize(j);
   m_velocities.resize(j);
   m_masses.resize(j);
   m_volume_ids.resize(j);

}


// ---------------------------------------------------------
///
/// ss are 0-1 alpha values
///
// ---------------------------------------------------------
   
void DynamicSurface::get_segment_collisions( const Vec2d& segment_point_a, 
                                             const Vec2d& segment_point_b, 
                                             std::vector<double>& ss, 
                                             std::vector<unsigned int>& edges ) const
{
   Vec2d query_low, query_high;
   minmax( segment_point_a, segment_point_b, query_low, query_high );
   query_low -= Vec2d( m_proximity_epsilon );
   query_high += Vec2d( m_proximity_epsilon );
   
   std::vector<unsigned int> overlapping_edges;
   m_broad_phase->get_potential_edge_collisions( query_low, query_high, overlapping_edges );
   
   for ( unsigned int j = 0; j < overlapping_edges.size(); ++j )
   {
      const Vec2ui& edge = m_mesh.edges[overlapping_edges[j]];
            
      Vec2ui sorted_edge = edge;
      if ( sorted_edge[1] < sorted_edge[0] ) { swap( sorted_edge[0], sorted_edge[1] ); }
      
      double segment_alpha, edge_alpha;

      const Vec2d edge_point_a = m_positions[ sorted_edge[0] ];
      const Vec2d edge_point_b = m_positions[ sorted_edge[1] ];
      
      unsigned int index_a = m_positions.size();
      unsigned int index_b = m_positions.size() + 1;
      
      if ( segment_segment_intersection( segment_point_a, index_a,
                                        segment_point_b, index_b,
                                        edge_point_a, sorted_edge[0], 
                                        edge_point_b, sorted_edge[1], 
                                        segment_alpha, edge_alpha ) )
      {
         ss.push_back( 1.0 - segment_alpha );
         edges.push_back( overlapping_edges[j] );
      }
   }
   
   
}
   
// ---------------------------------------------------------
///
/// Apply an impulse between a point and an edge
///
// ---------------------------------------------------------

void DynamicSurface::apply_edge_point_impulse( const Vec2ui& edge, 
                                               unsigned int v,
                                               double edge_alpha, 
                                               Vec2d& direction, 
                                               double magnitude )
{

   Vec2d& v0 = m_velocities[v];
   Vec2d& v1 = m_velocities[edge[0]];
   Vec2d& v2 = m_velocities[edge[1]];
   
   double inv_m0 = m_masses[v] < 100 ? 1.0 : 0.0;
   double inv_m1 = m_masses[edge[0]] < 100 ? 1.0: 0.0;
   double inv_m2 = m_masses[edge[1]] < 100 ? 1.0 : 0.0;

   double J = magnitude / (inv_m0 + edge_alpha*edge_alpha*inv_m1 + (1-edge_alpha)*(1-edge_alpha)*inv_m2);

   v0 += (J*inv_m0) * direction;
   v1 -= (J*edge_alpha*inv_m1) * direction;
   v2 -= (J*(1-edge_alpha)*inv_m2) * direction;
}
 


// ---------------------------------------------------------
///
/// Detect all edge-point proximities and apply repulsion impulses
///
// ---------------------------------------------------------

void DynamicSurface::handle_edge_point_proximities( double dt )
{
   std::cout << "---------------------- El Topo: Proximities ----------------------" << std::endl;
   
   rebuild_static_broad_phase();
   
   unsigned int point_edge_proximities = 0;
   
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      const Vec2ui& edge = m_mesh.edges[i];
      
      if ( edge[0] == edge[1] )    { continue; }

      for ( unsigned int v = 0; v < m_positions.size(); ++v )
      {
         if(edge[0] != v && edge[1] != v)
         {
            double distance, edge_alpha;
            Vec2d normal;
            
            point_segment_distance( false, m_positions[v], v, m_positions[edge[0]], edge[0], m_positions[edge[1]], edge[1],
                                    distance, edge_alpha, normal, 1.0 );
            
            if(distance < m_proximity_epsilon)
            {
               
               double relvel = dot(normal, m_velocities[v] - (edge_alpha*m_velocities[edge[0]] + (1-edge_alpha)*m_velocities[edge[1]]));
               double desired_relative_velocity = ( m_proximity_epsilon - distance ) / dt;
               double impulse = (desired_relative_velocity - relvel);
               
               apply_edge_point_impulse( edge, v, edge_alpha, normal, impulse);
               
               ++point_edge_proximities;
               
            }
         }      
      }
   }
   
}


void DynamicSurface::add_point_candidates(unsigned int v, CollisionCandidateSet& collision_candidates)
{
   Vec2d vmin, vmax;
   vertex_continuous_bounds(v, vmin, vmax);
   vmin -= Vec2d( m_proximity_epsilon );
   vmax += Vec2d( m_proximity_epsilon );
   
   std::vector<unsigned int> candidate_edges;
   m_broad_phase->get_potential_edge_collisions( vmin, vmax, candidate_edges);
   
   for(unsigned int j = 0; j < candidate_edges.size(); j++)
   {
      add_to_collision_candidates( collision_candidates, Vec2ui( candidate_edges[j], v ) );
   }
}
   
void DynamicSurface::add_edge_candidates(unsigned int e, CollisionCandidateSet& collision_candidates)
{
   Vec2d vmin, vmax;
   edge_continuous_bounds(e, vmin, vmax);
   vmin -= Vec2d( m_proximity_epsilon );
   vmax += Vec2d( m_proximity_epsilon );
   
   std::vector<unsigned int> candidate_vertices;
   m_broad_phase->get_potential_vertex_collisions( vmin, vmax, candidate_vertices);
   
   for(unsigned int j = 0; j < candidate_vertices.size(); j++)
   {
      add_to_collision_candidates( collision_candidates, Vec2ui( e, candidate_vertices[j] ) );
   }
   
}
   
void DynamicSurface::add_point_update_candidates(unsigned int v, CollisionCandidateSet& collision_candidates)
{
   add_point_candidates(v, collision_candidates);
   
   std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[v];
   
   for(unsigned int i = 0; i < incident_edges.size(); i++)
      add_edge_candidates(incident_edges[i], collision_candidates);   
}
   
   
// ---------------------------------------------------------
///
/// Perform one sweep of impulse collision handling, only for "deformable" vertices against "solid" edges
///
// ---------------------------------------------------------

void DynamicSurface::handle_point_vs_solid_edge_collisions( double dt )
{
   std::cout << "---------------------- El Topo: Deformable-solid point-edge collisions ----------------------" << std::endl;
   
   for(unsigned int i = 0; i < m_mesh.edges.size(); i++)
   {
      const Vec2ui& edge = m_mesh.edges[i];

      for(unsigned int v = 0; v < m_positions.size(); v++)
      {     

         if ( m_masses[v] < 100 && m_masses[edge[0]] > 100 && m_masses[edge[1]] > 100 )
         {
         
            if( edge[0] != v && edge[1] != v )
            {
               double time, edge_alpha, rel_disp;
               Vec2d normal;
                            
               Vec2ui sorted_edge = edge;
               if ( sorted_edge[1] < sorted_edge[0] ) { swap( sorted_edge[0], sorted_edge[1] ); }
               
               if ( point_segment_collision( m_positions[v], m_newpositions[v], v,
                                             m_positions[sorted_edge[0]], m_newpositions[sorted_edge[0]], sorted_edge[0],
                                             m_positions[sorted_edge[1]], m_newpositions[sorted_edge[1]], sorted_edge[1],
                                             edge_alpha,
                                             normal, time, rel_disp ) )                                 
               {
                  
                  double relvel = rel_disp / dt;
                  
                  apply_edge_point_impulse(edge, v, edge_alpha, normal, -relvel);
                                  
                  m_newpositions[v] = m_positions[v] + dt*m_velocities[v];
                  m_newpositions[edge[0]] = m_positions[edge[0]] + dt*m_velocities[edge[0]];
                  m_newpositions[edge[1]] = m_positions[edge[1]] + dt*m_velocities[edge[1]];
                  
               }
                  
            }
         }
      }      
   }
   
}


// ---------------------------------------------------------
///
/// Detect all continuous collisions and apply impulses to prevent them.
/// Return true if all collisions were resolved.
///
// ---------------------------------------------------------

bool DynamicSurface::handle_collisions(double dt)
{
   
   const unsigned int MAX_PASS = 3;
   const unsigned int MAX_CANDIDATES = (unsigned int) 1e+6; 
   
   CollisionCandidateSet update_collision_candidates;
   
   if ( MAX_PASS == 0 )
   {
      return false;
   }
   
   bool collision_found = true;
   bool candidate_overflow = false;

   for ( unsigned int pass = 0; ( collision_found && (pass < MAX_PASS) ); ++pass )
   {
      collision_found = false;
      
      std::cout << "---------------------- El Topo: Collisions ----------------------" << std::endl;
      
      rebuild_continuous_broad_phase();
      
      for(unsigned int i = 0; i < m_mesh.edges.size(); i++)
      {
         const Vec2ui& edge = m_mesh.edges[i];
         
         for(unsigned int v = 0; v < m_positions.size(); v++)
         {     
            
            if( edge[0] != v && edge[1] != v )
            {
               double time, edge_alpha, rel_disp;
               Vec2d normal;
               
               Vec2ui sorted_edge = edge;
               if ( sorted_edge[1] < sorted_edge[0] ) { swap( sorted_edge[0], sorted_edge[1] ); }
               
               if ( point_segment_collision( m_positions[v], m_newpositions[v], v,
                                            m_positions[sorted_edge[0]], m_newpositions[sorted_edge[0]], sorted_edge[0],
                                            m_positions[sorted_edge[1]], m_newpositions[sorted_edge[1]], sorted_edge[1],
                                            edge_alpha,
                                            normal, time, rel_disp ) )                                 
               {

                  if ( m_verbose )
                  {
                     std::cout << "collision normal: " << normal << std::endl;
                     std::cout << m_positions[v] << " " << m_newpositions[v] << std::endl;
                     std::cout << m_positions[sorted_edge[0]] << " " << m_newpositions[sorted_edge[0]] << std::endl;
                     std::cout << m_positions[sorted_edge[1]] << " " << m_newpositions[sorted_edge[1]] << std::endl;
                  }
                  
                  double relvel = rel_disp / dt;

                  if ( m_verbose )
                  {
                     std::cout << "impulse: " << -relvel << std::endl;
                     std::cout << "edge_alpha: " << edge_alpha << std::endl;                  
                  }
                  
                  apply_edge_point_impulse(sorted_edge, v, edge_alpha, normal, -relvel);
                  
                  m_newpositions[v] = m_positions[v] + dt*m_velocities[v];
                  m_newpositions[edge[0]] = m_positions[edge[0]] + dt*m_velocities[edge[0]];
                  m_newpositions[edge[1]] = m_positions[edge[1]] + dt*m_velocities[edge[1]];

                  if ( m_verbose )
                  {
                     std::cout << "post-correction: " << std::endl;
                     std::cout << m_positions[v] << " " << m_newpositions[v] << std::endl;
                     std::cout << m_positions[sorted_edge[0]] << " " << m_newpositions[sorted_edge[0]] << std::endl;
                     std::cout << m_positions[sorted_edge[1]] << " " << m_newpositions[sorted_edge[1]] << std::endl;
                  }
                  
                  if ( pass == MAX_PASS - 1 )
                  {
                     if ( update_collision_candidates.size() < MAX_CANDIDATES )
                     {
                        add_point_update_candidates(v, update_collision_candidates);
                        add_point_update_candidates(edge[0], update_collision_candidates);
                        add_point_update_candidates(edge[1], update_collision_candidates);
                     }
                     else
                     {
                        candidate_overflow = true;
                     }
                  }
                  
                  collision_found = true;
                  
               }            
            }
         }      
      }
                  
   }     // pass
   
   if ( m_verbose )
   {
      std::cout << "---------------------- unique-ifying update_collision_candidates ----------------------" << std::endl;
   }
   
   {
      CollisionCandidateSet::iterator new_end = std::unique( update_collision_candidates.begin(), update_collision_candidates.end() );
      update_collision_candidates.erase( new_end, update_collision_candidates.end() );
   }

   if ( m_verbose )
   {
      std::cout << "---------------------- handling new collisions ----------------------" << std::endl;
   }
   
   unsigned int n = update_collision_candidates.size();
   unsigned int c = 0;
   
   while( !update_collision_candidates.empty() && c < (5 * n) )
   {

      c++;
      CollisionCandidateSet::iterator iter = update_collision_candidates.begin();
      Vec2ui candidate = *iter;
      update_collision_candidates.erase(iter);
            
      unsigned int e = candidate[0];
      const Vec2ui& edge = m_mesh.edges[e];
      unsigned int v = candidate[1];
      
      if( edge[0] != v && edge[1] != v )
      {
         double time, edge_alpha, rel_disp;
         Vec2d normal;
         
         Vec2ui sorted_edge = edge;
         if ( sorted_edge[1] < sorted_edge[0] ) { swap( sorted_edge[0], sorted_edge[1] ); }
         
         if ( point_segment_collision( m_positions[v], m_newpositions[v], v,
                                      m_positions[sorted_edge[0]], m_newpositions[sorted_edge[0]], sorted_edge[0],
                                      m_positions[sorted_edge[1]], m_newpositions[sorted_edge[1]], sorted_edge[1],
                                      edge_alpha,
                                      normal, time, rel_disp ) )                                 
         {
            
            double relvel = rel_disp / dt;
            
            apply_edge_point_impulse(edge, v, edge_alpha, normal, -relvel);
            
            m_newpositions[v] = m_positions[v] + dt*m_velocities[v];
            m_newpositions[edge[0]] = m_positions[edge[0]] + dt*m_velocities[edge[0]];
            m_newpositions[edge[1]] = m_positions[edge[1]] + dt*m_velocities[edge[1]];
            
            if ( update_collision_candidates.size() < MAX_CANDIDATES )
            {
               add_point_update_candidates(v, update_collision_candidates);
               add_point_update_candidates(edge[0], update_collision_candidates);
               add_point_update_candidates(edge[1], update_collision_candidates);
            }
            else
            {
               candidate_overflow = true;
            }                           
         }
      }
      
   }
   
   return ( !candidate_overflow ) && ( update_collision_candidates.empty() );
   
}


// ---------------------------------------------------------
///
/// Detect all continuous collisions
///
// ---------------------------------------------------------

bool DynamicSurface::detect_collisions( std::vector<Collision>& collisions )
{
   if ( m_verbose )
   {
      std::cout << "---------------------- detecting collisions ----------------------" << std::endl;
   }
   
   static const unsigned int MAX_COLLISIONS = 5000;
        
   rebuild_continuous_broad_phase();

   for(unsigned int i = 0; i < m_mesh.edges.size(); i++)
   {
      const Vec2ui& edge = m_mesh.edges[i];
      
      for(unsigned int v = 0; v < m_positions.size(); v++)
      {     
         
         if( edge[0] != v && edge[1] != v )
         {
            double time, edge_alpha, rel_disp;
            Vec2d normal;
            
            Vec2ui sorted_edge = edge;
            if ( sorted_edge[1] < sorted_edge[0] ) { swap( sorted_edge[0], sorted_edge[1] ); }
            
            if ( point_segment_collision( m_positions[v], m_newpositions[v], v,
                                         m_positions[sorted_edge[0]], m_newpositions[sorted_edge[0]], sorted_edge[0],
                                         m_positions[sorted_edge[1]], m_newpositions[sorted_edge[1]], sorted_edge[1],
                                         edge_alpha,
                                         normal, time, rel_disp ) )                                 
            {
               
               Collision new_collision( Vec3ui( v, sorted_edge[0], sorted_edge[1] ), normal, edge_alpha );
               
               collisions.push_back( new_collision );              
               
               if ( collisions.size() > MAX_COLLISIONS ) 
               {
                  return false; 
               }
            }            
         }
      }      
   }

   
   return true;
   
}


// ---------------------------------------------------------
///
/// Detect continuous collisions among elements in the given ImpactZones, and adjacent to the given ImpactZones.
///
// ---------------------------------------------------------

void DynamicSurface::detect_new_collisions( const std::vector<ImpactZone> impact_zones, std::vector<Collision>& collisions )
{
   if ( m_verbose )
   {
      std::cout << "---------------------- detecting new collisions ----------------------" << std::endl;
   }
   
   std::vector<unsigned int> zone_vertices;
   std::vector<unsigned int> zone_edges;

   // Get all vertices in the impact zone
   
   for ( unsigned int i = 0; i < impact_zones.size(); ++i )
   {
      for ( unsigned int j = 0; j < impact_zones[i].collisions.size(); ++j )
      {
         add_unique( zone_vertices, impact_zones[i].collisions[j].vertex_indices[0] );
         add_unique( zone_vertices, impact_zones[i].collisions[j].vertex_indices[1] );
         add_unique( zone_vertices, impact_zones[i].collisions[j].vertex_indices[2] );
         
         //update_continuous_broad_phase( impact_zones[i].collisions[j].vertex_indices[0] );
         //update_continuous_broad_phase( impact_zones[i].collisions[j].vertex_indices[1] );
         //update_continuous_broad_phase( impact_zones[i].collisions[j].vertex_indices[2] );
         
      }
   }
   
   // Get all edges in the impact zone
   
   for ( unsigned int i = 0; i < zone_vertices.size(); ++i )
   {
      for ( unsigned int j = 0; j < m_mesh.vtxedge[zone_vertices[i]].size(); ++j )
      {
         add_unique( zone_edges, m_mesh.vtxedge[zone_vertices[i]][j] );
      }
   }

   // Check zone vertices vs all edges
   
   for ( unsigned int i = 0; i < zone_vertices.size(); ++i )
   {
      unsigned int vertex_index = zone_vertices[i];

      Vec2d query_low, query_high;
      vertex_continuous_bounds( zone_vertices[i], query_low, query_high );
      std::vector<unsigned int> overlapping_edges;
      m_broad_phase->get_potential_edge_collisions( query_low, query_high, overlapping_edges );
      
      for ( unsigned int j = 0; j < overlapping_edges.size(); ++j )
      {
         const Vec2ui& edge = m_mesh.edges[overlapping_edges[j]];
         
         if( edge[0] != vertex_index && edge[1] != vertex_index )
         {
            double time, edge_alpha, rel_disp;
            Vec2d normal;
            
            Vec2ui sorted_edge = edge;
            if ( sorted_edge[1] < sorted_edge[0] ) { swap( sorted_edge[0], sorted_edge[1] ); }
            
            if ( point_segment_collision( m_positions[vertex_index], m_newpositions[vertex_index], vertex_index,
                                          m_positions[sorted_edge[0]], m_newpositions[sorted_edge[0]], sorted_edge[0],
                                          m_positions[sorted_edge[1]], m_newpositions[sorted_edge[1]], sorted_edge[1],
                                          edge_alpha,
                                          normal, time, rel_disp ) )                                 
            {
               Collision new_collision( Vec3ui( vertex_index, sorted_edge[0], sorted_edge[1] ), normal, edge_alpha );
               add_unique_collision( collisions, new_collision );
            } 
         }
      }
   }
   
   // Check zone edges vs all vertices
   
   for ( unsigned int i = 0; i < zone_edges.size(); ++i )
   {
      const Vec2ui& edge = m_mesh.edges[zone_edges[i]];
      
      Vec2d query_low, query_high;
      edge_continuous_bounds( zone_edges[i], query_low, query_high );
      std::vector<unsigned int> overlapping_vertices;
      m_broad_phase->get_potential_vertex_collisions( query_low, query_high, overlapping_vertices );
      
      for ( unsigned int j = 0; j < overlapping_vertices.size(); ++j )
      {                      
         unsigned int vertex_index = overlapping_vertices[j];
         
         assert( m_mesh.vtxedge[vertex_index].size() != 0 );
         
         if( edge[0] != vertex_index && edge[1] != vertex_index )
         {
            double time, edge_alpha, rel_disp;
            Vec2d normal;
            
            Vec2ui sorted_edge = edge;
            if ( sorted_edge[1] < sorted_edge[0] ) { swap( sorted_edge[0], sorted_edge[1] ); }
            
            if ( point_segment_collision( m_positions[vertex_index], m_newpositions[vertex_index], vertex_index,
                                         m_positions[sorted_edge[0]], m_newpositions[sorted_edge[0]], sorted_edge[0],
                                         m_positions[sorted_edge[1]], m_newpositions[sorted_edge[1]], sorted_edge[1],
                                         edge_alpha,
                                         normal, time, rel_disp ) )                                 
            {
               Collision new_collision( Vec3ui( vertex_index, sorted_edge[0], sorted_edge[1] ), normal, edge_alpha );
               add_unique_collision( collisions, new_collision );
            } 
            
         }
      }
   }
   
}



// ---------------------------------------------------------
///
/// Combine impact zones which have overlapping vertex stencils
///
// ---------------------------------------------------------

void DynamicSurface::merge_impact_zones( std::vector<ImpactZone>& new_impact_zones, std::vector<ImpactZone>& master_impact_zones )
{

   bool merge_ocurred = true;
      
   for ( unsigned int i = 0; i < master_impact_zones.size(); ++i )
   {
      master_impact_zones[i].all_solved = true;
   }

   for ( unsigned int i = 0; i < new_impact_zones.size(); ++i )
   {
      new_impact_zones[i].all_solved = false;
   }

   
   while ( merge_ocurred )
   {
            
      merge_ocurred = false;
           
      for ( unsigned int i = 0; i < new_impact_zones.size(); ++i )
      {
         bool i_is_disjoint = true;
         
         for ( unsigned int j = 0; j < master_impact_zones.size(); ++j )
         {
            if ( master_impact_zones[j].share_vertices( new_impact_zones[i] ) )
            {
               
               master_impact_zones[j].all_solved &= new_impact_zones[i].all_solved;
               
               // steal all of i's collisions
               for ( unsigned int c = 0; c < new_impact_zones[i].collisions.size(); ++c )
               {

                  bool same_collision_exists = false;
                  
                  for ( unsigned int m = 0; m < master_impact_zones[j].collisions.size(); ++m )
                  {                 
                     if ( master_impact_zones[j].collisions[m].same_vertices( new_impact_zones[i].collisions[c] ) )
                     {

                        same_collision_exists = true;
                        break;
                     }
                     
                  }
                  
                  
                  master_impact_zones[j].collisions.push_back( new_impact_zones[i].collisions[c] );
                  
               }
               
               i_is_disjoint = false;
               merge_ocurred = true;
               break;
            }
         }     // end for(j)
         
         if ( i_is_disjoint )
         {
            // copy the impact zone
            
            ImpactZone new_zone;
            for ( unsigned int c = 0; c < new_impact_zones[i].collisions.size(); ++c )
            {
               new_zone.collisions.push_back( new_impact_zones[i].collisions[c] );
            }
            
            new_zone.all_solved = new_impact_zones[i].all_solved;
            
            master_impact_zones.push_back( new_zone );
         }
      }     // end for(i)
      
      new_impact_zones = master_impact_zones;
      master_impact_zones.clear();
      
   }  // while

   master_impact_zones = new_impact_zones;
     
}


   
// ---------------------------------------------------------
///
/// Iteratively project out relative normal velocities for a set of collisions in an impact zone until all collisions are solved.
///
// ---------------------------------------------------------

bool DynamicSurface::iterated_inelastic_projection( ImpactZone& iz, double dt )
{
   assert( m_masses.size() == m_positions.size() );
   
   static const unsigned int MAX_PROJECTION_ITERATIONS = 10;
   
   for ( unsigned int i = 0; i < MAX_PROJECTION_ITERATIONS; ++i )
   {
      //m_verbose = true;
      bool solver_ok = inelastic_projection( iz );
      m_verbose = false;

      if ( !solver_ok )
      {
         return false;
      }
      
      bool collision_still_exists = false;
      
      for ( unsigned int c = 0; c < iz.collisions.size(); ++c )
      {
         // run collision detection on this pair again
         
         Collision& collision = iz.collisions[c];
         const Vec3ui& vs = collision.vertex_indices;
            
         m_newpositions[vs[0]] = m_positions[vs[0]] + dt * m_velocities[vs[0]];
         m_newpositions[vs[1]] = m_positions[vs[1]] + dt * m_velocities[vs[1]];
         m_newpositions[vs[2]] = m_positions[vs[2]] + dt * m_velocities[vs[2]];
            
         double time, edge_alpha, rel_disp;
         Vec2d normal;

         assert( vs[1] < vs[2] );
         
         if ( point_segment_collision( m_positions[vs[0]], m_newpositions[vs[0]], vs[0],
                                       m_positions[vs[1]], m_newpositions[vs[1]], vs[1],
                                       m_positions[vs[2]], m_newpositions[vs[2]], vs[2],
                                       edge_alpha,
                                       normal, time, rel_disp ) )                                 
         {
            
            collision.normal = normal;
            collision.edge_alpha = edge_alpha;
            collision_still_exists = true;            
         }
         
      } // for collisions
      
      if ( false == collision_still_exists )  
      {
         return true; 
      }
      
   } // for iterations
 
   if ( m_verbose ) 
   {
      std::cout << "reached max iterations for this zone" << std::endl;
   }
   
   return false;
   
}


// ---------------------------------------------------------
///
/// Project out relative normal velocities for a set of collisions in an impact zone.
///
// ---------------------------------------------------------

bool DynamicSurface::inelastic_projection( const ImpactZone& iz )
{
   
   std::cout << " ----- inelastic projection ----- " << std::endl;
   
   static double IMPULSE_MULTIPLIER = 1.0;
   
   const unsigned int k = iz.collisions.size();    // notation from [Harmon et al 2008]: k == number of collisions
   
   std::vector<unsigned int> zone_vertices;
   iz.get_all_vertices( zone_vertices );
   
   const unsigned int n = zone_vertices.size();       // n == number of distinct colliding vertices
   
   if ( m_verbose )
   {
      std::cout << "GCT: " << 2*n << "x" << k << std::endl;
   }
   
   SparseMatrixDynamicCSR GCT( 2*n, k );
   GCT.set_zero();

   // construct matrix grad C transpose
   for ( unsigned int i = 0; i < k; ++i )
   {
      // set col i
      const Collision& coll = iz.collisions[i];

      Vec3d alphas( 1.0, -coll.edge_alpha, -(1.0 - coll.edge_alpha) );
      
      for ( unsigned int v = 0; v < 3; ++v )
      {
         // block row j ( == block column j of grad C )
         unsigned int j = coll.vertex_indices[v];
         
         std::vector<unsigned int>::iterator zone_vertex_iter = find( zone_vertices.begin(), zone_vertices.end(), j );
         
         assert( zone_vertex_iter != zone_vertices.end() );
         
         unsigned int mat_j = (unsigned int) ( zone_vertex_iter - zone_vertices.begin() );
         
         GCT(mat_j*2, i)   = alphas[v] * coll.normal[0];
         GCT(mat_j*2+1, i) = alphas[v] * coll.normal[1];
         
      }
   }

   Array1d inv_masses;
   inv_masses.reserve(2*n);
   Array1d column_velocities;
   column_velocities.reserve(2*n);
   
   for ( unsigned int i = 0; i < n; ++i )
   {
      if ( m_masses[zone_vertices[i]] < 100.0 )
      {
         inv_masses.push_back( 1.0 / m_masses[zone_vertices[i]] );
         inv_masses.push_back( 1.0 / m_masses[zone_vertices[i]] );
      }
      else
      {
         // infinite mass (scripted objects)         
         inv_masses.push_back( 0.0 );
         inv_masses.push_back( 0.0 );         
      }

      column_velocities.push_back( m_velocities[zone_vertices[i]][0] );
      column_velocities.push_back( m_velocities[zone_vertices[i]][1] );
   }

   //
   // minimize | M^(-1/2) * GC^T x - M^(1/2) * v |^2
   //
   
   // solution vector
   Array1d x(k);
         
   KrylovSolverStatus solver_result;
   
   // normal equations: GC * M^(-1) GCT * x = GC * v
   //                   A * x = b

   SparseMatrixDynamicCSR A( k, k );
   A.set_zero();
   AtDB( GCT, inv_masses.data, GCT, A ); 
   
   Array1d b(k);
   GCT.apply_transpose( column_velocities.data, b.data );   

   if ( m_verbose )  
   { 
      std::cout << "system built" << std::endl; 
      std::cout << "A: " << std::endl;
      for ( int i = 0; i < A.m; ++i )
      {
         for ( int j = 0; j < A.n; ++j )
         {
            std::cout << A(i,j) << " ";
         }
         std::cout << std::endl;
      }
      std::cout << "b: " << std::endl;
      for ( unsigned int i = 0; i < b.size(); ++i )
      {
         std::cout << b[i] << std::endl;
      }
   }
   
   MINRES_CR_Solver solver;   
   SparseMatrixStaticCSR solver_matrix( A );    // convert dynamic to static
   solver.max_iterations = 100;
   solver_result = solver.solve( solver_matrix, b.data, x.data ); 

   if ( solver_result != KRYLOV_CONVERGED )
   {
      std::cout << "CR solver failed: ";
      
      if ( solver_result == KRYLOV_BREAKDOWN )
      {
         std::cout << "KRYLOV_BREAKDOWN" << std::endl;
      }
      else
      {
         std::cout << "KRYLOV_EXCEEDED_MAX_ITERATIONS" << std::endl;
      }
      
      return false;          
   } 
   
   // apply impulses 
   Array1d applied_impulses(2*n);
   GCT.apply( x.data, applied_impulses.data );
     
   for ( unsigned int i = 0; i < applied_impulses.size(); ++i )
   {
      column_velocities[i] -= IMPULSE_MULTIPLIER * inv_masses[i] * applied_impulses[i];
   }
   
   for ( unsigned int i = 0; i < n; ++i )
   {
      Vec2d new_velocity( column_velocities[2*i], column_velocities[2*i + 1] );
      
      m_velocities[zone_vertices[i]][0] = column_velocities[2*i];
      m_velocities[zone_vertices[i]][1] = column_velocities[2*i + 1];
   }
   
   return true;
   
}


// ---------------------------------------------------------
///
/// Handle all collisions simultaneously by iteratively solving individual impact zones until no new collisions are detected.
///
// ---------------------------------------------------------

bool DynamicSurface::handle_collisions_simultaneous(double dt)
{
   if ( m_verbose )
   {
      std::cout << "---------------------- handling simultaneous collisions ----------------------" << std::endl;
   }

   // copy
   std::vector<Vec2d> old_velocities = m_velocities;

   std::vector<ImpactZone> impact_zones;

   bool finished_detecting_collisions = false;
   
   std::vector<Collision> total_collisions;
   finished_detecting_collisions = detect_collisions(total_collisions);
   
   while ( false == total_collisions.empty() )
   {      
      // insert each new collision constraint into its own impact zone
      std::vector<ImpactZone> new_impact_zones;
      for ( unsigned int i = 0; i < total_collisions.size(); ++i )
      {
         ImpactZone new_zone;
         new_zone.collisions.push_back( total_collisions[i] );
         new_impact_zones.push_back( new_zone );
      }
      
      assert( new_impact_zones.size() == total_collisions.size() );

      // merge all impact zones that share vertices
      merge_impact_zones( new_impact_zones, impact_zones );

      for ( int i = 0; i < (int) impact_zones.size(); ++i )
      {
         if ( impact_zones[i].all_solved ) 
         {
            impact_zones.erase( impact_zones.begin() + i );
            --i;
         }
      }

      unsigned int total_num_collisions = 0;
      for ( int i = 0; i < (int) impact_zones.size(); ++i )
      {
         assert( false == impact_zones[i].all_solved );
         total_num_collisions += impact_zones[i].collisions.size();
      }            
      
      if ( m_verbose )
      {
         std::cout << "---------------------- solving impact zones ----------------------" << std::endl;
      }
      
      // for each impact zone
      for ( unsigned int i = 0; i < impact_zones.size(); ++i )
      {
         
         std::vector<unsigned int> zone_vertices;
         
         // reset impact zone to pre-response m_velocities
         for ( unsigned int j = 0; j < impact_zones[i].collisions.size(); ++j )
         {
            Vec3ui& vs = impact_zones[i].collisions[j].vertex_indices;
            
            m_velocities[vs[0]] = old_velocities[vs[0]];
            m_velocities[vs[1]] = old_velocities[vs[1]];
            m_velocities[vs[2]] = old_velocities[vs[2]];
            
            zone_vertices.push_back( vs[0] );
            zone_vertices.push_back( vs[1] );
            zone_vertices.push_back( vs[2] );
            
         }
         
         // apply inelastic projection
         
         bool solver_ok = iterated_inelastic_projection( impact_zones[i], dt );
         
         if ( false == solver_ok )
         {
            if ( m_verbose )
            {
               std::cout << "Impoact zone solver problem.  Returning." << std::endl;
            }
            return false;
         }
         
         // reset predicted m_positions
         for ( unsigned int j = 0; j < impact_zones[i].collisions.size(); ++j )
         {
            const Vec3ui& vs = impact_zones[i].collisions[j].vertex_indices;            

            m_newpositions[vs[0]] = m_positions[vs[0]] + dt * m_velocities[vs[0]];
            m_newpositions[vs[1]] = m_positions[vs[1]] + dt * m_velocities[vs[1]];
            m_newpositions[vs[2]] = m_positions[vs[2]] + dt * m_velocities[vs[2]];
         } 
      
      }  // for IZs
      
      total_collisions.clear();
      
      if ( !finished_detecting_collisions )
      {
         if ( m_verbose )
         {
            std::cout << "attempting to finish global collision detection" << std::endl;
         }
         finished_detecting_collisions = detect_collisions( total_collisions );
         impact_zones.clear();
      }
      else
      {
         detect_new_collisions( impact_zones, total_collisions );
      }
   
   }
   
   return true;
}


   
   
bool DynamicSurface::compute_rigid_motion( const ImpactZone& iz )
{
   Vec2d centre_of_mass(0,0);
   Vec2d avg_velocity(0,0);
   double total_mass = 0.0;
   
   std::vector<unsigned int> vertices;
   iz.get_all_vertices( vertices );
   
   for ( unsigned int i = 0; i < vertices.size(); ++i )
   {
      unsigned int v = vertices[i];
      centre_of_mass += m_masses[v] * m_positions[v];
      avg_velocity += m_masses[v] * m_velocities[v];
      total_mass += m_masses[v];
   }
   centre_of_mass /= total_mass;
   avg_velocity /= total_mass;
   
   double angular_momentum = 0.0;
   double moment_of_inertia = 0.0;
   for ( unsigned int i = 0; i < vertices.size(); ++i )
   {
      unsigned int v = vertices[i];
      angular_momentum += m_masses[v] * cross( m_positions[v] - centre_of_mass, m_velocities[v] - avg_velocity );
      moment_of_inertia += m_masses[v] * mag2( m_positions[v] - centre_of_mass );
   }
   
   if ( fabs(moment_of_inertia) < 1e-10 )
   {
      return false;
   }
   
   double new_angular_velocity = 1.0 / moment_of_inertia * angular_momentum;

   for ( unsigned int i = 0; i < vertices.size(); ++i )
   {
      unsigned int v = vertices[i];
      m_velocities[v] = avg_velocity + new_angular_velocity * perp( m_positions[v] - centre_of_mass );
   }
   
   return true;
   
}

   
bool DynamicSurface::compute_average_motion( const ImpactZone& iz )
{
   Vec2d centre_of_mass(0,0);
   Vec2d avg_velocity(0,0);
   double total_mass = 0.0;
   
   std::vector<unsigned int> vertices;
   iz.get_all_vertices( vertices );
   
   for ( unsigned int i = 0; i < vertices.size(); ++i )
   {
      unsigned int v = vertices[i];
      avg_velocity += m_masses[v] * m_velocities[v];
      total_mass += m_masses[v];
   }
   
   avg_velocity /= total_mass;

   for ( unsigned int i = 0; i < vertices.size(); ++i )
   {
      unsigned int v = vertices[i];
      m_velocities[v] = avg_velocity;
   }
   
   return true;
   
}
   
   
   
bool DynamicSurface::rigid_impact_zones( double dt )
{
   if ( m_verbose )
   {
      std::cout << "---------------------- rigid impact zones ----------------------" << std::endl;
   }
   
   std::vector<ImpactZone> impact_zones;
   
   bool finished_detecting_collisions = false;
   
   std::vector<Collision> total_collisions;
   finished_detecting_collisions = detect_collisions(total_collisions);
   
   if ( m_verbose ) { std::cout << total_collisions.size() << " collisions detected " << std::endl; }
   
   while ( false == total_collisions.empty() )
   {      
      // insert each new collision constraint into its own impact zone
      std::vector<ImpactZone> new_impact_zones;
      for ( unsigned int i = 0; i < total_collisions.size(); ++i )
      {
         ImpactZone new_zone;
         new_zone.collisions.push_back( total_collisions[i] );
         new_impact_zones.push_back( new_zone );
      }
      
      assert( new_impact_zones.size() == total_collisions.size() );
      
      // merge all impact zones that share vertices
      merge_impact_zones( new_impact_zones, impact_zones );
      
      for ( int i = 0; i < (int) impact_zones.size(); ++i )
      {
         if ( impact_zones[i].all_solved ) 
         {
            impact_zones.erase( impact_zones.begin() + i );
            --i;
         }
      }
      
      unsigned int total_num_collisions = 0;
      for ( int i = 0; i < (int) impact_zones.size(); ++i )
      {
         assert( false == impact_zones[i].all_solved );
         total_num_collisions += impact_zones[i].collisions.size();
      }            
      
      if ( m_verbose )
      {
         std::cout << impact_zones.size() << " impact zones built.  Solving with RIZ solver..." << std::endl;
      }
      
      // for each impact zone
      for ( unsigned int i = 0; i < impact_zones.size(); ++i )
      {
         // solve for rigid motion
         
         bool solver_ok = compute_rigid_motion( impact_zones[i] );
         
         if ( false == solver_ok )
         {
            std::cout << "WARNING: rigid impact zone solver problem" << std::endl;
            return false;
         }
         
         // update predicted m_positions
         for ( unsigned int j = 0; j < impact_zones[i].collisions.size(); ++j )
         {
            const Vec3ui& vs = impact_zones[i].collisions[j].vertex_indices;            
            m_newpositions[vs[0]] = m_positions[vs[0]] + dt * m_velocities[vs[0]];
            m_newpositions[vs[1]] = m_positions[vs[1]] + dt * m_velocities[vs[1]];
            m_newpositions[vs[2]] = m_positions[vs[2]] + dt * m_velocities[vs[2]];
         } 
         
         
         bool collision_still_exists = false;
         
         for ( unsigned int c = 0; c < impact_zones[i].collisions.size(); ++c )
         {
            // run collision detection on this pair again
            
            Collision& collision = impact_zones[i].collisions[c];
            const Vec3ui& vs = collision.vertex_indices;
            
            m_newpositions[vs[0]] = m_positions[vs[0]] + dt * m_velocities[vs[0]];
            m_newpositions[vs[1]] = m_positions[vs[1]] + dt * m_velocities[vs[1]];
            m_newpositions[vs[2]] = m_positions[vs[2]] + dt * m_velocities[vs[2]];
            
            double time, edge_alpha, rel_disp;
            Vec2d normal;
            
            assert( vs[1] < vs[2] );
            
            if ( point_segment_collision( m_positions[vs[0]], m_newpositions[vs[0]], vs[0],
                                         m_positions[vs[1]], m_newpositions[vs[1]], vs[1],
                                         m_positions[vs[2]], m_newpositions[vs[2]], vs[2],
                                         edge_alpha,
                                         normal, time, rel_disp ) )                                 
            {
               
               collision.normal = normal;
               collision.edge_alpha = edge_alpha;
               collision_still_exists = true;                           
            }
            
         } // for collisions


         if ( collision_still_exists )
         {

            std::cout << "Getting desperate!" << std::endl;
            
            // compute one velocity and assign it to all vertices
            
            compute_average_motion( impact_zones[i] );

            
            // NOTE: Collisions may still be registered here if vertices are 
            // in a "colliding" state at the beginning of the timestep. This can
            // happen if mesh improvement moves vertices close to each other.
            
            collision_still_exists = false;
            
            for ( unsigned int c = 0; c < impact_zones[i].collisions.size(); ++c )
            {
               Collision& collision = impact_zones[i].collisions[c];
               const Vec3ui& vs = collision.vertex_indices;
               m_newpositions[vs[0]] = m_positions[vs[0]] + dt * m_velocities[vs[0]];
               m_newpositions[vs[1]] = m_positions[vs[1]] + dt * m_velocities[vs[1]];
               m_newpositions[vs[2]] = m_positions[vs[2]] + dt * m_velocities[vs[2]];
               double time, edge_alpha, rel_disp;
               Vec2d normal;
               
               assert( vs[1] < vs[2] );
               
               if ( point_segment_collision( m_positions[vs[0]], m_newpositions[vs[0]], vs[0],
                                            m_positions[vs[1]], m_newpositions[vs[1]], vs[1],
                                            m_positions[vs[2]], m_newpositions[vs[2]], vs[2],
                                            edge_alpha,
                                            normal, time, rel_disp ) )                                 
               {
                  std::cout << "El Topo: ERROR: " << std::endl;
                  std::cout << "Collision persists even when all vertices have the same velocity" << std::endl;
                  std::cout << "This is likely a problem with rounding error not being treated properly" << std::endl;
                  exit(1);
               }
               
            } // for collisions
            
            
               
         }
         else
         {
            std::cout << "RIZ: success" << std::endl;
         }
         
         
      }  // for IZs
      
      total_collisions.clear();
      
      if ( !finished_detecting_collisions )
      {
         if ( m_verbose ) { std::cout << "attempting to finish global collision detection" << std::endl; }
         
         finished_detecting_collisions = detect_collisions( total_collisions );
         impact_zones.clear();
      }
      else
      {
         detect_new_collisions( impact_zones, total_collisions );
      }

      if ( m_verbose ) { std::cout << total_collisions.size() << " new collisions detected " << std::endl; }
      
   }
   
   return true;
   
}
   
   
void DynamicSurface::rebuild_static_broad_phase( )
{
   m_broad_phase->update_broad_phase_static( *this );
}

void DynamicSurface::rebuild_continuous_broad_phase( )
{
   m_broad_phase->update_broad_phase_continuous( *this );
}

void DynamicSurface::update_static_broad_phase( unsigned int vertex_index )
{
   const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[ vertex_index ];
   
   Vec2d low, high;
   vertex_static_bounds( vertex_index, low, high );
   m_broad_phase->update_vertex( vertex_index, low, high );
        
   for ( unsigned int e = 0; e < incident_edges.size(); ++e )
   {
      edge_static_bounds( incident_edges[e], low, high );
      m_broad_phase->update_edge( incident_edges[e], low, high );
   }
}


void DynamicSurface::update_continuous_broad_phase( unsigned int vertex_index )
{
   const std::vector<unsigned int>& incident_edges = m_mesh.vtxedge[ vertex_index ];
   
   Vec2d low, high;
   vertex_continuous_bounds( vertex_index, low, high );
   m_broad_phase->update_vertex( vertex_index, low, high );
   
   for ( unsigned int e = 0; e < incident_edges.size(); ++e )
   {
      edge_continuous_bounds( incident_edges[e], low, high );
      m_broad_phase->update_edge( incident_edges[e], low, high );
   }      
}


// ---------------------------------------------------------
///
/// Compute the (padded) AABB of a vertex
///
// ---------------------------------------------------------

void DynamicSurface::vertex_static_bounds(unsigned int v, Vec2d &xmin, Vec2d &xmax) const
{
   if ( m_mesh.vtxedge[v].size() == 0 )
   {
      xmin = Vec2d( m_proximity_epsilon );
      xmax = -Vec2d( m_proximity_epsilon );
   }
   else
   {
      xmin = m_positions[v] - Vec2d( m_proximity_epsilon );
      xmax = m_positions[v] + Vec2d( m_proximity_epsilon );
   }
}

// ---------------------------------------------------------
///
/// Compute the AABB of an edge
///
// ---------------------------------------------------------

void DynamicSurface::edge_static_bounds(unsigned int e, Vec2d &xmin, Vec2d &xmax) const
{
   const Vec2ui& edge = m_mesh.edges[e];
   if ( edge[0] == edge[1] )
   {
      xmin = Vec2d( m_proximity_epsilon );
      xmax = -Vec2d( m_proximity_epsilon );
   }
   else
   {            
      minmax(m_positions[edge[0]], m_positions[edge[1]], xmin, xmax);
      xmin -= Vec2d( m_proximity_epsilon );
      xmax += Vec2d( m_proximity_epsilon );
   }
}


// ---------------------------------------------------------
///
/// Compute the AABB of a continuous vertex
///
// ---------------------------------------------------------

void DynamicSurface::vertex_continuous_bounds(unsigned int v, Vec2d &xmin, Vec2d &xmax) const
{
   if ( m_mesh.vtxedge[v].size() == 0 )
   {
      xmin = Vec2d( m_proximity_epsilon );
      xmax = -Vec2d( m_proximity_epsilon );
   }
   else
   {
      minmax(m_positions[v], m_newpositions[v], xmin, xmax);
      xmin -= Vec2d( m_proximity_epsilon );
      xmax += Vec2d( m_proximity_epsilon );
   }
}

// ---------------------------------------------------------
///
/// Compute the AABB of a continuous edge
///
// ---------------------------------------------------------

void DynamicSurface::edge_continuous_bounds(unsigned int e, Vec2d &xmin, Vec2d &xmax) const
{

   const Vec2ui& edge = m_mesh.edges[e];   

   if ( edge[0] == edge[1] )
   {
      xmin = Vec2d( m_proximity_epsilon );
      xmax = -Vec2d( m_proximity_epsilon );
   }
   else
   {      
      minmax(m_positions[edge[0]], m_newpositions[edge[0]], m_positions[edge[1]], m_newpositions[edge[1]], xmin, xmax);
      xmin -= Vec2d( m_proximity_epsilon );
      xmax += Vec2d( m_proximity_epsilon );
   }
}


   
// ---------------------------------------------------------
///
/// Advance mesh by one time step 
///
// ---------------------------------------------------------

void DynamicSurface::integrate( double dt )
{     
      
   assert( m_positions.size() == m_velocities.size() );
   
   assert_mesh_is_intersection_free();
  
   // Handle proximities

   // Set velocities
   for(unsigned int i = 0; i < m_positions.size(); i++)
   {
      m_velocities[i] = ( m_newpositions[i] - m_positions[i] ) / dt;
   }
   
   if ( m_collision_safety )
   {
      handle_edge_point_proximities(dt);
   }
   
   // update new positions with perturbed velocities
   for(unsigned int i = 0; i < m_positions.size(); i++)
   {
      m_newpositions[i] = m_positions[i] + dt * m_velocities[i];
   }

   rebuild_continuous_broad_phase();

   if ( m_collision_safety )
   {        
      // Handle continuous collisions

      bool all_collisions_handled = false;

      handle_point_vs_solid_edge_collisions( dt );
      all_collisions_handled = handle_collisions( dt );

      // failsafe impact zones
      
      if ( !all_collisions_handled )
      {
         all_collisions_handled = handle_collisions_simultaneous(dt);            
      }

      if ( !all_collisions_handled )
      {
         rigid_impact_zones( dt );
      }

      std::vector<Collision> collisions;
      detect_collisions(collisions);
      assert( collisions.size() == 0 );
      
   }

   // Set m_positions
   for(unsigned int i = 0; i < m_positions.size(); i++)
   {
      m_positions[i] = m_newpositions[i];
   } 

   std::cout << "finished integrating" << std::endl;
   
   assert_mesh_is_intersection_free();

      
}


   
// ---------------------------------------------------------
///
/// Fire an assert if any edge is intersecting another edge
///
// ---------------------------------------------------------

void DynamicSurface::assert_mesh_is_intersection_free( )
{
   
   if ( m_verbose ) { std::cout << "checking for intersections... " << std::endl; }
      
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      const Vec2ui& edge_a = m_mesh.edges[i];
      
      if ( edge_a[0] == edge_a[1] )    { continue; }
      
      assert( m_mesh.get_edge( edge_a[0], edge_a[1] ) != m_mesh.edges.size() );

      for ( unsigned int j = 0; j < m_mesh.edges.size(); ++j )
      {
         
         const Vec2ui& edge_b = m_mesh.edges[ j ];
         
         if ( edge_b[0] == edge_b[1] )    { continue; }
         
         if (    edge_b[0] == edge_a[0] || edge_b[0] == edge_a[1] 
              || edge_b[1] == edge_a[0] || edge_b[1] == edge_a[1] )
         {
            continue;
         }
         
         double s01, s23;
         
         if ( segment_segment_intersection( m_positions[edge_a[0]], edge_a[0], 
                                            m_positions[edge_a[1]], edge_a[1],
                                            m_positions[edge_b[0]], edge_b[0], 
                                            m_positions[edge_b[1]], edge_b[1],
                                            s01, s23 ) )
         {   
            
            if ( m_collision_safety )
            {
               Vec2d x0 = m_positions[edge_a[0]];
               Vec2d x1 = m_positions[edge_a[1]];
               Vec2d x2 = m_positions[edge_b[0]];
               Vec2d x3 = m_positions[edge_b[1]];

               double x10=x1[0]-x0[0], y10=x1[1]-x0[1];
               double x31=x3[0]-x1[0], y31=x3[1]-x1[1];
               double x32=x3[0]-x2[0], y32=x3[1]-x2[1];
               double det=x32*y10-x10*y32;
               
               if ( fabs(det) < 1e-8 )
               {
                  // degenerate case
                  continue;
               }
               
               std::cout << "Intersection!  Edge " << edge_a << " vs edge " << edge_b << std::endl;     
               
               std::cout << "x0: " << x0 << std::endl;
               std::cout << "x1: " << x1 << std::endl;
               std::cout << "x2: " << x2 << std::endl;
               std::cout << "x3: " << x3 << std::endl;
               std::cout << "params: " << s01 << ", " << s23 << std::endl;
               
               double unclamped_s01=(x31*y32-x32*y31)/det;
               double unclamped_s23=(x31*y10-x10*y31)/det;

               std::cout << "det: " << det << std::endl;
               std::cout << "unclamped_s01: " << unclamped_s01 << std::endl;
               std::cout << "unclamped_s23: " << unclamped_s23 << std::endl;
               
               assert(0);
            }
            
         }
      }
   }
    
}

// ---------------------------------------------------------
///
/// Using m_newpositions as the geometry, fire an assert if any edge is intersecting any triangles
///
// ---------------------------------------------------------

void DynamicSurface::assert_predicted_mesh_is_intersection_free( )
{
   std::cout << "WARNING: checking intersections at end time disabled" << std::endl;
}

   
} // namespace eltopo2d



