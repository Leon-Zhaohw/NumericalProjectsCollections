
// ---------------------------------------------------------
//
//  dynamicsurface.h
//  
//  
//  
//
// ---------------------------------------------------------

#ifndef DYNAMICSURFACE2D_H
#define DYNAMICSURFACE2D_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <deque>
#include <map>

#include "mat.h"
#include "edgemesh.h"


namespace eltopo2d
{
   
// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

class BroadPhaseGrid;
   
// A potentially colliding pair of primitives: (edge, vertex)
typedef std::deque<Vec2ui> CollisionCandidateSet;

// ---------------------------------------------------------
//  Interface declarations
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// A collision
///
// --------------------------------------------------------

struct Collision
{
   
   Collision( const Vec3ui& in_vertex_indices, const Vec2d& in_normal, const double in_edge_alpha ) :
      vertex_indices( in_vertex_indices ),
      normal( in_normal ),
      edge_alpha( in_edge_alpha )
   {
   }

   // One or more vertices is shared between this Collision and other
   inline bool overlap_vertices( const Collision& other ) const;
   
   // ALL vertices are shared between this Collision and other
   inline bool same_vertices( const Collision& other ) const;
      
   // Which vertices are involved in the collision: vertex, edge, edge
   Vec3ui vertex_indices;
   
   // Collision normal
   Vec2d normal;
   
   // Barycentric coordinate of the point of intersection
   double edge_alpha;
   
};

// --------------------------------------------------------
///
/// Used in the simultaneous handling of collisions: a set of connected elements which are in collision
///
// --------------------------------------------------------

struct ImpactZone
{
   ImpactZone() :
      collisions(),
      all_solved( false )
   {}
   
   // Get the set of all vertices in this impact zone
   void get_all_vertices( std::vector<unsigned int>& vertices ) const;
      
   // Whether this ImpactZones shares vertices with other
   bool share_vertices( const ImpactZone& other ) const;

   // Set of collisions with connected vertices
   std::vector<Collision> collisions;  
   
   // Whether all collisions in this zone have been solved (i.e. no longer colliding)
   bool all_solved;
   
};




// --------------------------------------------------------
///
/// A surface mesh.  Essentially consists of a NonDestructiveEdgeMesh object coupled with a set of vertex locations in 3D space.
///
// --------------------------------------------------------

class DynamicSurface
{

public:
 
   DynamicSurface(    unsigned int num_vertices,
                      const double* vertices,
                      unsigned int num_edges,
                      const int* edges,
                      const double* masses,                  
                      double in_proximity_epsilon,
                      bool in_collision_safety,
                      bool in_verbose );
   
   DynamicSurface( const std::vector<Vec2d>& vs, 
                   const std::vector<Vec2ui>& es, 
                   const std::vector<double>& masses,
                   double in_proximity_epsilon,
                   bool in_collision_safety,
                   bool in_verbose );

private:
   
   // Disallowed, do not implement
   DynamicSurface( const DynamicSurface& );
   DynamicSurface& operator=( const DynamicSurface& );
   
public:
   
   virtual ~DynamicSurface() {}

   // ---------------------------------------------------------
   // Simulation step

   /// (Implemented in SurfTrack)
   /// 
   virtual void improve_mesh( ) {}

   /// Advance from current state to collision-free state as close as possible to predicted state.
   /// 
   virtual void integrate(double dt);
   
   // ---------------------------------------------------------
   // Mesh bookkeeping

   unsigned int add_edge(const Vec2ui& t);
   void remove_edge(unsigned int t);  

   unsigned int add_vertex( const Vec2d& new_vertex_position, 
                            const Vec2d& new_vertex_velocity, 
                            double new_vertex_mass, 
                            unsigned int volume_id = 0 );
   
   void remove_vertex(unsigned int v);

   void clear_deleted_vertices( );

   // ---------------------------------------------------------
   // Utility
   
   inline Vec2d get_edge_normal(unsigned int edge) const;
   inline Vec2d get_edge_normal(const Vec2ui& edge) const;
   inline Vec2d get_edge_normal(unsigned int v0, unsigned int v1) const;
   
   inline Vec2d get_vertex_normal( unsigned int vertex ) const;
   
   inline double get_edge_length( unsigned int edge_index ) const;
   inline double get_total_edge_length( ) const;
   inline double get_average_edge_length( ) const;
   inline double get_average_non_solid_edge_length( ) const;
   
   /// Determine volume IDs for all vertices
   void partition_volumes();
   
   inline double get_area() const;
   inline double get_length() const;
   inline double get_non_solid_length( ) const;
             
   void compute_all_edge_normals( std::vector<Vec2d>& midpoints, std::vector<Vec2d>& normals );
   
   double get_curvature( unsigned int vertex_id ) const;

   // ---------------------------------------------------------
   // Query

   double get_distance_to_surface( const Vec2d& x ) const;
   bool point_is_inside_volume( const Vec2d& x ) const;

   // ss are 0-1 alphas, not the geometric distance along the segment
   void get_segment_collisions( const Vec2d& segment_point_a, 
                                const Vec2d& segment_point_b, 
                                std::vector<double>& ss, 
                                std::vector<unsigned int>& edges ) const;
   
   // ---------------------------------------------------------
   // Proximity / collision detection and resolution

   void apply_edge_point_impulse( const Vec2ui& e, unsigned int v, double edge_alpha, Vec2d& direction, double magnitude );

   void handle_edge_point_proximities( double dt );

   void add_point_candidates(unsigned int v, CollisionCandidateSet& collision_candidates);  
   void add_point_update_candidates(unsigned int v, CollisionCandidateSet& collision_candidates);
   void add_edge_candidates(unsigned int e, CollisionCandidateSet& collision_candidates);

   // Impulse-based collision resolution on individual collisions [Bridson 2002]
   // return true if we think we've handled all collisions
   bool handle_collisions(double dt);

   // run one sweep of collision handling, considering only collisions between movable vertices and nonmovable triangles
   void handle_point_vs_solid_edge_collisions( double dt );
   
   // ---------------------------------------------------------
   // Simulataneous collision detection [Harmon 2008]
   
   // detect all collisions
   bool detect_collisions( std::vector<Collision>& collisions );
   
   // detect collisions among vertices present in impact_zones
   void detect_new_collisions( const std::vector<ImpactZone> impact_zones, std::vector<Collision>& collisions );

   // merge impact zones with common vertices
   void merge_impact_zones( std::vector<ImpactZone>& impact_zones, std::vector<ImpactZone>& new_impact_zones );
   
   // iteratively run collision detection and inelastic projection on an active set of collisions
   bool iterated_inelastic_projection( ImpactZone& iz, double dt );

   // attempt to set normal velocity to zero for all collisions in the impact zone
   bool inelastic_projection( const ImpactZone& iz );
  
   // detect and solve all collisions
   bool handle_collisions_simultaneous(double dt);
   
   
   // ---------------------------------------------------------
   // Rigid impact zones
    
   bool compute_rigid_motion( const ImpactZone& iz );
   
   bool compute_average_motion( const ImpactZone& iz );
   
   bool rigid_impact_zones( double dt );
   
   
   // ---------------------------------------------------------
   // Broadphase

   void rebuild_static_broad_phase( ); 
   void rebuild_continuous_broad_phase( );
   
   void update_static_broad_phase( unsigned int vertex_index );   void update_continuous_broad_phase( unsigned int vertex_index );
   
   void vertex_static_bounds(unsigned int v, Vec2d &xmin, Vec2d &xmax) const;
   void edge_static_bounds(unsigned int e, Vec2d &xmin, Vec2d &xmax) const;
   
   void vertex_continuous_bounds(unsigned int v, Vec2d &xmin, Vec2d &xmax) const;
   void edge_continuous_bounds(unsigned int e, Vec2d &xmin, Vec2d &xmax) const;
   
   
   // ---------------------------------------------------------
   // Intersection detection 
   
   void assert_mesh_is_intersection_free( );             // uses m_positions
   void assert_predicted_mesh_is_intersection_free();    // uses m_newpositions

   // ---------------------------------------------------------
   // Data members
   
   // Elements closer than this are proximal
   double m_proximity_epsilon;

   // Dump lots of details to stdout
   bool m_verbose;
   
   // Ensure that no mesh elements intersect, during mesh moving and mesh maintenance
   bool m_collision_safety;
   
   // Vertex positions, predicted locations, velocities and masses
   std::vector<Vec2d> m_positions, m_newpositions, m_velocities;
   std::vector<double> m_masses;
   
   // associate each vertex with one particular volume
   std::vector<unsigned int> m_volume_ids;
   
   // each volume is a vector of vertex indices
   std::vector< std::vector< unsigned int> > m_volumes;
   
   // The mesh "graph"
   EdgeMesh m_mesh;
   
   BroadPhaseGrid *m_broad_phase;
   
};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Determine if another collision has any vertices in common with this collision.
///
// --------------------------------------------------------

inline bool Collision::overlap_vertices( const Collision& other ) const
{
   for ( unsigned short i = 0; i < 3; ++i )
   {
      if ( vertex_indices[i] == other.vertex_indices[0] || 
           vertex_indices[i] == other.vertex_indices[1] || 
           vertex_indices[i] == other.vertex_indices[2] )
      {
         return true;
      }
   }
   
   return false;
}

// --------------------------------------------------------
///
/// Determine if another collision has all the same vertices as this collision.
///
// --------------------------------------------------------

inline bool Collision::same_vertices( const Collision& other ) const
{
   bool found[3];
   for ( unsigned short i = 0; i < 3; ++i )
   {
      if ( vertex_indices[i] == other.vertex_indices[0] || 
           vertex_indices[i] == other.vertex_indices[1] || 
           vertex_indices[i] == other.vertex_indices[2] )
      {
         found[i] = true;
      }
      else
      {
         found[i] = false;
      }
   }
   
   return ( found[0] && found[1] && found[2] );
}

// --------------------------------------------------------
///
/// Extract the set of all vertices in all collisions in an ImpactZone
///
// --------------------------------------------------------

inline void ImpactZone::get_all_vertices( std::vector<unsigned int>& vertices ) const
{
   vertices.clear();
   for ( unsigned int i = 0; i < collisions.size(); ++i )
   {
      add_unique( vertices, collisions[i].vertex_indices[0] );
      add_unique( vertices, collisions[i].vertex_indices[1] );
      add_unique( vertices, collisions[i].vertex_indices[2] );
   }
}


// --------------------------------------------------------
///
/// Determine whether another ImpactZone shares any vertices with this ImpactZone
///
// --------------------------------------------------------

inline bool ImpactZone::share_vertices( const ImpactZone& other ) const
{
   for ( unsigned int i = 0; i < collisions.size(); ++i )
   {
      for ( unsigned int j = 0; j < other.collisions.size(); ++j )
      {
         if ( collisions[i].overlap_vertices( other.collisions[j] ) )
         {
            return true;
         }
      }
   }
   
   return false;
}


// --------------------------------------------------------
///
/// Add a collision to the list
///
// --------------------------------------------------------

inline void add_to_collision_candidates( CollisionCandidateSet& collision_candidates, const Vec2ui& new_collision )
{  
   collision_candidates.push_back( new_collision );
   return;   
}


inline Vec2d DynamicSurface::get_edge_normal(unsigned int edge) const
{
   return get_edge_normal( m_mesh.edges[edge][0], m_mesh.edges[edge][1] ); 
}
   
inline Vec2d DynamicSurface::get_edge_normal(const Vec2ui& edge) const
{
   return get_edge_normal( edge[0], edge[1] );
}

inline Vec2d DynamicSurface::get_edge_normal(unsigned int v0, unsigned int v1) const
{
   Vec2d ev = m_positions[v1] - m_positions[v0];
   ev /= mag(ev);
   return Vec2d( -ev[1], ev[0] );
}
   
// --------------------------------------------------------
///
/// Compute surface normal at the specified vertex (unweighted average of incident triangle normals).
///
// --------------------------------------------------------

inline Vec2d DynamicSurface::get_vertex_normal( unsigned int vertex ) const
{
   Vec2d normal(0);
   for ( unsigned int i = 0; i < m_mesh.vtxedge[vertex].size(); ++i )
   {
      unsigned int v0 = m_mesh.edges[ m_mesh.vtxedge[vertex][i] ][0];
      unsigned int v1 = m_mesh.edges[ m_mesh.vtxedge[vertex][i] ][1];
      Vec2d ev = m_positions[v1] - m_positions[v0];
      double vecmag = mag(ev);
      if(vecmag >= 1e-7)
         ev = ev/vecmag;
      normal += perp( ev );
   }
   
   double len = mag(normal);
   
   if ( len < 1e-7 )
   {
      return Vec2d(0,0);
   }
   else
   {
      return normal / len;
   }
   
}

// --------------------------------------------------------
///
/// Compute length of the specified edge
///
// --------------------------------------------------------

inline double DynamicSurface::get_edge_length( unsigned int edge_index ) const
{
   return mag( m_positions[ m_mesh.edges[edge_index][1] ] - m_positions[ m_mesh.edges[edge_index][0] ] );
}

// --------------------------------------------------------
///
/// Compute average length over all mesh edges
///
// --------------------------------------------------------

inline double DynamicSurface::get_average_edge_length() const
{
   double sum_lengths = 0;
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      const Vec2ui& e = m_mesh.edges[i]; 
      if ( e[0] == e[1] )  { continue; }
      sum_lengths += mag( m_positions[e[1]] - m_positions[e[0]] ); 
   }
   return sum_lengths / (double) m_mesh.edges.size();   
}

// --------------------------------------------------------
///
/// Compute average length over edges on non-solid meshes
///
// --------------------------------------------------------

inline double DynamicSurface::get_average_non_solid_edge_length() const
{
   double sum_lengths = 0;
   unsigned int counted_edges = 0;
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      const Vec2ui& e = m_mesh.edges[i]; 
      if ( e[0] == e[1] )  { continue; }
      if ( m_masses[e[0]] > 100.0 || m_masses[e[1]] > 100.0 ) { continue; }
      sum_lengths += mag( m_positions[e[1]] - m_positions[e[0]] ); 
      ++counted_edges;
   }
   return sum_lengths / (double) counted_edges;   
}

// --------------------------------------------------------
///
/// Compute the area enclosed by this surface
///
// --------------------------------------------------------

inline double DynamicSurface::get_area( ) const
{
   double area = 0.0;
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      if ( m_mesh.edges[i][0] == m_mesh.edges[i][1] ) { continue; }
      Vec2d dx = m_positions[ m_mesh.edges[i][1] ] - m_positions[ m_mesh.edges[i][0] ];
      area += dot( m_positions[ m_mesh.edges[i][0] ], perp(dx) );
   }
   return area;
}

   
// --------------------------------------------------------
   
inline double DynamicSurface::get_length( ) const
{
   double length = 0.0;
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      if ( m_mesh.edges[i][0] == m_mesh.edges[i][1] ) { continue; }
      Vec2d dx = m_positions[ m_mesh.edges[i][1] ] - m_positions[ m_mesh.edges[i][0] ];
      length += mag(dx);
   }
   return length;
}

// --------------------------------------------------------
   
inline double DynamicSurface::get_non_solid_length( ) const
{
   double length = 0.0;
   for ( unsigned int i = 0; i < m_mesh.edges.size(); ++i )
   {
      if ( m_mesh.edges[i][0] == m_mesh.edges[i][1] ) { continue; }
      if ( m_masses[m_mesh.edges[i][0]] > 1.5 || m_masses[m_mesh.edges[i][1]] > 1.5 ) { continue; }
      Vec2d dx = m_positions[ m_mesh.edges[i][1] ] - m_positions[ m_mesh.edges[i][0] ];
      length += mag(dx);
   }
   return length;
}


} // namespace



#endif
