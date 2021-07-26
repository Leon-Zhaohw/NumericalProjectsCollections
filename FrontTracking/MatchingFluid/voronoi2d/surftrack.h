// ---------------------------------------------------------
//
//  
//
//  
//
// ---------------------------------------------------------

#ifndef SURFTRACK2D_H
#define SURFTRACK2D_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include "dynamicsurface.h"

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Interface declarations
// ---------------------------------------------------------

namespace eltopo2d
{
   
   
// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

   
struct SurfTrackInitializationParameters
{
   SurfTrackInitializationParameters() :
      m_proximity_epsilon( 1e-4 ),
      m_verbose( false ),
      m_collision_safety( true ),
      m_max_area_change( 0.001 ),
      m_min_edge_length( 0.05 ),
      m_max_edge_length( 0.15 ),
      m_merge_proximity_epsilon( 1e-3 ),
      m_allow_topology_changes( true ),
      m_perform_improvement( true )
   {}
   
   double m_proximity_epsilon;
   bool m_verbose;
   bool m_collision_safety;
   double m_max_area_change;
   double m_min_edge_length;
   double m_max_edge_length;
   double m_merge_proximity_epsilon;
   bool m_allow_topology_changes;
   bool m_perform_improvement;
   
};
  
    

// ---------------------------------------------------------
///
/// Used to build a list of edges sorted in order of increasing length.
/// 
// ---------------------------------------------------------
   
struct SortableEdge
{
   unsigned int edge_index;
   double edge_length;
   
   SortableEdge( unsigned int ei, double el ) : edge_index(ei), edge_length(el) {}
   
   bool operator<( const SortableEdge& other ) const
   {
      return (this->edge_length < other.edge_length);
   }
};
  
// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------
   
struct SortablePointEdgeProximity
{
   SortablePointEdgeProximity( unsigned int p, unsigned int e, double d ) :
   point_index( p ),
   edge_index( e ),
   distance( d )
   {}
   
   unsigned int point_index;
   unsigned int edge_index;
   double distance;
   
   bool operator<( const SortablePointEdgeProximity& other ) const
   {
      return distance < other.distance;
   }
};

   
// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

   
class SurfTrack : public DynamicSurface
{
   
public:
      

   SurfTrack( unsigned int num_vertices,
              const double* vertices,
              unsigned int num_edges,
              const int* edges,
              const double* masses,
              const SurfTrackInitializationParameters& initial_parameters );
   
   /// Create a SurfTrack object from a set of vertices and edges using the specified paramaters
   ///
   SurfTrack( const std::vector<Vec2d>& vs, 
              const std::vector<Vec2ui>& es, 
              const std::vector<double>& masses,
              const SurfTrackInitializationParameters& initial_parameters );
   
   
private:
   
   SurfTrack( );
   
   // Disallow copying and assignment by declaring private
   //
   SurfTrack( const SurfTrack& );
   SurfTrack& operator=( const SurfTrack& );
   
   
public:
   
   /// Display the surface in OpenGL using the specified options
   ///
   void render( unsigned int options );
      
   /// advance one time step (just calls the version in DynamicSurface)
   ///
   inline void integrate(double dt);
   
   /// run mesh maintenance operations
   ///
   void improve_mesh( );
   
   /// run merging
   ///
   void topology_changes( );
   
   // ---------------------------------------------------------
   // mesh maintenance operations
   // ---------------------------------------------------------
   
   /// Split an edge, using subdivision_scheme to determine the new vertex location, if safe to do so.
   ///
   bool split_edge( unsigned int edge );
   
   /// Split all long edges
   ///
   bool split_pass();
   
   /// Delete an edge by moving its source vertex to its destination vertex
   ///
   bool collapse_edge(unsigned int edge );
   
   /// Collapse all short edges
   ///
   bool collapse_pass();
   
   // ---------------------------------------------------------
   // topological merging
   // ---------------------------------------------------------
   
   /// Attempt to merge
   ///
   bool merge( unsigned int vertex_index, unsigned int edge_index );   
   
   /// Zipper
   ///
   void merge_pass();
   
   // ---------------------------------------------------------
   // mesh cleanup
   // ---------------------------------------------------------
   
   /// 
   ///
   void trim_non_manifold( const std::vector<unsigned int>& edge_indices );
   
   /// 
   ///
   void trim_non_manifold();
   
   /// Find vertices with disconnected neighbourhoods, and pull them apart
   ///
   void separate_singular_vertices();
   
   // ---------------------------------------------------------
   // mesh maintenance helpers
   // ---------------------------------------------------------
   
   // split
   bool split_edge_pseudo_motion_introduces_collision( const Vec2d& new_vertex_position, 
                                                       const Vec2d& new_vertex_smooth_position, 
                                                       unsigned int edge, 
                                                       unsigned int vertex_a, 
                                                       unsigned int vertex_b );

   
   // collapse
   bool collapse_edge_introduces_collision( unsigned int vertex_to_delete, unsigned int vertex_to_keep, unsigned int edge, const Vec2d& vertex_new_position );
   bool collapse_edge_introduces_normal_inversion( unsigned int vertex_to_delete, unsigned int vertex_to_keep, unsigned int edge, const Vec2d& vertex_new_position ) const;
   bool collapse_edge_introduces_area_change( unsigned int vertex_to_delete, unsigned int vertex_to_keep, unsigned int edge, const Vec2d& vertex_new_position );   
   
   // merge
   bool get_zipper_edges( unsigned int vertex_index, unsigned int edge_index, std::vector<Vec2ui>& new_edges );
   bool zippering_introduces_collision( const std::vector<Vec2ui>& new_edges, const std::vector<unsigned int>& deleted_edges );
   bool zippering_introduces_area_change( unsigned int vertex_index, unsigned int edge_index );         
   
   // ---------------------------------------------------------
   // Member variables
   // ---------------------------------------------------------
     
   /// Maximum volume change allowed when flipping or collapsing an edge
   double m_max_area_change;
   
   /// Mimimum edge length.  Edges shorter than this will be collapsed.
   double m_min_edge_length;   
   
   /// Maximum edge length.  Edges longer than this will be subdivided.
   double m_max_edge_length;   
   
   /// Elements within this distance will trigger a merge attempt
   double m_merge_proximity_epsilon;
   
   /// Edges which are involved in connectivity changes which may introduce degeneracies
   std::vector<unsigned int> m_dirty_edges;
   
   /// Whether to allow merging and separation
   bool m_allow_topology_changes;
   
   /// Whether to perform adaptivity operations
   bool m_perform_improvement;
   
};


// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

  
inline void SurfTrack::integrate(double dt)
{
   DynamicSurface::integrate(dt);
}

} // namespace

#endif

