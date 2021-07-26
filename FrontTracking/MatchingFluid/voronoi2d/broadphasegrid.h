// ---------------------------------------------------------
//
//  broadphasegrid.h
//  
// ---------------------------------------------------------

#ifndef BROADPHASEGRID2D_H
#define BROADPHASEGRID2D_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include "accelerationgrid.h"

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

namespace eltopo2d
{

   
class DynamicSurface;

// ---------------------------------------------------------
//  Interface declarations
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Broad phase collision detector using three regular grids: one grid each for vertices, edges and triangles.
///
// --------------------------------------------------------

class BroadPhaseGrid 
{
public:
   
   BroadPhaseGrid() :
      vertex_grid(), edge_grid()
   {}

   /// Rebuild the broad phase using current vertex positions
   ///
   void update_broad_phase_static( const DynamicSurface& surface );
   
   /// Rebuild the broad phase using current and predicted vertex positions
   ///
   void update_broad_phase_continuous( const DynamicSurface& surface );

   inline void add_vertex( unsigned int index, const Vec2d& aabb_low, const Vec2d& aabb_high ); 
   inline void add_edge( unsigned int index, const Vec2d& aabb_low, const Vec2d& aabb_high ); 
   
   inline void update_vertex( unsigned int index, const Vec2d& aabb_low, const Vec2d& aabb_high ); 
   inline void update_edge( unsigned int index, const Vec2d& aabb_low, const Vec2d& aabb_high ); 
   
   inline void remove_vertex( unsigned int index ); 
   inline void remove_edge( unsigned int index ); 
   
   /// Get the set of vertices whose bounding volumes overlap the specified bounding volume
   ///
   inline void get_potential_vertex_collisions( const Vec2d& aabb_low, 
                                                const Vec2d& aabb_high, 
                                                std::vector<unsigned int>& overlapping_vertices );
   
   /// Get the set of edges whose bounding volumes overlap the specified bounding volume
   ///
   inline void get_potential_edge_collisions( const Vec2d& aabb_low, 
                                              const Vec2d& aabb_high, 
                                              std::vector<unsigned int>& overlapping_edges );
   
   /// Rebuild one of the grids
   ///
   void build_acceleration_grid( AccelerationGrid& grid, 
                                 std::vector<Vec2d>& xmins, 
                                 std::vector<Vec2d>& xmaxs, 
                                 double length_scale, 
                                 double grid_padding );
   
   /// Regular grids
   ///
   AccelerationGrid vertex_grid;
   AccelerationGrid edge_grid;
   
};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Add a vertex to the broad phase
///
// --------------------------------------------------------

inline void BroadPhaseGrid::add_vertex( unsigned int index, const Vec2d& aabb_low, const Vec2d& aabb_high )
{
   vertex_grid.add_element( index, aabb_low, aabb_high );
}

// --------------------------------------------------------
///
/// Add an edge to the broad phase
///
// --------------------------------------------------------

inline void BroadPhaseGrid::add_edge( unsigned int index, const Vec2d& aabb_low, const Vec2d& aabb_high )
{
   edge_grid.add_element( index, aabb_low, aabb_high );
}


inline void BroadPhaseGrid::update_vertex( unsigned int index, const Vec2d& aabb_low, const Vec2d& aabb_high )
{
   vertex_grid.update_element( index, aabb_low, aabb_high );
}

inline void BroadPhaseGrid::update_edge( unsigned int index, const Vec2d& aabb_low, const Vec2d& aabb_high )
{
   edge_grid.update_element( index, aabb_low, aabb_high );   
}


// --------------------------------------------------------
///
/// Remove a vertex from the broad phase
///
// --------------------------------------------------------

inline void BroadPhaseGrid::remove_vertex( unsigned int index )
{
   vertex_grid.remove_element( index );
}

// --------------------------------------------------------
///
/// Remove an edge from the broad phase
///
// --------------------------------------------------------

inline void BroadPhaseGrid::remove_edge( unsigned int index )
{
   edge_grid.remove_element( index );
}


// --------------------------------------------------------
///
/// Query the broad phase to get the set of all vertices overlapping the given AABB
///
// --------------------------------------------------------

inline void BroadPhaseGrid::get_potential_vertex_collisions( const Vec2d& aabb_low, const Vec2d& aabb_high, std::vector<unsigned int>& overlapping_vertices )
{
   vertex_grid.find_overlapping_elements( aabb_low, aabb_high, overlapping_vertices );
}

// --------------------------------------------------------
///
/// Query the broad phase to get the set of all edges overlapping the given AABB
///
// --------------------------------------------------------

inline void BroadPhaseGrid::get_potential_edge_collisions( const Vec2d& aabb_low, const Vec2d& aabb_high, std::vector<unsigned int>& overlapping_edges )
{
   edge_grid.find_overlapping_elements( aabb_low, aabb_high, overlapping_edges );
}



}  // namespace

#endif




