// ---------------------------------------------------------
//
//  broadphasegrid.cpp
//  
//  Broad phase collision detection culling using three regular, volumetric grids.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include "broadphasegrid.h"
#include "dynamicsurface.h"

// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

namespace eltopo2d
{


// --------------------------------------------------------
///
/// Construct one grid from the given set of AABBs, using the given length scale as the cell size, with the given padding
///
// --------------------------------------------------------

void BroadPhaseGrid::build_acceleration_grid( AccelerationGrid& grid, 
                                              std::vector<Vec2d>& xmins, 
                                              std::vector<Vec2d>& xmaxs, 
                                              double length_scale, 
                                              double grid_padding )
{

   Vec2d xmax = xmaxs[0];
   Vec2d xmin = xmins[0];
   double maxdistance = 0;
   
   unsigned int n = xmins.size();
   
   for(unsigned int i = 0; i < n; i++)
   {
      update_minmax(xmins[i], xmin, xmax);
      update_minmax(xmaxs[i], xmin, xmax);
      maxdistance = std::max(maxdistance, mag(xmaxs[i] - xmins[i]));
   }
   
   for(unsigned int i = 0; i < 2; i++)
   {
      xmin[i] -= 2*maxdistance + grid_padding;
      xmax[i] += 2*maxdistance + grid_padding;
   }
   
   Vec2ui dims(1,1);
          
   if(mag(xmax-xmin) > grid_padding)
   {
      for(unsigned int i = 0; i < 2; i++)
      {
         unsigned int d = (unsigned int)ceil((xmax[i] - xmin[i])/length_scale);
         
         if(d < 1) d = 1;
         if(d > n) d = n;
         dims[i] = d;
      }
   }
      
   grid.set(dims, xmin, xmax);
   
   for(unsigned int i = 0; i < xmins.size(); ++i )
   {
      // don't add inside-out AABBs
      if ( xmins[i][0] > xmaxs[i][0] )  { continue; }
      grid.add_element(i, xmins[i], xmaxs[i]);
   }
}


// --------------------------------------------------------
///
/// Rebuild acceleration grids according to the given triangle mesh
///
// --------------------------------------------------------

void BroadPhaseGrid::update_broad_phase_static( const DynamicSurface& surface )
{
   double grid_scale = surface.get_average_edge_length();
   
   {
      unsigned int num_vertices = surface.m_positions.size();
      std::vector<Vec2d> vertex_xmins(num_vertices), vertex_xmaxs(num_vertices);
      for(unsigned int i = 0; i < num_vertices; i++)
      {
         surface.vertex_static_bounds(i, vertex_xmins[i], vertex_xmaxs[i]);
      }
      build_acceleration_grid( vertex_grid, vertex_xmins, vertex_xmaxs, grid_scale, surface.m_proximity_epsilon );
   }
   
   {
      unsigned int num_edges = surface.m_mesh.edges.size();
      std::vector<Vec2d> edge_xmins(num_edges), edge_xmaxs(num_edges);
      for(unsigned int i = 0; i < num_edges; i++)
      {
         surface.edge_static_bounds(i, edge_xmins[i], edge_xmaxs[i]);
      }      
      build_acceleration_grid( edge_grid, edge_xmins, edge_xmaxs, grid_scale, surface.m_proximity_epsilon );
   }
      
}



// --------------------------------------------------------
///
/// Rebuild acceleration grids according to the given triangle mesh
///
// --------------------------------------------------------

void BroadPhaseGrid::update_broad_phase_continuous( const DynamicSurface& surface )
{
   double grid_scale = surface.get_average_edge_length();
   
   {
      unsigned int num_vertices = surface.m_positions.size();
      std::vector<Vec2d> vertex_xmins(num_vertices), vertex_xmaxs(num_vertices);
      for(unsigned int i = 0; i < num_vertices; i++)
      {           
         surface.vertex_continuous_bounds(i, vertex_xmins[i], vertex_xmaxs[i]);
      }
      build_acceleration_grid( vertex_grid, vertex_xmins, vertex_xmaxs, grid_scale, surface.m_proximity_epsilon );
   }
   
   {
      unsigned int num_edges = surface.m_mesh.edges.size();
      std::vector<Vec2d> edge_xmins(num_edges), edge_xmaxs(num_edges);
      for(unsigned int i = 0; i < num_edges; i++)
      {
         surface.edge_continuous_bounds(i, edge_xmins[i], edge_xmaxs[i]);
      }
      build_acceleration_grid( edge_grid, edge_xmins, edge_xmaxs, grid_scale, surface.m_proximity_epsilon );
   }
      
}



}  // namespace

