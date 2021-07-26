// ---------------------------------------------------------
//
//  
//  
//  
//
// ---------------------------------------------------------

#ifndef ACCELERATIONGRID2D_H
#define ACCELERATIONGRID2D_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include "array2.h"
#include "vec.h"

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Interface declarations
// ---------------------------------------------------------

namespace eltopo2d
{
   
// --------------------------------------------------------
///
/// Regular grid collision culling structure
///
// --------------------------------------------------------

class AccelerationGrid
{

public:
      
   AccelerationGrid();
   ~AccelerationGrid();
   
private:
   AccelerationGrid(AccelerationGrid& other);
   
public:
   
   /// Define the grid given, the extents of the domain and the number of voxels along each dimension
   ///
   void set( const Vec2ui& dims, const Vec2d& xmin, const Vec2d& xmax );
   
   /// Generate a set of voxel indices from a pair of AABB extents
   ///
   void boundstoindices( const Vec2d& xmin, const Vec2d& xmax, Vec2i& xmini, Vec2i& xmaxi);
   
   /// Add an object with the specified index and AABB to the grid
   ///
   void add_element(unsigned int idx, const Vec2d& xmin, const Vec2d& xmax);
   
   /// Remove an object with the specified index from the grid
   ///
   void remove_element(unsigned int idx);
   
   /// Reset the specified object's AABB
   ///
   void update_element(unsigned int idx, const Vec2d& xmin, const Vec2d& xmax);
      
   /// Remove all elements from the grid
   ///
   void clear();
      
   /// Return the set of elements which have AABBs overlapping the query AABB.
   ///
   void find_overlapping_elements( const Vec2d& xmin, const Vec2d& xmax, std::vector<unsigned int>& results );
   
   // TEMP statistics-gathering functions
   
   float num_elements_per_cell();
   float num_cells_per_element();
   unsigned int max_elements_per_cell();
   
   
   /// Each cell contains an array of indices specifying the elements whose AABBs overlap the cell
   ///
   Array2< std::vector<unsigned int>* > cells;
   
   /// For each element, a list of index pairs, each pair specifying a cell which overlaps the element. 
   ///
   std::vector<std::vector<Vec2ui> > elementidxs;
   
   /// Element AABBs
   ///
   std::vector<Vec2d> elementxmins, elementxmaxs;
   
   /// For each element, the timestamp of the last query that examined the element
   ///
   std::vector<unsigned int> elementquery;
   
   /// Timestamp of the last query
   ///
   unsigned int lastquery;
   
   /// Lower/upper corners of the entire grid
   ///
   Vec2d gridxmin, gridxmax;
   
   /// Cell dimensions
   ///
   Vec2d cellsize;
   
   /// Inverse cell dimensions
   ///
   Vec2d invcellsize;

};

}  // namespace

inline bool aabb_overlap( const Vec2d& xmin, const Vec2d& xmax, const Vec2d& oxmin, const Vec2d& oxmax )
{
   if( (xmin[0] <= oxmax[0] && xmin[1] <= oxmax[1]) &&
       (xmax[0] >= oxmin[0] && xmax[1] >= oxmin[1]) )
   {
      return true;
   }
   
   return false;
   
}
   


#endif
