// ---------------------------------------------------------
//
//  accelerationgrid.cpp
//  
//  A grid-based collision test culling structure.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include "accelerationgrid.h"

#include <vector>

#include <array2.h>
#include <util.h>
#include <vec.h>
#include <wallclocktime.h>
#include <limits>

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
/// Default constructor
///
// --------------------------------------------------------

AccelerationGrid::AccelerationGrid() :
   cells(0,0),
   elementidxs(0),
   elementxmins(0),
   elementxmaxs(0),
   elementquery(0),
   lastquery(0),
   gridxmin(0,0),
   gridxmax(0,0),
   cellsize(0,0),
   invcellsize(0,0)
{
   Vec2ui dims(1,1);
   Vec2d xmin(0,0), xmax(1,1);
   set(dims, xmin, xmax);
}


// --------------------------------------------------------
///
/// Destructor: clear all grids
///
// --------------------------------------------------------

AccelerationGrid::~AccelerationGrid()
{
   clear();
}

// --------------------------------------------------------
///
/// Define the grid, given the extents of the domain and the number of desired voxels along each dimension
///
// --------------------------------------------------------

void AccelerationGrid::set( const Vec2ui& dims, const Vec2d& xmin, const Vec2d& xmax )
{
   gridxmin = xmin;
   gridxmax = xmax;
   
   for(unsigned int i = 0; i < 2; i++)
   {
      cellsize[i] = (gridxmax[i]-gridxmin[i])/dims[i];
      invcellsize[i] = 1.0 / cellsize[i];
   }
   
   clear();
   
   cells.resize(dims[0], dims[1]);    
   for(unsigned int i = 0; i < cells.a.size(); i++)
   {
      cells.a[i] = 0;
   }   
}

// --------------------------------------------------------
///
/// Generate a set of voxel indices from a pair of AABB extents
///
// --------------------------------------------------------

void AccelerationGrid::boundstoindices(const Vec2d& xmin, const Vec2d& xmax, Vec2i& xmini, Vec2i& xmaxi)
{
   
   xmini[0] = (int) floor((xmin[0] - gridxmin[0]) * invcellsize[0]);
   xmini[1] = (int) floor((xmin[1] - gridxmin[1]) * invcellsize[1]);
   
   xmaxi[0] = (int) floor((xmax[0] - gridxmin[0]) * invcellsize[0]);
   xmaxi[1] = (int) floor((xmax[1] - gridxmin[1]) * invcellsize[1]);
      
   if(xmini[0] < 0) xmini[0] = 0;
   if(xmini[1] < 0) xmini[1] = 0;
   
   if(xmaxi[0] >= cells.ni) xmaxi[0] = cells.ni-1;
   if(xmaxi[1] >= cells.nj) xmaxi[1] = cells.nj-1;
}

// --------------------------------------------------------
///
/// Add an object with the specified index and AABB to the grid
///
// --------------------------------------------------------

void AccelerationGrid::add_element(unsigned int idx, const Vec2d& xmin, const Vec2d& xmax)
{
   
   if(elementidxs.size() <= idx)
   {
      elementidxs.resize(idx+1);
      elementxmins.resize(idx+1);
      elementxmaxs.resize(idx+1);
      elementquery.resize(idx+1);
   }
   
   elementxmins[idx] = xmin;
   elementxmaxs[idx] = xmax;
   elementquery[idx] = 0;
   
   Vec2i xmini, xmaxi;
   boundstoindices(xmin, xmax, xmini, xmaxi);
   
   elementidxs[idx].clear();

   for(int i = xmini[0]; i <= xmaxi[0]; i++)
   {
      for(int j = xmini[1]; j <= xmaxi[1]; j++)
      {
         std::vector<unsigned int>*& cell = cells(i, j);
         if(!cell)
            cell = new std::vector<unsigned int>();
         
         cell->push_back(idx);
         
         elementidxs[idx].push_back(Vec2ui(i, j));
      }
   }
}

// --------------------------------------------------------
///
/// Remove an object with the specified index from the grid
///
// --------------------------------------------------------

void AccelerationGrid::remove_element(unsigned int idx)
{
   for(unsigned int c = 0; c < elementidxs[idx].size(); c++)
   {
      Vec2ui cellcoords = elementidxs[idx][c];
      std::vector<unsigned int>* cell = cells(cellcoords[0], cellcoords[1]);
      
      std::vector<unsigned int>::iterator it = cell->begin();
      while(*it != idx)
      {
         it++;
      }
      
      cell->erase(it);
   }
   
   elementidxs[idx].clear();
}

// --------------------------------------------------------
///
/// Reset the specified object's AABB
///
// --------------------------------------------------------

void AccelerationGrid::update_element(unsigned int idx, const Vec2d& xmin, const Vec2d& xmax)
{
   remove_element(idx);
   add_element(idx, xmin, xmax);
}

// --------------------------------------------------------
///
/// Remove all elements from the grid
///
// --------------------------------------------------------

void AccelerationGrid::clear()
{

   for(unsigned int i = 0; i < cells.a.size(); i++)
   {
      std::vector<unsigned int>*& cell = cells.a[i];  
      if(cell)
      {
         delete cell;
         cell = 0;
      }
   }
   
   elementidxs.clear();
   elementxmins.clear();
   elementxmaxs.clear();
   elementquery.clear();
   lastquery = 0;
   
}

// --------------------------------------------------------
///
/// Return the set of elements which have AABBs overlapping the query AABB.
///
// --------------------------------------------------------

void AccelerationGrid::find_overlapping_elements( const Vec2d& xmin, const Vec2d& xmax, std::vector<unsigned int>& results ) 
{
   if(lastquery == std::numeric_limits<unsigned int>::max())
   {
      std::vector<unsigned int>::iterator iter = elementquery.begin();
      for( ; iter != elementquery.end(); ++iter )
      {
         *iter = 0;
      }
      lastquery = 0;
   }

   ++lastquery;
   
   results.clear();
      
   Vec2i xmini, xmaxi;
   boundstoindices(xmin, xmax, xmini, xmaxi);
   
   for(int i = xmini[0]; i <= xmaxi[0]; ++i)
   {
      for(int j = xmini[1]; j <= xmaxi[1]; ++j)
      {
         std::vector<unsigned int>* cell = cells(i, j);
         
         if(cell)
         {
            for( std::vector<unsigned int>::const_iterator citer = cell->begin(); citer != cell->end(); ++citer)
            {
               unsigned int oidx = *citer;
               
               // Check if the object has already been found during this query
               
               if(elementquery[oidx] < lastquery)
               {
                  
                  // Object has not been found.  Set elementquery so that it will not be tested again during this query.
                  
                  elementquery[oidx] = lastquery;
               
                  const Vec2d& oxmin = elementxmins[oidx];
                  const Vec2d& oxmax = elementxmaxs[oidx];
                  
                  if( (xmin[0] <= oxmax[0] && xmin[1] <= oxmax[1]) &&
                      (xmax[0] >= oxmin[0] && xmax[1] >= oxmin[1]) )
                  {
                     results.push_back(oidx);
                  }
               
               }
            }
            
         }
      }
   }
}


// --------------------------------------------------------
///
/// TEMP statistics-gathering functions
///
// --------------------------------------------------------

float AccelerationGrid::num_elements_per_cell()
{
   unsigned int num_non_empty_cells = 0;
   unsigned int counted_num_elements = 0;
   
   for ( unsigned int i = 0; i < cells.a.size(); ++i )
   {
      if ( cells.a[i] != NULL && cells.a[i]->size() != 0 )
      {
         ++num_non_empty_cells;
         counted_num_elements += cells.a[i]->size();
      }
   }
   
   return (float) counted_num_elements / (float) num_non_empty_cells;
}


float AccelerationGrid::num_cells_per_element()
{
   unsigned int num_elements = elementidxs.size();
   unsigned int num_cells_referenced = 0;
   
   for ( unsigned int i = 0; i < elementidxs.size(); ++i )
   {
      num_cells_referenced += elementidxs[i].size();
   }
   
   return (float) num_cells_referenced / (float) num_elements;
}


unsigned int AccelerationGrid::max_elements_per_cell()
{
   unsigned int max_elements = 0;
   
   for ( unsigned int i = 0; i < cells.a.size(); ++i )
   {
      if ( cells.a[i] != NULL )
      {
         max_elements = max( max_elements, (unsigned int) cells.a[i]->size() );
      }
   }
   
   return max_elements;
}


}  // namespace

