
// ---------------------------------------------------------
//
//  
//  
//
// ---------------------------------------------------------

#ifndef EDGEMESH_H
#define EDGEMESH_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <vec.h>

// ---------------------------------------------------------
//  Interface declarations
// ---------------------------------------------------------

namespace eltopo2d
{
   
// --------------------------------------------------------
///
/// Connectivity information for an edge mesh.  Contains no information on the vertex locations in space.
///
// --------------------------------------------------------

struct EdgeMesh
{  
   EdgeMesh() :
      edges(0), vtxedge(0)
   {}
   
   void clear();
   
   void clear_connectivity();
   void update_connectivity( unsigned int nv );
   
   void check_orientation();
   
   /// Find the index of an edge in the list of edges
   ///
   unsigned int get_edge(unsigned int vtx0, unsigned int vtx1) const;  
   inline unsigned int get_edge( const Vec2ui& e ) const;  
   
   unsigned int nondestructive_add_edge(unsigned int vtx0, unsigned int vtx1);
   void nondestructive_remove_edge( unsigned int edge_index );
   
   unsigned int nondestructive_add_vertex( );
   void nondestructive_remove_vertex(unsigned int vtx);
   
   
   /// Remove primitives which have been deleted (defrag data structures)
   ///
   void nondestructive_clear_unused();
   
   /// Dump some information about this NonDestructiveTriMesh to stdout
   ///
   void print_state();
   
   // ---------------------------------------------------------
   // Data members
   
   /// List of triangles: the fundamental data
   ///
   std::vector<Vec2ui> edges;
      
   /// Edges incident on vertices (given a vertex, which edges is it incident on)
   ///
   std::vector<std::vector<unsigned int> > vtxedge; 
      
};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

inline unsigned int EdgeMesh::get_edge( const Vec2ui& e ) const
{
   return get_edge( e[0], e[1] );
}
   
   
} // namespace

#endif
