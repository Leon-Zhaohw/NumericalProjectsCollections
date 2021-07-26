
// ---------------------------------------------------------
//
//  
//  
//  
// 
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include "edgemesh.h"

#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Extern globals
// ---------------------------------------------------------


namespace eltopo2d
{

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Clear all mesh information
///
// --------------------------------------------------------

void EdgeMesh::clear()
{
   edges.clear();
   clear_connectivity();
}




// --------------------------------------------------------
///
/// Add a vertex, update connectivity.  Return index of new vertex.
///
// --------------------------------------------------------

unsigned int EdgeMesh::nondestructive_add_vertex( )
{  
   vtxedge.resize( vtxedge.size() + 1 );
   return vtxedge.size() - 1;
}


// --------------------------------------------------------
///
/// Remove a vertex, update connectivity
///
// --------------------------------------------------------

void EdgeMesh::nondestructive_remove_vertex(unsigned int vtx)
{
    vtxedge[vtx].clear();   //edges incident on vertices   
}


// --------------------------------------------------------
///
/// Add an edge to the list.  Return the index of the new edge.
///
// --------------------------------------------------------

unsigned int EdgeMesh::nondestructive_add_edge(unsigned int vtx0, unsigned int vtx1)
{
   int edge_index = edges.size();
   edges.push_back(Vec2ui(vtx0, vtx1));
   
   vtxedge[vtx0].push_back(edge_index);
   vtxedge[vtx1].push_back(edge_index);
   
   return edge_index;
}

   
// --------------------------------------------------------
///
/// Mark an edge as deleted, update connectivity
///
// --------------------------------------------------------

void EdgeMesh::nondestructive_remove_edge( unsigned int edge_index )
{
   // vertex 0
   {
      std::vector<unsigned int>& vertex_to_edge_map = vtxedge[ edges[edge_index][0] ];
      for ( unsigned int i=0; i < vertex_to_edge_map.size(); ++i)
      {
         if ( vertex_to_edge_map[i] == edge_index )
         {
            vertex_to_edge_map.erase( vertex_to_edge_map.begin() + i );
            --i;
         }
      }
   }

   // vertex 1
   {
      std::vector<unsigned int>& vertex_to_edge_map = vtxedge[ edges[edge_index][1] ];
      for ( unsigned int i=0; i < vertex_to_edge_map.size(); ++i)
      {
         if ( vertex_to_edge_map[i] == edge_index )
         {
            vertex_to_edge_map.erase( vertex_to_edge_map.begin() + i );
            --i;
         }
      }
   }

   // TEMP
   
   for ( unsigned int i = 0; i < vtxedge.size(); ++i )
   {
      for ( unsigned int j = 0; j < vtxedge[i].size(); ++j )
      {
         assert( vtxedge[i][j] != edge_index );
      }
   }
        
   edges[edge_index][0] = 0;
   edges[edge_index][1] = 0;  
}


// --------------------------------------------------------
///
/// Find edge specified by two vertices.  Return edges.size if the edge is not found.
///
// --------------------------------------------------------

unsigned int EdgeMesh::get_edge(unsigned int vtx0, unsigned int vtx1) const
{
   assert( vtx0 < vtxedge.size() );
   assert( vtx1 < vtxedge.size() );
   assert( vtx0 != vtx1 );
   
   const std::vector<unsigned int>& incident_edges = vtxedge[vtx0];
   
   for(unsigned int e0 = 0; e0 < incident_edges.size(); e0++)
   {
      unsigned int current_edge_index = incident_edges[e0];
      assert( edges[current_edge_index][0] == vtx0 || edges[current_edge_index][1] == vtx0 );
      
      if ( edges[current_edge_index][0] == vtx1 || edges[current_edge_index][1] == vtx1 )
      {
         return current_edge_index;
      }
      
   }
   
   return edges.size();
}


// --------------------------------------------------------
///
/// Remove primitives which have been deleted by nondestructive_remove_whatever
///
// --------------------------------------------------------

void EdgeMesh::nondestructive_clear_unused()
{
   clear_connectivity();
   
   // This seems to be faster, on average, than clearing and rebuilding the list of edges.
   for( int i = 0; i < (int) edges.size(); i++ )
   {
      if( edges[i][0] == edges[i][1] )
      {
         edges.erase( edges.begin() + i );
         --i;
      }
   }
      
}

// --------------------------------------------------------
///
/// Remove auxiliary connectivity information
///
// --------------------------------------------------------

void EdgeMesh::clear_connectivity()
{
   vtxedge.clear();
}


// --------------------------------------------------------
///
/// Clear and rebuild connectivity information
///
// --------------------------------------------------------

void EdgeMesh::update_connectivity( unsigned int nv )
{
   clear_connectivity();
   
   vtxedge.resize(nv);
   
   for ( unsigned int e = 0; e < edges.size(); ++e )
   {
      const Vec2ui& ce = edges[e];
      if ( ce[0] != ce[1] )
      {
         vtxedge[ce[0]].push_back(e);
         vtxedge[ce[1]].push_back(e);
      }
   }
   
}

   
// --------------------------------------------------------
///
/// 
///
// --------------------------------------------------------
  
void EdgeMesh::check_orientation()
{
   for ( unsigned int i = 0; i < vtxedge.size(); ++i )
   {
      if ( vtxedge[i].size() == 0 ) { continue; }
      
      assert( vtxedge[i].size() == 2 );
      
      const Vec2ui& ea = edges[vtxedge[i][0]];
      const Vec2ui& eb = edges[vtxedge[i][1]];      
      
      if ( ea[0] == i ) 
      {
         assert( eb[0] != i );
         assert( eb[1] == i );
      }
      else
      {
         assert( ea[1] == i );
         assert( eb[0] == i );
         assert( eb[1] != i );
      }
   }
}
   
// --------------------------------------------------------
///
/// Dump some information about this EdgeMesh to stdout
///
// --------------------------------------------------------

void EdgeMesh::print_state()
{   
   int usededges = 0;
   for(int i = 0; i < (int)edges.size(); i++)
   {
      if(edges[i][0] != edges[i][1])
      {
         usededges++;
      }
   }
   
   unsigned int nv = vtxedge.size();
   
   int usedvertices = 0;
   bool *uv = new bool[nv];
   for( unsigned int i = 0; i < nv; i++)
   {
      uv[i] = false;
   }
   
   for( unsigned int i = 0; i < edges.size(); i++)
   {
      if( edges[i][0] != edges[i][1] )
      {
         if(edges[i][0] >= nv || edges[i][1] >= nv)
         {
            usedvertices = -1;
            break;
         };
         
         uv[edges[i][0]] = true;
         uv[edges[i][1]] = true;
      }
   }
   
   if(usedvertices != -1)
   {
      for( unsigned int i = 0; i < nv; i++)
      {
         if(uv[i])
         {
            usedvertices++;
         }
      }
   }
   
   std::cout << "edges: " << usededges << "/" << edges.size() << std::endl;
   std::cout << "vertices: " << usedvertices << "/" << nv << std::endl;
   
   delete[] uv;
   
}


}  // namespace


