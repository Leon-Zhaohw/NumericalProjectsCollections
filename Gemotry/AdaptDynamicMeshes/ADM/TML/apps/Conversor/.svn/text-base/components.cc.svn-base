#include <tml.h>
#include <iostream>
#include <fstream>
#include <stack>

// enumerate components from 0 to (num_components-1)
int compute_num_components(TML::TriMesh<>& mesh)
{
  int num = -1;
  TML::Vertex* w = NULL;
  TML::Halfedge* h;
  TML::TriMesh<>::VertexIter vi;
  std::stack<TML::Vertex*> Q;
  for (vi=mesh.verts_begin(); vi!=mesh.verts_end(); vi++)
    {
      if ( (*vi)->mark()!=-1 ) continue;
      num++;
      Q.push( *vi );
      while (!Q.empty())
        {
          w = Q.top(); Q.pop();
          w->set_mark(num);
          for (h=w->star_first(); h!=NULL; h=w->star_next(h))
            {
              if (h->org()->mark()!=-1) continue;
              Q.push( h->org() );
            }
        }
    }
  return (num+1);
}

//----------------------------------------------------------//

int main(int argc, char** argv)
{
  if (argc!=3)
    throw TML::Error("USAGE: ./components <input> <output>");
  
  TML::TriMesh<> mesh(argv[1]);
  
  int num_comp = compute_num_components(mesh);
  
  int i, imax, max;
  TML::TriMesh<>::VertexIter vi;

  // computing size of each component
  int len[num_comp];
  for (i=0; i<num_comp; i++) len[i] = 0;
  for (vi=mesh.verts_begin(); vi!=mesh.verts_end(); vi++) len[ (*vi)->mark() ]++;
  
  // finding the biggest component
  imax = 0; max = len[0];
  for (i=1; i<num_comp; i++) 
    if (len[i]>max) { imax = i; max = len[i]; }
  
  ////////////////////////////////////
  // removing the biggest component //
  ////////////////////////////////////
  TML::TriMeshIO new_mesh;
  
  // getting vertices
  int counter = 0;
  TML::R3 p;
  for (vi=mesh.verts_begin(); vi!=mesh.verts_end(); vi++)
    {
      if ( (*vi)->mark()!=imax ) continue;
      (*vi)->set_id( counter );  counter++;
      p = (*vi)->pos(); new_mesh.add(p);
    }
  
  // getting facets
  TML::TriMesh<>::FacetIter fi;
  for (fi=mesh.facets_begin(); fi!=mesh.facets_end(); fi++)
    {
      if ( (*fi)->vertex(0)->mark()!=(*fi)->vertex(1)->mark() ||
	   (*fi)->vertex(1)->mark()!=(*fi)->vertex(2)->mark() || 
	   (*fi)->vertex(2)->mark()!=(*fi)->vertex(0)->mark()  )
	throw TML::Error("Something wrong ... components are connected :(");
      
      if ( (*fi)->vertex(0)->mark()!=imax ) continue;

      TML::Polygon t(3);
      for (i=0; i<3; i++) t.vertex( i, (*fi)->vertex(i)->id() );
      new_mesh.add( t );
    }
  
  return new_mesh.Write(argv[2]);
}
