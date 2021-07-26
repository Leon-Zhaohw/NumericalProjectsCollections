#include <tml.h>
#include <vector>
#include <iostream>
#include <cassert>

int main(int argc, char** argv)
{
  if (argc!=3)
    throw TML::Error("USAGE: ./boundary <input> <output>");

  TML::TriMesh<> mesh(argv[1]);
  
  TML::TriMeshIO new_mesh;
  assert( new_mesh.Read(argv[1]) );
  
  int i, n, counter;
  TML::Vertex* w;
  TML::R3 center;
  TML::TriMesh<>::VertexIter vi;
  std::vector<TML::Vertex*> vrts;
  std::vector<TML::Vertex*>::iterator it0, it1;

  n = mesh.num_verts();
  for (i=0, vi=mesh.verts_begin(); vi!=mesh.verts_end(); vi++, i++)
    {
      if ( !(*vi)->is_bdry() ) continue;
      
      if ( (*vi)->mark()!=-1 ) continue;
      
      //std::cout << "COMECANDO " << i << std::endl;

      // getting boundary cicle
      counter = 0;
      center = TML::R3();      
      w = *vi;
      do{
	counter++;
	center += w->pos();

	assert( w->mark()!=0 );
	w->set_mark(0);
	vrts.push_back( w );
	w = w->star_first()->org();
      }while (w!=(*vi));
      center /= double(counter);

      //std::cout << "CICLO: " << counter << std::endl;

      // creating new triangles
      new_mesh.add(center);
      it0 = vrts.begin(); 
      it1 = it0; it1++;
      for ( ; it1!=vrts.end(); it0++, it1++)
	{
	  TML::Polygon t(3);
	  t.vertex( 0, n );
	  t.vertex( 1, (*it0)->id() );
	  t.vertex( 2, (*it1)->id() );
	  new_mesh.add( t );
	}
      it1 = vrts.begin();
      TML::Polygon t(3);
      t.vertex( 0, n );
      t.vertex( 1, (*it0)->id() );
      t.vertex( 2, (*it1)->id() );
      new_mesh.add( t );
      vrts.clear();
      n++;

      //std::cout << "FECHADO!" << std::endl;
    }

  std::cout << "I add " << n - mesh.num_verts() << " vertices and "
	    << new_mesh.get_facets().size() - mesh.num_facets() << " triangles." << std::endl;
  return new_mesh.Write(argv[2]);
}
