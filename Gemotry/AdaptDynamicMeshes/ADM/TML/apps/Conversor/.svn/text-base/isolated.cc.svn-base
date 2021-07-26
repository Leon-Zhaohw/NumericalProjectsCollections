#include <tml.h>
#include <iostream>

int main(int argc, char** argv)
{
  if (argc!=3)
    throw TML::Error("USAGE: ./isolated <input> <output>");

  TML::TriMesh<> mesh(argv[1]);
  TML::TriMeshIO new_mesh;

  int i, num_verts = 0;
  TML::TriMesh<>::FacetIter fi;
  for (fi=mesh.facets_begin(); fi!=mesh.facets_end(); fi++)
    {
      TML::R3 p;
      TML::Polygon t(3);
      for (i=0; i<3; i++)
	{
	  if ( (*fi)->vertex(i)->mark()!=-1 ) 
	    {
	      t.vertex( i, (*fi)->vertex(i)->mark() );
	      continue;
	    }

	  (*fi)->vertex(i)->set_mark(num_verts);
	  num_verts++;

	  t.vertex( i, (*fi)->vertex(i)->mark() );
	   
	  p = (*fi)->vertex(i)->pos();
	  new_mesh.add(p);
	}
      new_mesh.add(t);
    }

  return new_mesh.Write( argv[2] );
}
