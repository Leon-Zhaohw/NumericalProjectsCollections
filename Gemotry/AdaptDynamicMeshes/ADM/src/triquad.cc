#include "adm.h"

using namespace ADM;

void mark_link(A48Edge* e)
{
  Halfedge* h0 = e->hedge(0);
  Halfedge* h1 = e->hedge(1);
  if (h0->facet() != NULL) 
    {
      h0->next()->edge()->set_mark(1);
      h0->prev()->edge()->set_mark(1);
    }
  if (h1->facet() != NULL) 
    {
      h1->next()->edge()->set_mark(1);
      h1->prev()->edge()->set_mark(1);
    }
}

void AdaptiveMesh::make_triquad()
{
  MxHeap q;

  for (VertexIter p=verts_begin(); p!=verts_end(); p++)
    (*p)->set_level(0);
  
  for (EdgeIter ei=edges_begin(); ei!=edges_end(); ei++) 
    {
      (*ei)->set_mark(0);
      q.insert( (MxHeapable*)(*ei), (*ei)->length() );
    }
  
  A48Edge* e;
  while ( (e=(A48Edge*)q.extract()) != NULL )
    {
      if ( e->mark()==0 ) 
	{
	  mark_link(e);
	  if (e->hedge(0)->facet() != NULL)
	    split(e->hedge(0));
	  else
	    split(e->hedge(1));
	}
    }
  
  for (FacetIter f=facets_begin(); f!=facets_end(); f++)
    if ((*f)->level() == 0)
      split(*f);

  for (VertexIter p=verts_begin(); p!=verts_end(); p++)
    (*p)->set_level(0);
}

bool AdaptiveMesh::is_triquad()
{
  for (EdgeIter ei=edges_begin(); ei!=edges_end(); ei++) 
    {
      A48Edge* e = *ei;
      if (!e->is_bdry()) 
	{
	  bool f0s = ( der_cast(e->hedge(0)->facet())->subd_edge() == e->hedge(0) );
	  bool f1s = ( der_cast(e->hedge(1)->facet())->subd_edge() == e->hedge(1) );
	  if ((f0s && !f1s) || (!f0s && f1s))
	    return false;
	}
    }
  return true;
}
