#include "adm.h"

using namespace ADM;

A48Vertex* AdaptiveMesh::refine(Halfedge* e)
{
  A48Facet* f; 
  Halfedge* r[2]; 
  int n = 0;
  
  f = der_cast(e->facet());
  if ( f!=NULL && f->subd_edge()!=e) 
    r[n++] = f->subd_edge();

  f = der_cast(e->mate()->facet());
  if ( f!=NULL && f->subd_edge()!=e->mate() ) 
    r[n++] = f->subd_edge();
  
  for (int k=0; k<n; k++)
    refine(r[k]);
  
  update_ref_front( der_cast(e->edge()) );
  A48Vertex* w = split(e);
  update_ref_front(w);
  return w;
}

Halfedge* AdaptiveMesh::simplify(A48Vertex* w)
{
  int n = 0, weld_deg = (w->is_bdry())? 3 : 4;
  do {
    int lmax = w->level(); 
    Halfedge* e; 
    A48Vertex *u, *v;
    for (e=w->star_first(), n=0; e!=NULL; e=w->star_next(e), n++) 
      {
	u = der_cast(e->org());
	if (u->level() > lmax) 
	  {
	    lmax = u->level(); 
	    v = u;
	  }
      }
    if (lmax > w->level())
      simplify(v);
  } while (n > weld_deg);
  
  update_simpl_front(w);
  Halfedge* e = weld(w);
  if (e==NULL) TML::Error("weld returns NULL");
  update_simpl_front( der_cast(e->edge()) );
  return e;
}

void AdaptiveMesh::adapt_refine(double t)
{
  for (EdgeIter ei=edges_begin(); ei!=edges_end(); ei++)
    if ((*ei)->is_in_heap())
      rf_.update((MxHeapable*)(*ei), criteria_->ref_rank(*ei));
  
  A48Edge* e;
  while ( (e=(A48Edge*)rf_.extract()) != NULL ) 
    {
      if (e->heap_key() < t) 
	{
	  rf_.insert((MxHeapable*)e, criteria_->ref_rank(e));
	  return;
	}
      refine(e->hedge(0));
    }
} 

void AdaptiveMesh::adapt_simplify(double t)
{
  for (VertexIter vi=verts_begin(); vi!=verts_end(); vi++) 
    if ((*vi)->is_in_heap())
      sf_.update((MxHeapable*)(*vi), criteria_->simp_rank(*vi));

  A48Vertex* v;
  while ( (v=(A48Vertex*)sf_.extract()) != NULL ) 
    {
      if (v->heap_key() < t) 
	{
	  sf_.insert((MxHeapable*)v, criteria_->simp_rank(v));
	  return;
	}
      simplify(v);
    }
}

void AdaptiveMesh::adapt_refine_step(double t)
{
  for (EdgeIter ei=edges_begin(); ei!= edges_end(); ei++)
    if ((*ei)->is_in_heap())
      rf_.update((MxHeapable*)(*ei), criteria_->ref_rank(*ei));

  A48Edge* e = (A48Edge*)rf_.extract();
  if (e!=NULL)
    if (e->heap_key() < t)
      {
	rf_.insert((MxHeapable*)e, criteria_->ref_rank(e));
	std::cout << "--END refinement" << std::endl;
      }
    else
      refine(e->hedge(0));
}

void AdaptiveMesh::adapt_simplify_step(double t)
{
  for (VertexIter vi=verts_begin(); vi!=verts_end(); vi++)
    if ((*vi)->is_in_heap())
      sf_.update((MxHeapable*)(*vi), criteria_->simp_rank(*vi));
  
  A48Vertex* v = (A48Vertex*)sf_.extract();
  if (v!=NULL)
    if (v->heap_key() < t)
      sf_.insert((MxHeapable*)v, criteria_->simp_rank(v));
    else
      simplify(v);
}
