#include "adm.h"
  
using namespace ADM;

void AdaptiveMesh::new_front()
{
  while (rf_.extract() != NULL) ;
  
  for (EdgeIter ei= edges_begin(); ei!=edges_end(); ei++) 
    {
      if ((*ei)->is_subdiv())
	rf_.insert((MxHeapable*)(*ei), criteria_->ref_rank(*ei));
    }
  
  while (sf_.extract() != NULL) ;
  
  for (VertexIter w=verts_begin(); w!=verts_end(); w++)
    if ((*w)->is_weld())
      sf_.insert((MxHeapable*)(*w), criteria_->simp_rank(*w));
}

void AdaptiveMesh::update_ref_front(A48Edge* e)
{
  rf_.remove(e);

  if (e->hedge(0)->next() != NULL) 
    {
      A48Facet*  f = der_cast(e->hedge(0)->facet());
      A48Vertex* v = der_cast(e->hedge(0)->next()->dst());
      if ( !v->is_weld( base_cast(f) ) ) 
	sf_.remove(v);
    }
  
  if (e->hedge(1)->next() != NULL) 
    {
      A48Facet*  f = der_cast(e->hedge(1)->facet());
      A48Vertex* v = der_cast(e->hedge(1)->next()->dst());
      if ( !v->is_weld( base_cast(f) ) )
	sf_.remove(v);
    }
}

void AdaptiveMesh::update_ref_front(A48Vertex* v)
{
  sf_.insert((MxHeapable*)v, criteria_->simp_rank(v));
  for (Halfedge* s=v->star_first(); s!=NULL; s=v->star_next(s))
    if (s->prev() != NULL) 
      {
	A48Edge* e = der_cast(s->prev()->edge());
	if (e->is_subdiv()) 
	  rf_.insert((MxHeapable*)e, criteria_->ref_rank(e));
      }
}

void AdaptiveMesh::update_simpl_front(A48Vertex* v)
{
  sf_.remove(v);
  for (Halfedge* s=v->star_first(); s!=NULL; s=v->star_next(s))
    if (s->prev() != NULL) 
      {
	Halfedge* m = s->prev()->mate();
	if ( !halfedge_is_subdiv(m) ) 
	  {
	    A48Edge* e = der_cast(s->prev()->edge());
	    rf_.remove(e);
	  }
      }
}

void AdaptiveMesh::update_simpl_front(A48Edge* e)
{
  rf_.insert((MxHeapable*)e, criteria_->ref_rank(e));

  if (e->hedge(0)->next() != NULL) 
    {
      A48Vertex* w = der_cast(e->hedge(0)->next()->dst());
      if (w->is_weld())
	sf_.insert((MxHeapable*)w, criteria_->simp_rank(w));
    }

  if (e->hedge(1)->next() != NULL) 
    {
      A48Vertex* w = der_cast(e->hedge(1)->next()->dst());
      if (w->is_weld())
	sf_.insert((MxHeapable*)w, criteria_->simp_rank(w));
    }
}
