#include "adm.h"

using namespace ADM;

Halfedge* ADM::halfedge_reuse(Halfedge* h, Vertex* v0, Vertex* v1)
{
  h->set_org(v0);     h->mate()->set_org(v1);
  h->set_next(NULL);  h->mate()->set_next(NULL);
  h->set_facet(NULL); h->mate()->set_facet(NULL);
  return h;
}

bool ADM::halfedge_is_subdiv(Halfedge* h)
{ return ( !h->is_bdry() && der_cast(h->facet())->subd_edge()==h ); }

/*------------------------------------------------------------*/

int A48Edge::level()
{ return std::max( der_cast(h_[0].org())->level(), der_cast(h_[1].org())->level() ); }

/*------------------------------------------------------------*/

bool A48Vertex::is_weld(Facet* exclude)
{
  int n = 0, k = 0;
  for (Halfedge* e=star_first(); e!=NULL; e=star_next(e), k++) 
    {
      A48Facet* f = der_cast(e->facet());
      if ( f!=NULL && der_cast(f->weld_vertex())==this && f!=der_cast(exclude)) n++;
    }
  return (n > 0 && k > 2);
}     

/*------------------------------------------------------------*/

bool A48Facet::is_inbase()
{ 
  return ( der_cast(vertex(0))->level()==0 && 
	   der_cast(vertex(1))->level()==0 && 
	   der_cast(vertex(2))->level()==0 ); 
}

int A48Facet::level()
{ 
  return std::max( der_cast(vertex(0))->level(), 
		   std::max( der_cast(vertex(1))->level(), 
			     der_cast(vertex(2))->level() ) ); 
}
