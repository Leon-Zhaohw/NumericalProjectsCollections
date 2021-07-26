#ifndef _A48_FACET_
#define _A48_FACET_

namespace ADM
{
  class A48Facet : public Facet
    {
    public:
      A48Facet(Halfedge* h0, Halfedge* h1, Halfedge* h2) 
	: Facet(h0,h1,h2) { }
      
      virtual ~A48Facet() { }
      
      Halfedge*  subd_edge() { return hedge(0);  } 
      Vertex* weld_vertex()  { return vertex(0); }
      
      bool is_inbase();
      int  level();
    };

  // Casting TML types to A48 types
  inline A48Facet*  der_cast(Facet* f)  { return static_cast<A48Facet*>(f);  }
}

#endif
