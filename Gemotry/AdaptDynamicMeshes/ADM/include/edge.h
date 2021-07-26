#ifndef _A48_EDGE_
#define _A48_EDGE_ 1

namespace ADM
{
  Halfedge* halfedge_reuse(Halfedge* h, Vertex* v0, Vertex* v1);
  
  bool halfedge_is_subdiv(Halfedge* h);

  /*------------------------------------------------------------*/
  
  class A48Edge : public Edge
    {
    public:
      A48Edge(Vertex* p0, Vertex* p1) 
	: Edge(p0,p1) { }
      
      virtual ~A48Edge() { }
      
      bool is_subdiv()
	{ return ( halfedge_is_subdiv( &(h_[0]) ) || halfedge_is_subdiv( &(h_[1]) ) ); }
      
      int level();
    };
  
  // Casting TML types to A48 types
  inline A48Edge* der_cast(Edge* e)   { return static_cast<A48Edge*>(e);   }
}

#endif
