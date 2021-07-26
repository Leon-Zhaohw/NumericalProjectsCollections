#ifndef _EDGE_
#define _EDGE_ 1


namespace TML
{
  class Halfedge : public Markable, public MxHeapable
    {
    protected:
      Vertex*   o_;
      Edge*     e_;
      Facet*    f_;
      Halfedge* n_;
      
    public:
      Halfedge() : o_(NULL), e_(NULL), f_(NULL), n_(NULL) { }
      virtual ~Halfedge() { }
      
      // Basic Methods
      Vertex* org() { return o_; };
      Vertex* dst() { return mate()->org(); };
      
      Edge*  edge() { return e_; };
      
      Facet* facet() { return f_; };
      
      Halfedge* prev() { return (n_)? n_->next() : NULL; };
      Halfedge* next() { return n_; };
      Halfedge* mate();
      
      bool is_bdry() { return (f_==NULL); }

      void set_facet(Facet* f)   { f_ = f; };
      void set_org(Vertex* v)    { o_ = v; };
      void set_next(Halfedge* h) { n_ = h; };
      void set_edge(Edge* e)     { e_ = e; };

      // Geometric Methods
      R3     direction();
      double angle();
      double cotangent();
      bool   is_obtuse();
      double curvature();
    };
  
  /*-------------------------------------------------------------*/
  
  class Edge : public Markable, public MxHeapable
    {
    protected:
      Halfedge h_[2];
      
    public:
      Edge(Vertex* p0, Vertex* p1) 
	{
	  h_[0].set_org(p0);     h_[1].set_org(p1);
	  h_[0].set_facet(NULL); h_[1].set_facet(NULL);;
	  h_[0].set_next(NULL);  h_[1].set_next(NULL);
	  h_[0].set_edge(this);  h_[1].set_edge(this);
	}
      virtual ~Edge() { }
      
      // Basic Methods
      Halfedge* hedge(int i)
	{
	  switch (i) 
	    {
	    case 0: return &(h_[0]);
	    case 1: return &(h_[1]);
	    };
	  throw Error("Edge::hedge(i)");
	}

      Vertex* org() { return h_[0].org(); }
      Vertex* dst() { return h_[1].org(); }
      
      bool is_bdry() { return (h_[0].is_bdry() || h_[1].is_bdry()); }
      
      // Geometric Methods
      virtual void set_attr() { }
      
      double length2();
      double length();
      double dihedral();
      double cotangent();
      R3     middle();
    };
}

#endif
