#ifndef _FACET_
#define _FACET_ 1

namespace TML
{
  class Facet : public Markable, public MxHeapable
    {
    protected:
      int       id_;
      Halfedge* e_;

    public:
      Facet(Halfedge* h0, Halfedge* h1, Halfedge* h2) : id_(0), e_(NULL)
	{ reuse(h0, h1, h2); }
      virtual ~Facet() { }

      // Basic Methods
      int  id()          { return id_; }
      void set_id(int i) { id_ = i; }

      Halfedge* hedge(int k);
      Vertex*   vertex(int k);
      
      void set_hedge(int k, Halfedge* h);
      void set_vertex(int k, Vertex* v);
      
      void link_star_verts();
      
      Facet* reuse(Halfedge *e0, Halfedge *e1, Halfedge *e2);
    
      // Geometric Methods
      virtual void set_attr(int i) { set_id(i); }
      virtual double area();
      virtual R3     normal();
      virtual R3     barycenter();
      virtual R3     circumcenter();
    };

 /*-------------------------------------------------------------*/

  class Polygon
    {
    protected:
      int n_, *verts_;
      
    public:
      Polygon() : n_(0), verts_(NULL)      { }
      Polygon(int n) : n_(n), verts_(NULL) { set(n); }
      ~Polygon() { if (!verts_) delete[] verts_; }
      
      void set(int n)
	{
	  n_ = n;
	  if (verts_) delete[] verts_;
	  verts_ = new int[n_];
	  for (int i=0; i<n_; i++) verts_[i] = -1;
	}
      
      inline int size()        const { return n_; }
      inline int vertex(int i) const { return verts_[i]; }
      
      inline int& vertex(int i)        { return verts_[i]; }
      inline void vertex(int i, int v) { verts_[i] = v; }
    };
}

#endif
