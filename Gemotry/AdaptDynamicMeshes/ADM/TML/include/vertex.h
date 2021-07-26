#ifndef _VERTEX_
#define _VERTEX_ 1

namespace TML
{
  class Vertex : public Markable, public MxHeapable
    {
    protected:
      int       id_;
      R3        g_;
      Halfedge* s_;
   
    public:
      Vertex(R3 p=R3()) : id_(0), g_(p), s_(NULL) { }
      virtual ~Vertex() { }

      // Basic Methods
      bool is_bdry() 
	{ 
	  if (s_==NULL) throw Error("Vertex::is_brdy: Isolated vertex");
	  return s_->edge()->is_bdry(); 
	}
      
      Halfedge* star_first() const 
	{ 
	  if (s_==NULL) throw Error("Vertex::star_first: Isolated vertex");
	  return s_; 
	}
      Halfedge* star_next(Halfedge* h) const 
	{
	  if (h->is_bdry()) return NULL; // other side of boundary
	  else { Halfedge* n = h->next()->mate(); return (n == s_)? NULL : n; }
	}
      inline int degree()
	{
	  int n = 0;
	  for (Halfedge* h=star_first(); h!=NULL; h=star_next(h)) n++;
	  return n;
	}
      
      void set_star(Halfedge* h) { s_ = h; }
      
      int id()  { return id_; }
      R3  pos() { return g_;  }
      void set_id(int i) { id_ = i; }
      void set_pos(R3 p) { g_ = p;  }
      
      // Geometric Methods
      virtual void set_attr(int i) { set_id(i); }
      
      double min_edge_length();
      double avg_edge_length();
      R3     star_center();

      double compute_area();
      R3     compute_normal();
      double compute_mean_curvature();
      double compute_gaus_curvature();

      virtual double area() { return compute_area(); }
      virtual R3 normal()   { return compute_normal(); }
      virtual double mean_curvature() { return compute_mean_curvature(); }
      virtual double gaus_curvature() { return compute_gaus_curvature(); }
      
      void curvature(double& k1, R3& d1, double& k2, R3& d2);
      
    protected:
      void add_tensor(Matrix& T, const double& coeff, const R3& e);
      void reorder(Array& E, Matrix& V, int* indices);
    };
}

#endif
