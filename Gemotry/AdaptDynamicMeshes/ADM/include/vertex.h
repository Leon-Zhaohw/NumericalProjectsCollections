#ifndef _A48_VERTEX_
#define _A48_VERTEX_ 1

namespace ADM
{
  class A48Vertex : public Vertex
    {
    protected:
      int lv_; 
      
      double kmax_, kmin_, kmean_, kgaus_;
      R3 dmax_, dmin_, normal_;
      
    public:
      A48Vertex(R3 p=R3()) : Vertex(p), lv_(0), 
	kmax_(0.), kmin_(0.), kmean_(0.), kgaus_(0.),
	dmax_(), dmin_(), normal_()
	{ }
      
      virtual ~A48Vertex() { }

      /*------------------------------------------------------------*/
      
      inline int level()        { return lv_;     }
      inline double max_curv()  { return kmax_;   }
      inline double min_curv()  { return kmin_;   }
      inline double mean_curv() { return kmean_;  }
      inline double gaus_curv() { return kgaus_;  }
      inline R3 max_direc()     { return dmax_;   }
      inline R3 min_direc()     { return dmin_;   }
      inline R3 normal()        { return normal_; }
    
      inline void set_level(int l)        { lv_ = l;     }
      inline void set_max_curv(double k)  { kmax_ = k;   }
      inline void set_min_curv(double k)  { kmin_ = k;   }
      inline void set_mean_curv(double k) { kmean_ = k;  }
      inline void set_gaus_curv(double k) { kgaus_ = k;  }
      inline void set_max_direc(R3 d)     { dmax_ = d;   }
      inline void set_min_direc(R3 d)     { dmin_ = d;   }
      inline void set_normal(R3 n)        { normal_ = n; }
    
      /*------------------------------------------------------------*/
      
      bool is_inbase() 
	{ return ( level()==0 ); }
      
      bool is_weld(Facet* exclude=NULL);
    };

  // Casting TML types to A48 types
  inline A48Vertex* der_cast(Vertex* v) { return static_cast<A48Vertex*>(v); }
}

#endif
