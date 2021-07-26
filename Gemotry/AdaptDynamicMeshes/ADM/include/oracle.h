#ifndef _ORACLE_
#define _ORACLE_ 1

namespace ADM
{
  // Abstract class -> interface to describe surface.
  class Oracle
    {
    public:
      // deform oracle's parameters
      virtual void deform() = 0; 
      
      // update vertex after deformation
      virtual void update(A48Vertex* v) = 0; 

      // project vertex back to surface
      virtual void project(A48Vertex* v, R3 dir=R3()) = 0; 

      // compute approximation error at p
      virtual double error(const R3& p) = 0;

      /*------------------------------------------------------------*/
      
      virtual void normal(A48Vertex* v)
	{ v->set_normal( v->compute_normal() ); }

      virtual void curvature(A48Vertex* v)
	{
	  double k1, k2;
	  R3     d1, d2;
	  v->curvature(k1,d1,k2,d2);
	  
	  v->set_max_curv(k1);     v->set_min_curv(k2);
	  v->set_max_direc(d1);    v->set_min_direc(d2);
	  v->set_gaus_curv(k1*k2); v->set_mean_curv(0.5*(k1+k2));
	}
    };
}

#endif
