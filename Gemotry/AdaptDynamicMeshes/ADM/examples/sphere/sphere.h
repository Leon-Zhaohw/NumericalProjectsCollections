#ifndef _SPHERE_ 
#define _SPHERE_ 1

#include <adm.h>

class Sphere : public ADM::Oracle
{
 protected: 
  double R;

 public:
  Sphere(double r) : R(r) { }
  
  virtual ~Sphere() { }
  
  /*------------------------------------------------------------*/

  void deform() 
    { if (R<5.) R += 0.5; }
  
  void update(ADM::A48Vertex* v)
    { project(v); }

  void project(ADM::A48Vertex* v, ADM::R3 n=ADM::R3())
    { v->set_pos( position(v->pos()) ); }
  
  double error(const ADM::R3& p)
    { return (position(p) - p).length(); }

  /*------------------------------------------------------------*/
  
  void normal(ADM::A48Vertex* v)
    { v->set_normal( normal(v->pos()) ); }
  
  void curvature(ADM::A48Vertex* v)
    {
      ADM::R3 d1, d2;
      double  k1, k2;
      
      curvature( v->pos(), k1, k2, d1, d2 );

      v->set_max_curv(k1);     v->set_min_curv(k2);
      v->set_max_direc(d1);    v->set_min_direc(d2);
      v->set_gaus_curv(k1*k2); v->set_mean_curv(0.5*(k1+k2));
    }
  
  /*------------------------------------------------------------*/
  
  ADM::R3 position(const ADM::R3& p) 
    {
      ADM::R3 a = p;
      a.normalize();
      return R*a;
    }
  
  ADM::R3 normal(const ADM::R3& p)
    {
      ADM::R3 n = p;
      n.normalize();
      return n;
    }
  
  void curvature(const ADM::R3& p, double& k1, double& k2, ADM::R3& d1, ADM::R3& d2)
    {
      k1 = k2 = 1 / R;
      
      //Checking extremes
      if  ( (std::fabs(p.x)<FEQ_EPS && std::fabs(p.y)<FEQ_EPS) || 
	    (std::fabs(p.x)<FEQ_EPS && std::fabs(p.z)<FEQ_EPS) ) 
	d1 = ADM::R3(1.,0.,0.);
      else 
	if ( std::fabs(p.y)<FEQ_EPS && std::fabs(p.z)<FEQ_EPS )  
	  d1 = ADM::R3(0.,1.,0.);
	else
	  {
	    ADM::R3 t = ADM::R3(-p.y,p.x,0);
	    t.normalize();
	    d1 = t;
	  }
      d2 = cross( normal(p), d1 );
    }
};

#endif
