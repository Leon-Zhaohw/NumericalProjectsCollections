#ifndef _PARAMETRIC_SURF_
#define _PARAMETRIC_SURF_ 1

#include <adm.h>

class ParametricSurf : public ADM::Oracle
{
 protected:
  ADM::Matrix Jac, Wein, EigVec;
  ADM::Array  EigVal;

 public:
  ParametricSurf() 
    : Jac(3,2), Wein(2,2), EigVec(2,2), EigVal(2)
    { }
  
  virtual ~ParametricSurf() { }
  
  /*------------------------------------------------------------*/
  
  inline void update(ADM::A48Vertex* v) { project(v); }
  
  inline void project(ADM::A48Vertex* v, ADM::R3 n=ADM::R3())
    { v->set_pos( position(v->pos()) ); }
  
  inline double error(const ADM::R3& p)
    { return (position(p) - p).length(); }
  
  /*------------------------------------------------------------*/
  
  virtual void normal(ADM::A48Vertex* v)
    { v->set_normal( parametric_normal(v->pos()) ); }
  
  virtual void curvature(ADM::A48Vertex* v)
    {
      ADM::R3 d1, d2;
      double  k1, k2, km, kg;
      
      parametric_curvature( v->pos(), k1, k2, d1, d2, km, kg);

      v->set_max_curv(k1);   v->set_min_curv(k2);
      v->set_max_direc(d1);  v->set_min_direc(d2);
      v->set_gaus_curv(kg);  v->set_mean_curv(km); 
    }
  
  /*------------------------------------------------------------*/
  
  ADM::R3 position(const ADM::R3& p)
    {
      double u, v;
      parameterization( p, u, v );
      return coordinates( u, v );
    }
  
  ADM::R3 parametric_normal(const ADM::R3& p)
    {
      double u, v;
      parameterization( p, u, v );

      ADM::R3 d_u = deriv_u(u,v);
      ADM::R3 d_v = deriv_v(u,v);
      ADM::R3 n = cross( d_u, d_v );
      n.normalize();

      return n;
    }
  
  void parametric_curvature(const ADM::R3& p,
			    double& k1, double& k2,
			    ADM::R3& d1, ADM::R3& d2,
			    double& km, double& kg)
    {
      double u, v;
      parameterization( p, u, v );
      
      // Jacobian matrix
      ADM::R3 d_u = deriv_u(u,v), d_v = deriv_v(u,v);
      Jac[0][0] = d_u.x;  Jac[0][1] = d_v.x;
      Jac[1][0] = d_u.y;  Jac[1][1] = d_v.y;
      Jac[2][0] = d_u.z;  Jac[2][1] = d_v.z;
      
      // First fundamental form
      double E = dot( d_u, d_u );
      double F = dot( d_u, d_v );
      double G = dot( d_v, d_v );
      
      // Second fundamental form
      ADM::R3 d_uu=deriv_uu(u,v), d_uv=deriv_uv(u,v), d_vv=deriv_vv(u,v);
      ADM::R3 n = cross( d_u,d_v ); n.normalize();
      double e = dot( n, d_uu );
      double f = dot( n, d_uv );
      double g = dot( n, d_vv );
      
      double den = E*G - F*F;
      kg = (e*g - f*f) / den;              // gaussian curv
      km = -0.5*(e*G - 2*f*F + g*E) / den; // mean curv
      
      // dN matrix
      Wein[0][0] = (f*F - e*G)/den; Wein[0][1] = (g*F - f*G)/den;
      Wein[1][0] = (e*F - f*E)/den; Wein[1][1] = (f*F - g*E)/den;
      
      ADM::Eigen eig(Wein);
      eig.getRealEigenvalues(EigVal);
      eig.getV(EigVec);
      
      int vmax, vmin;
      if (EigVal[0]>=EigVal[1]) { vmax = 0; vmin = 1; }
      else                      { vmax = 1; vmin = 0; }
      
      k1 = EigVal[vmax]; // max curv
      k2 = EigVal[vmin]; // min curv
      
      ADM::Matrix P(3,2);
      P = TNT::matmult( Jac, EigVec );
      
      // max direc
      d1 = ADM::R3( P[0][vmax], P[1][vmax], P[2][vmax] );
      d1.normalize();
      
      // min direc
      d2 = ADM::R3( P[0][vmin], P[1][vmin], P[2][vmin] );
      d2.normalize();
    }

  /*------------------------------------------------------------*/
  
  virtual void parameterization( const ADM::R3& p, double& u, double& v ) = 0;
  virtual ADM::R3 coordinates(const double& u, const double& v) = 0;
  virtual ADM::R3 deriv_u (const double& u, const double& v) = 0;
  virtual ADM::R3 deriv_v (const double& u, const double& v) = 0;
  virtual ADM::R3 deriv_uu(const double& u, const double& v) = 0;
  virtual ADM::R3 deriv_uv(const double& u, const double& v) = 0;
  virtual ADM::R3 deriv_vv(const double& u, const double& v) = 0;  
};

#endif
