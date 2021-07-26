#ifndef _PARAMSINC_
#define _PARAMSINC_ 1

#include "surfaces/parametric.h"

class ParamSinc : public virtual ParametricSurf
{
 protected:
  double A, B, K, X0, Y0;
  
 public:
  ParamSinc(double a, double b, double k, double x0, double y0) 
    : ParametricSurf(), A(a), B(b), K(k), X0(x0), Y0(y0) { }

  virtual ~ParamSinc() { }
  
  /*------------------------------------------------------------*/

  void normal(ADM::A48Vertex* v)
    {
      double u0, u1;
      ADM::R3 p = v->pos();
      parameterization( p, u0, u1 );
      
      if ( R(u0,u1)<FEQ_EPS )
	v->set_normal( v->compute_normal() );
      else
	v->set_normal( parametric_normal(p) );
    }
  
  void curvature(ADM::A48Vertex* v)
    {
      double u0, u1;
      ADM::R3 p = v->pos();
      parameterization( p, u0, u1 );
      
      ADM::R3 d1, d2;
      double  k1, k2, km, kg;
      
      if ( R(u0,u1)<FEQ_EPS )
	{
	  v->curvature(k1,d1,k2,d2);
	  km = 0.5*(k1+k2);
	  kg = k1*k2;
	} else {
	  parametric_curvature( v->pos(), k1, k2, d1, d2, km, kg);
	}
      
      v->set_max_curv(k1);   v->set_min_curv(k2);
      v->set_max_direc(d1);  v->set_min_direc(d2);
      v->set_gaus_curv(kg);  v->set_mean_curv(km);
    }

  /*------------------------------------------------------------*/
  
  void parameterization(const ADM::R3& p, double& u, double& v)
    { u = p.x - X0; v = -( p.z - Y0 ); }

  ADM::R3 coordinates(const double& u, const double& v)
    { double r = R(u,v); return ADM::R3(u+X0, F(r), -v+Y0); }

  ADM::R3 deriv_u(const double& u, const double& v)
    { double r = R(u,v); return ADM::R3(1., F_r(r)*R_u(u,v), 0.); }

  ADM::R3 deriv_v(const double& u, const double& v)
    { double r = R(u,v); return ADM::R3(0., F_r(r)*R_v(u,v), -1.); }

  ADM::R3 deriv_uu(const double& u, const double& v)
    { double r = R(u,v); return ADM::R3(0., pow(R_u(u,v),2)*F_rr(r) + R_uu(u,v)*F_r(r), 0.); }

  ADM::R3 deriv_uv(const double& u, const double& v)
    { double r = R(u,v); return ADM::R3(0., R_u(u,v)*R_v(u,v)*F_rr(r) + R_uv(u,v)*F_r(r), 0.); }

  ADM::R3 deriv_vv(const double& u, const double& v)
    { double r = R(u,v); return ADM::R3(0., pow(R_v(u,v),2)*F_rr(r) + R_vv(u,v)*F_r(r), 0.); }
 
  /*------------------------------------------------------------*/

  double R(const double& u, const double& v)   
    { return sqrt(u*u + v*v); }
  
  double R_u(const double& u, const double& v) 
    { double r = R(u,v); return (r<FEQ_EPS)? 0. : u/r; }
  
  double R_v(const double& u, const double& v) 
    { double r = R(u,v); return (r<FEQ_EPS)? 0. : v/r; }
  
  double R_uu(const double& u, const double& v)
    { double r = R(u,v); return (r<FEQ_EPS)? 0. : (1/r - pow(u,2)/pow(r,3)); }

  double R_uv(const double& u, const double& v)
    { double r = R(u,v); return (r<FEQ_EPS)? 0. : - u*v/pow(r,3); }

  double R_vv(const double& u, const double& v)
    { double r = R(u,v); return (r<FEQ_EPS)? 0. : (1/r - pow(v,2)/pow(r,3)); }
  
  /*------------------------------------------------------------*/
  
  double G(const double& r) 
    { return (r<FEQ_EPS)? 1. : std::sin(B*r) / (B*r); }
  
  double G_r(const double& r)
    { return (r<FEQ_EPS)? 0. : ( std::cos(B*r) - G(r) ) / r; }
  
  double G_rr(const double& r) 
    { return (r<FEQ_EPS)? 0. : ( -B*B*G(r) -2*G_r(r)/r ); }

  /*------------------------------------------------------------*/

  double F(const double& r) 
    { return A*G(r)*std::exp(-K*r*r); }
  
  double F_r(const double& r)  
    { double e = std::exp(-K*r*r); return A*(G_r(r)*e - 2.*K*r*G(r)*e); }
 
  double F_rr(const double& r)
    { 
      double e = std::exp(-K*r*r);
      return A*(G_rr(r)*e - 4.*K*r*G_r(r)*e - 2.*K*G(r)*(e - 2.*K*r*r*e));
    }
};

/*********************************************************************************/

class HeightSinc : public virtual ParamSinc
{
 private:
  bool flag;

 public:
  HeightSinc(double a, double b, double k, double x0, double y0) 
    : ParamSinc(a,b,k,x0,y0), flag(true) { }

  ~HeightSinc() { }
  
  void deform()
    {
      if (A >  0.8) flag = false;
      if (A < -0.8) flag = true;
      
      if (flag) A += 0.02;
      else      A -= 0.02;

      std::cout << A << std::endl;
    }
};

class BumpSinc : public virtual ParamSinc
{
 private:
  bool flag;

 public:
  BumpSinc(double a, double b, double k, double x0, double y0) 
    : ParamSinc(a,b,k,x0,y0), flag(true) { }
  
  ~BumpSinc() { }
  
  void deform()
    {
      if (X0 >  0.8) flag = false;
      if (X0 < -0.8) flag = true;
      
      if (flag) X0 += 0.02;
      else      X0 -= 0.02;
    }
};

#endif
