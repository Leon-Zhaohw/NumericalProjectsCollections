#ifndef _PARAMSIN_
#define _PARAMSIN_ 1

#include "surfaces/parametric.h"

class ParamSin : public ParametricSurf
{
 private:
  double A, B, C, D;

 public:
  ParamSin(double a, double b, double c, double d) 
    : ParametricSurf(), A(a), B(b), C(c), D(d) { }
  
  ~ParamSin() { }
  
  /*------------------------------------------------------------*/

  void deform()
    {
      C += 5.*M_PI / 180.;
      if (C >= M_PI)      C -= 2*M_PI;
      else if (C < -M_PI) C += 2*M_PI;
    }
  
  /*------------------------------------------------------------*/
  
  void parameterization( const ADM::R3& p, double& u, double& v ) 
    { u = p.x; v = p.z; }
  
  ADM::R3 coordinates(const double& u, const double& v) 
    { return ADM::R3(u, A*std::sin(B*u+C)+D, v); }
  
  ADM::R3 deriv_u (const double& u, const double& v)           
    { return ADM::R3(1, A*B*cos(B*u+C), 0); }
  
  ADM::R3 deriv_v (const double& u, const double& v)            
    { return ADM::R3(0,0,-1); }
  
  ADM::R3 deriv_uu(const double& u, const double& v)             
    { return ADM::R3(0, -A*B*B*sin(B*u+C), 0); }
  
  ADM::R3 deriv_uv(const double& u, const double& v)            
    { return ADM::R3(); }
  
  ADM::R3 deriv_vv(const double& u, const double& v)             
    { return ADM::R3(); }
};

#endif
