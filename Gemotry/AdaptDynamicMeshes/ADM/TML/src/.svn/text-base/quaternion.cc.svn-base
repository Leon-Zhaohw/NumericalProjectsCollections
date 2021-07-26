#include "tml.h"

using namespace TML;

//QUATERNION ROTATION
R3 TML::RotParam(const R3& from, const R3& to, double& angle)
{
  R3 f=R3(from), t=R3(to);
  f.normalize();
  t.normalize();
  
  double inner = dot(f,t);
  if ( TML::zero(inner,1.) )
    { angle = 0.0; return R3(); }
  
  angle = acos(inner);
  R3 axis = cross(f,t);
  axis.normalize();
  return R3(axis);
}

R3 TML::Rotate(const R3& p, const R3& axis, double angle)
{
  if ( TML::zero(angle,0.) ) return R3(p);
  
  //quaternion -> q*p*q'
  double q0 = cos(angle/2);
  R3 q1 = axis * ( (double)sin(angle/2) );
  
  //a = q*p
  double a0 = - dot(q1,p);
  R3 a1 = p*q0 + cross(q1,p);
  
  //r = a*q'
  q1 *= -1;
  return ( q1*a0 + a1*q0 + cross(a1,q1) );
}
