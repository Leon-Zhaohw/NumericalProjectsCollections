#include "tml.h"

using namespace TML;

Halfedge* Halfedge::mate() 
{ return (this == e_->hedge(0))? e_->hedge(1) : e_->hedge(0); }

R3 Halfedge::direction()
{ 
  R3 d = dst()->pos() - org()->pos(); 
  d.normalize();
  return d;
}

// angle at org() 
double Halfedge::angle()
{
  R3 v1 = direction();
  R3 v2 = prev()->mate()->direction();
  R3 n  = cross(v1,v2);
  return std::atan2( n.length(), dot(v1,v2) );
}

// cotg at angle opposite to halfedge
double Halfedge::cotangent()
{
  if (is_bdry()) return 0.;
  R3 v1 = prev()->direction();
  R3 v2 = next()->mate()->direction();
  R3 n  = cross(v1,v2);
  return dot(v1,v2) / n.length();
}

// check if angle at org() is obtuse
bool Halfedge::is_obtuse()
{ return ( dot( direction(), prev()->mate()->direction() ) < 0. ); }

// curvature at halfedge from org() to dst()
double Halfedge::curvature()
{
  R3 n = dst()->normal();
  R3 v = dst()->pos() - org()->pos();
  return ( 2.*dot(v,n) / v.length2() );
}

/*-------------------------------------------------------------*/

double Edge::length2()
{ return (org()->pos()-dst()->pos()).length2(); }

double Edge::length()
{ return (org()->pos()-dst()->pos()).length(); }

double Edge::dihedral()
{
  if ( is_bdry() ) throw Error("Edge::dihedral");
  R3 n1 = h_[0].facet()->normal();
  R3 n2 = h_[1].facet()->normal();
  R3 v  = cross(n1,n2);
  R3 d  = h_[0].direction();
  return std::atan2( dot(v,d), dot(n1,n2) );
}

double Edge::cotangent()
{ return (h_[0].cotangent() + h_[1].cotangent()); }

R3 Edge::middle()
{ return 0.5*( org()->pos() + dst()->pos() ); }
