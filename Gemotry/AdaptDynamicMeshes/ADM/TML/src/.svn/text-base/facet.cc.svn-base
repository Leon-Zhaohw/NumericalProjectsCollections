#include "tml.h"

using namespace TML;

#define NEXT3(k) ((k < 2)? (k+1) : 0) // Returns next half-edge index.
#define PREV3(k) ((k > 0)? (k-1) : 2) // Returns previous half-edge index.

///////////////////
// Basic Methods //
///////////////////
Halfedge* Facet::hedge(int k)
{
  switch (k) {
  case 0: return e_;
  case 1: return e_->next();
  case 2: return e_->next()->next();
  }
  throw Error("Facet::hedge index");
}

Vertex* Facet::vertex(int k)
{ return hedge(PREV3(k))->org(); }

void Facet::set_hedge(int k, Halfedge* h)
{
  Halfedge *n = hedge(NEXT3(k));
  Halfedge *p = hedge(PREV3(k));
  h->set_next(n);
  p->set_next(h);
  if (k == 0) e_ = h;
}

void Facet::set_vertex(int k, Vertex* v)
{ hedge(PREV3(k))->set_org(v); }

void Facet::link_star_verts()
{
  for (int k=0; k<3; k++)
    vertex(k)->set_star( hedge(NEXT3(k)) );
}

Facet* Facet::reuse(Halfedge* e0, Halfedge* e1, Halfedge* e2)
{
  e_ = e0;
  e0->set_next(e1);    e1->set_next(e2);    e2->set_next(e0);
  e0->set_facet(this); e1->set_facet(this); e2->set_facet(this);
  return this;
}

///////////////////////
// Geometric Methods //
///////////////////////
double Facet::area()
{
  R3 a = cross( vertex(1)->pos() - vertex(0)->pos(),
		vertex(2)->pos() - vertex(1)->pos() );
  return 0.5*a.length();
}

R3 Facet::normal()
{
  R3 n = cross( vertex(1)->pos() - vertex(0)->pos(),
		vertex(2)->pos() - vertex(1)->pos() );
  n.normalize();
  return n;
}

R3 Facet::barycenter()
{ return (vertex(0)->pos() + vertex(1)->pos() + vertex(2)->pos()) / 3.; }

R3 Facet::circumcenter()
{
  double a = ( vertex(2)->pos() - vertex(1)->pos() ).length2();
  double b = ( vertex(0)->pos() - vertex(2)->pos() ).length2();
  double c = ( vertex(1)->pos() - vertex(0)->pos() ).length2();

  double u = a*(-a + b + c);
  double v = b*( a - b + c);
  double w = c*( a + b - c);
  double t = u + v + w;

  u/=t; v/=t; w/=t;
  return ( u*vertex(0)->pos() + v*vertex(1)->pos() + w*vertex(2)->pos() );
}
