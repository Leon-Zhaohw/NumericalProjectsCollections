#include "adm.h"

using namespace ADM;

bool AdaptCriteria::pathology(A48Vertex* v)
{
  for (Halfedge* h=v->star_first(); h!=NULL; h=v->star_next(h))
    {
      if ( h->is_bdry() ) continue;
      if ( pathology( der_cast(h->facet()) ) ) return true;
    }
  return false;
}

bool AdaptCriteria::pathology(A48Edge* e)
{
  for (int i=0; i<2; i++)
    {
      if ( e->hedge(i)->is_bdry() ) continue;
      if ( pathology( der_cast(e->hedge(i)->facet()) ) ) return true;
    }
  return false;
}

bool AdaptCriteria::pathology(A48Facet* f)
{
  double angle;
  for (int i=0; i<3; i++)
    {
      angle = f->hedge(i)->angle();
      if ( angle < 0. )
	{ std::cout << "Angle " << angle << std::endl; TML::Error("Angle"); }
      if ( angle<MIN_ANGLE || angle>MAX_ANGLE )
	return true;
    }
  return false;
}

/**************************************************************************/

double StochasticSampling::quad_section_sampling(int lv, R3 v0, R3 v1, R3 v2)
{
  if (lv==0) { R3 b = (v0+v1+v2)/3.; return surf_->error(b); }

  R3 subd[3];
  subd[0] = 0.5*(v1 + v2);
  subd[1] = 0.5*(v2 + v0);
  subd[2] = 0.5*(v0 + v1);
  lv--;
  return ( quad_section_sampling( lv, v0,      subd[2], subd[1]) +
           quad_section_sampling( lv, v1,      subd[0], subd[2]) +
           quad_section_sampling( lv, v2,      subd[1], subd[0]) +
           quad_section_sampling( lv, subd[0], subd[1], subd[2] ) );
}

double StochasticSampling::area(R3& v0, R3& v1, R3& v2)
{
  R3 d0 = v1 - v0;
  R3 d1 = v2 - v1;
  R3 a  = cross(d0,d1);
  return 0.5*a.length();
}

double StochasticSampling::sampling_ref_rank(A48Edge* e)
{
  A48Facet* f;
  double lv, error=0., ns=0.;
  
  if (e->hedge(0)->facet())
    {
      f = der_cast( e->hedge(0)->facet() );
      lv = std::floor( density_*f->area() );
      ns += std::pow( 4, lv );
      error += quad_section_sampling( int(lv), 
				      f->vertex(0)->pos(),
				      f->vertex(1)->pos(),
				      f->vertex(2)->pos() );
    }

  if (e->hedge(1)->facet())
    {
      f = der_cast( e->hedge(1)->facet() );
      lv = std::floor( density_*f->area() );
      ns += std::pow( 4, lv );
      error += quad_section_sampling( int(lv), 
				      f->vertex(0)->pos(),
				      f->vertex(1)->pos(),
				      f->vertex(2)->pos() );
    }
  
  if (ns==0.) TML::Error("StochasticSampling::sampling_ref_rank ns");
  
  return (error / ns);
}

double StochasticSampling::sampling_simp_rank(A48Vertex* v)
{
  int num = 0;
  R3 block[4];
  
  Halfedge* e;
  for (e=v->star_first(); e!=NULL; e=v->star_next(e))
    {
      if (der_cast(e->org())->level() < v->level())
        {
          block[num] = e->org()->pos();
          num++;
        }
    }
  
  double lv, error, ns;
  if (num==3)
    {
      lv = std::floor( density_*area(block[0],block[1],block[2]) );
      ns = std::pow( 4, lv );
      error = quad_section_sampling( int(lv), block[0], block[1], block[2] );
      return (ns / error);
    }

  if (num==4)
    {
      lv = std::floor( density_*area(block[0],block[2],block[1]) );
      ns = std::pow( 4, lv );
      error = quad_section_sampling( int(lv), block[0], block[2], block[1] );

      lv = std::floor( density_*area(block[0],block[2],block[3]) );
      ns += std::pow( 4, lv );
      error += quad_section_sampling( int(lv), block[0], block[2], block[3] );

      return (ns / error);
    }
  
  if (num>4) TML::Error("StochasticSampling::sampling_simp_rank num");

  return -1.;
}

