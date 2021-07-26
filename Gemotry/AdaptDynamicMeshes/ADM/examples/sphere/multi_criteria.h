#ifndef _MULTI_CRITERIA_
#define _MULTI_CRITERIA_ 1

#include <adm.h>

class MultiCriteria : public ADM::StochasticSampling
{
 protected:
  double min_len_;

 public:
  MultiCriteria( ADM::Oracle* s, double min_len )
    : ADM::StochasticSampling(s), min_len_(min_len) { }
  
  virtual ~MultiCriteria() { }
  
  double ref_rank(ADM::A48Edge* e)
    {
      ADM::R3 o = e->org()->pos();
      ADM::R3 d = e->dst()->pos();

      if (o.x>0. || d.x>0.)
	return ( e->length() > min_len_ )? 1000 : -1;
      else
	return sampling_ref_rank(e);
    }

  double simp_rank(ADM::A48Vertex* v)
    {
      if ( v->is_inbase() ) return -1;
      
      //getting basic block
      ADM::R3 block[4];
      int num = 0;
      for (ADM::Halfedge* e=v->star_first(); e!=NULL; e=v->star_next(e))
	{
	  if ( ADM::der_cast(e->org())->level() < v->level())
	    {
	      block[num] = e->org()->pos();
	      num++;
	    }
	}
      
      ADM::R3 o = block[0];
      ADM::R3 d = block[2];
      
      if (o.x>0. || d.x>0.)
	return ( (o-d).length()<min_len_ )? 1000 : -1;
      else
	return sampling_simp_rank(v);
    }
};

#endif
