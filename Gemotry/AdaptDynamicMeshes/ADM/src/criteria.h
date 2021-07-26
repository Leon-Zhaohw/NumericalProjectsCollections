#ifndef _CRITERIA_
#define _CRITERIA_

namespace ADM
{
#define MIN_ANGLE 0.174532925199433 //  10 degrees
#define MAX_ANGLE 2.792526803190927 // 160 degrees

  // Abstract class -> adaptation criteria interface
  class AdaptCriteria
    {
    protected:
      Oracle* surf_;

    public:
      void init(Oracle* s) { surf_ = s; }

      Oracle* surf() const { return surf_; }
      
      virtual bool pathology(A48Vertex* v);
      virtual bool pathology(A48Edge* e);
      virtual bool pathology(A48Facet* f);

      virtual double ref_rank (A48Edge* e)   = 0;
      virtual double simp_rank(A48Vertex* v) = 0;
    };

  /*------------------------------------------------------------*/
  
  // Static stochastic sampling over triangles
  class StochasticSampling : public AdaptCriteria
    {
    protected:
      double density_;

    public:
      StochasticSampling(Oracle* s) 
	: density_(2.) 
	{ init(s); }
      
      virtual ~StochasticSampling() { }
      
      double density() const       { return density_;} 
      void   set_density(double d) { density_ = d; }
      
      virtual double ref_rank(A48Edge* e)    
	{ return ( pathology(e) )? -1. : sampling_ref_rank(e);  }
	
      virtual double simp_rank(A48Vertex* v) 
	{ 
	  return ( ( v->is_inbase() )? -1 : 
		   ( pathology(v) )? 1e6 : sampling_simp_rank(v) );
	}
      
    protected:
      // Sampling without learning
      double sampling_ref_rank  (A48Edge* e);
      double sampling_simp_rank(A48Vertex* v);
      
      double area(R3& v0, R3& v1, R3& v2);
      double quad_section_sampling(int lv, R3 v0, R3 v1, R3 v2);
    };
}

#endif
