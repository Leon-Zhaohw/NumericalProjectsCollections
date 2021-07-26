#include "adm.h"

using namespace ADM;

AdaptiveMesh::AdaptiveMesh(const char* base_mesh, Oracle* s)
  : TML::TriMesh< A48Vertex, A48Edge, A48Facet >(base_mesh),
    surf_(s), criteria_(new StochasticSampling(s)), 
    sf_(), rf_(),
    threshold_(0.008), stiffness_(-2.), hysteresis_(0.), smooth_steps_(1)
{
  if (!is_triquad()) make_triquad();
  update_attr();
  align_base_mesh();
  new_front();
  
#ifdef CHECKER
  greedy_color();
#endif
}

AdaptiveMesh::AdaptiveMesh(const char* base_mesh, AdaptCriteria* cr)
  : TML::TriMesh< A48Vertex, A48Edge, A48Facet >(base_mesh),
    surf_(cr->surf()), criteria_(cr),
    sf_(), rf_(),
    threshold_(0.008), stiffness_(-2.), hysteresis_(0.), smooth_steps_(1)
{
  if (!is_triquad()) make_triquad();
  update_attr();
  align_base_mesh();
  new_front();

#ifdef CHECKER
  greedy_color();
#endif
}

/*------------------------------------------------------------*/
      
void AdaptiveMesh::align_base_mesh()
{
  for (VertexIter vi=verts_begin(); vi!=verts_end(); vi++)
    surf_->project(*vi);
  update_attr();
}

void AdaptiveMesh::adapt()
{
  strutural_phase();
  for (int i=0; i<smooth_steps_; i++)
    smoothing_phase();
}

void AdaptiveMesh::deform()
{
  surf_->deform();
  for (VertexIter vi=verts_begin(); vi!=verts_end(); vi++)
    surf_->update(*vi);
  update_attr();
}

void AdaptiveMesh::strutural_phase()
{
  simplify_phase();
  refine_phase();
}

void AdaptiveMesh::simplify_phase()
{ 
  adapt_simplify( 1./(threshold_-hysteresis_) );   
  update_attr();
}

void AdaptiveMesh::simplify_one_step()
{ 
  adapt_simplify_step( 1./(threshold_-hysteresis_) ); 
  update_attr();
}

void AdaptiveMesh::refine_phase()
{ 
  adapt_refine( threshold_+hysteresis_ ); 
  update_attr();
}

void AdaptiveMesh::refine_one_step()
{ 
  adapt_refine_step( threshold_+hysteresis_ ); 
  update_attr();
}

void AdaptiveMesh::smoothing_phase()
{
  double w1, w2, alpha, beta;
  R3 p, n, c, offset, move;
  
  double abs_max, abs_min, delta;
  stat_curvature( abs_max, abs_min );
  delta = abs_max - abs_min;
  
  for (VertexIter vi=verts_begin(); vi!=verts_end(); vi++)
    {
      if ( (*vi)->is_bdry() ) continue;
      
      p = (*vi)->pos();
      n = (*vi)->normal();
      c = (*vi)->star_center();
      
      // computing anistropic smoothing
      offset = (c-p) - dot(c-p,n)*n;
      
      if ( TML::zero(delta,0.,1e-6) )
	{ w1 = w2 = 1.; }
      else
	{
	  alpha = ( std::fabs((*vi)->max_curv()) - abs_min ) / delta; 
	  w1 = std::exp( stiffness_ * alpha );
	  alpha = ( std::fabs((*vi)->min_curv()) - abs_min ) / delta; 
	  w2 = std::exp( stiffness_ * alpha );
	}
      
      move = 
	(w1 * dot( offset, (*vi)->max_direc() )) * (*vi)->max_direc() + 
	(w2 * dot( offset, (*vi)->min_direc() )) * (*vi)->min_direc();
      
      (*vi)->set_pos( p + move );
      surf_->project(*vi);
      
      // checking folding
      move = (*vi)->pos() - p;
      avoid_folding( *vi, move, beta );
      if ( !TML::zero(beta,1.) )
	(*vi)->set_pos( p );
    }
  
  update_attr();
}

void AdaptiveMesh::stat_curvature(double& max_val, double& min_val)
{
  double abs1, abs2;

  VertexIter vi = verts_begin();
  abs1 = std::fabs( (*vi)->max_curv() );
  abs2 = std::fabs( (*vi)->min_curv() );
  max_val = std::max( abs1, abs2 );
  min_val = std::min( abs1, abs2 );
  vi++;
  
  for ( ; vi!=verts_end(); vi++)
    {
      abs1 = std::fabs( (*vi)->max_curv() );
      abs2 = std::fabs( (*vi)->min_curv() );
      max_val = std::max( max_val, std::max( abs1, abs2 ) );
      min_val = std::min( min_val, std::min( abs1, abs2 ) );
    }
}

void AdaptiveMesh::avoid_folding(A48Vertex* v, const R3& u, double& beta)
{
  int i;
  double b, c;
  R3 p0, p1, p2, c0, c1, bisector;
  Halfedge *hv, *hf;

  beta = 1;
  for (hv=v->star_first(); hv!=NULL; hv=v->star_next(hv))
    {
      if (hv->is_bdry()) continue;

      p0 = hv->dst()->pos();
      p1 = hv->prev()->org()->pos();
      p2 = hv->org()->pos();

      c0 = cross( p1-p0, p2-p0 );
      c1 = cross( p2-p1, u );

      // checking facet
      b = dot( c0, c1 );
      c = dot( c0, c0 );
      if ( !TML::zero(b,0.) && b<0. )
	beta = std::min( beta, -c/b );
      
      // checking edges
      for (i=0, hf=hv; i<3; i++, hf=hf->next())
	{
	  if (hf->mate()->is_bdry()) continue;

	  bisector = 0.5*( hf->facet()->normal() + hf->mate()->facet()->normal() );
	  b = dot( bisector, c1 );
	  c = dot( bisector, c0 );
	  if ( !TML::zero(b,0.) && b<0. )
	    beta = std::min( beta, -c/b );
	}
    }
}

void AdaptiveMesh::update_attr()
{
  set_edge_attr();
  set_facet_attr();
  
  int i;
  VertexIter vi;
  for (vi=verts_begin(), i=0; vi!=verts_end(); vi++, i++)
    {
      (*vi)->set_attr(i);
      surf_->normal(*vi);
      surf_->curvature(*vi);
    }
}

void AdaptiveMesh::sample(A48Edge* e, A48Vertex* v)
{
  R3 m, n, n1, n2;
  m = 0.5*( e->org()->pos() + e->dst()->pos() );
  n1 = ( e->hedge(0)->is_bdry() )? R3() : e->hedge(0)->facet()->normal();
  n2 = ( e->hedge(1)->is_bdry() )? R3() : e->hedge(1)->facet()->normal();
  n = 0.5*( n1 + n2 );
  n.normalize();

  v->set_pos( m );
  surf_->project( v, n );
}

void AdaptiveMesh::sample(A48Facet* f, A48Vertex* v)
{
  R3 b, n;
  b = f->barycenter();
  n = f->normal();
  
  v->set_pos( b );
  surf_->project( v, n );
}
