#ifndef _LEVEL_SET_
#define _LEVEL_SET_ 1

#include <adm.h>
#include <cpls3D.h>

class LevelSet : public ADM::Oracle
{
 protected:
  CPLS3D::PLS3D*     pls_;

 public:
  LevelSet( CPLS3D::PLS3D* pls ) : 
    pls_(pls)
    { }
  
  virtual ~LevelSet() 
    { }
  
  /*------------------------------------------------------------*/
  
  void deform()
    { 
      pls_->Update(); 
      std::cout << "deformation done " << pls_->getVelocity()->time() << std::endl;
    }

  void update(ADM::A48Vertex* v)
    { 
      // Offsetting vertex
      ADM::R3 p, offset;
      p      = v->pos();
      offset = RK4( p, pls_->getDT() );
      
      // Align new vertex position with level-set.
      p = project_point( p+offset, ADM::R3() );
      v->set_pos(p);
    }
  
  void project( ADM::A48Vertex* v, ADM::R3 dir )
    { v->set_pos( project_point( v->pos(), dir ) ); }
  
  double error(const ADM::R3& p)
    { return std::fabs( pls_->getLevelSet()->get( p.x, p.y, p.z ) ); }
  
  /*------------------------------------------------------------*/
  
  ADM::R3 project_point(ADM::R3 pos, ADM::R3 dir)
    {
      bool flag = ( dir==ADM::R3() )? true : false;
      
      ADM::R3 n, p = pos;
      
      CPLS3D::Grid3D* grid = &( pls_->getLevelSet()->Phi() );
      CPLS3D::Vec3D d, q;      
      
      double phi, alpha=0.5; // alpha < sqrt(3)/3
      while ( !TML::zero( (phi=grid->get(p.x,p.y,p.z)), 0., 1e-3 ) )
	{
	  if (flag)
	    {
	      q = CPLS3D::Vec3D( p.x, p.y, p.z );
	      grid->Normal( q, d );
	      n = - alpha * phi * ADM::R3(d[0],d[1],d[2]);
	    } else {
	      n = - alpha * phi * dir;
	    }
	  p += n;
	}
      
      return p;
    }
  
  /*------------------------------------------------------------*/
  
  // RK4 integration scheme
  ADM::R3 RK4(ADM::R3 p, double timestep)
    {
      CPLS3D::Vec3D q, y, k1, k2, k3, k4;
      
      q = CPLS3D::Vec3D( p.x, p.y, p.z );
      
      y = q;
      pls_->GetVelocity(y, k1);
      
      y = q + 0.5*timestep*k1;      
      pls_->GetVelocity(y, k2);

      y = q + 0.5*timestep*k2;
      pls_->GetVelocity(y, k3);

      y = q + timestep*k3;
      pls_->GetVelocity(y, k4);
      
      q = timestep*( k1 + 2*k2 + 2*k3 + k4 )/6;
      return ( ADM::R3(q[0],q[1],q[2]) );
    }
  
  /*------------------------------------------------------------*/
  
  /*
    void update_beta(ADM::A48Facet* face, ADM::R3& q)
    {
    ADM::R3 p10, p20, u10, u20, c0, c1, c2;
    double root1, root2;
    
    p10 = face->vertex(1)->pos() - face->vertex(0)->pos();
    p20 = face->vertex(2)->pos() - face->vertex(0)->pos();
    
    u10 = offset_[face->vertex(1)->id()] - offset_[face->vertex(0)->id()];
    u20 = offset_[face->vertex(2)->id()] - offset_[face->vertex(0)->id()];
    
    c0 = TML::cross( p10, p20 );
    c1 = TML::cross( p10, u20 ) - TML::cross( p20, u10 );
    c2 = TML::cross( u10, u20 );
    
    quadratic( TML::dot(q,c2), TML::dot(q,c1), TML::dot(q,c0), root1, root2 );
    
    // New beta is the smaller positive root found less than 1.
    if (root1<0.) root1 = 1.;
    if (root2<0.) root2 = 1.;
    beta_ = std::min( beta_, std::min(root1,root2) );
    }
    
    // Finding beta for facets
    void facet_check()
    {
    ADM::R3 p10, p20, c0;
    
    ADM::AdaptiveMesh::FacetIter fi;
    for (fi=mesh_->facets_begin(); fi!=mesh_->facets_end(); fi++)
    {
    p10 = (*fi)->vertex(1)->pos() - (*fi)->vertex(0)->pos();
    p20 = (*fi)->vertex(2)->pos() - (*fi)->vertex(0)->pos();
    c0 = TML::cross( p10, p20 );
    update_beta( *fi, c0 );
    }
    }
    
    // Finding beta for edges
    void edge_check()
    {
    ADM::R3 bisector;
    
    ADM::AdaptiveMesh::EdgeIter ei;
    for (ei=mesh_->edges_begin(); ei!=mesh_->edges_end(); ei++)
    {
    if ( (*ei)->is_bdry() ) continue;
    
    bisector = (*ei)->hedge(0)->facet()->normal() + (*ei)->hedge(1)->facet()->normal();
    bisector.normalize();
    
    update_beta( ADM::der_cast( (*ei)->hedge(0)->facet() ), bisector );
    update_beta( ADM::der_cast( (*ei)->hedge(1)->facet() ), bisector );
    }
    }
    
    // Solving quadratic equation
    void quadratic(double a, double b, double c, double& x1, double& x2)
    {
    if ( TML::zero(a,0.) )
    {
    if ( TML::zero(b,0.) )
    {
    // constant equation
    x1 = x2 = 1.;
    return;
    }
    
    // linear equation
    x1 = x2 = - c / b;
    return;
    }
    
    double delta = b*b - 4*a*c;
    if (delta<0.)
    {
    // complex roots
    x1 = x2 = 1.;
    return;
    }
    
    x1 = ( -b + std::sqrt(delta) ) / (2*a);
    x2 = ( -b - std::sqrt(delta) ) / (2*a);
    }
  */
};

#endif
