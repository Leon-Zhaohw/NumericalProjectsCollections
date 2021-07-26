#ifndef _IMPLICIT_SURF_
#define _IMPLICIT_SURF_ 1

#include <adm.h>

class ImplicitSurf : public ADM::Oracle
{
 protected:
  ADM::Matrix I, H, A, W, V;
  ADM::Array  E;
  
 public:
  ImplicitSurf() 
    : I(3,3), H(3,3), A(3,3), W(3,3), V(3,3), E(3) 
    { 
      for (int i=0; i<3; i++) 
	for (int j=0; j<3; j++) 
	  I[i][j] = (i==j)? 1. : 0.;
    }
  
  virtual ~ImplicitSurf() { }
        
  /*------------------------------------------------------------*/

  void update(ADM::A48Vertex* v)
    { project(v); }
  
  void project(ADM::A48Vertex* v, ADM::R3 n=ADM::R3())
    { v->set_pos( position(v->pos()) ); }
  
  double error(const ADM::R3& p)
    { return (position(p) - p).length(); }
      
  /*------------------------------------------------------------*/

  void normal(ADM::A48Vertex* v)
    { v->set_normal( normal(v->pos()) ); }
      
  void curvature(ADM::A48Vertex* v)
    {
      double  k1, k2;
      ADM::R3 d1, d2;
      curvature( v->pos(), k1, k2, d1, d2 );

      v->set_max_curv(k1);     v->set_min_curv(k2);
      v->set_max_direc(d1);    v->set_min_direc(d2);
      v->set_gaus_curv(k1*k2); v->set_mean_curv(0.5*(k1+k2));
    }
  
  /*------------------------------------------------------------*/

  ADM::R3 normal(const ADM::R3& p)
    {
      ADM::R3 g = gradient(p);
      g.normalize();
      return g;
    }

  void curvature(const ADM::R3& p, double& k1, double& k2, ADM::R3& d1, ADM::R3& d2)
    {
      /***********************************************************
       * Computing principal curvatures and directions following *
       * 'Level Set Surface Editing Operators'                   *
       * K. Museth, D. E. Breen, R. T. Whitaker, A. H. Barr      *
       ***********************************************************/
      ADM::R3 g = gradient(p);
      double len = g.length();
      ADM::R3 n = g/len;
	  
      A = I - extprod(n,n);
      H[0][0] = Fxx(p)/len; H[0][1] = Fxy(p)/len; H[0][2] = Fxz(p)/len;
      H[1][0] = Fxy(p)/len; H[1][1] = Fyy(p)/len; H[1][2] = Fyz(p)/len;
      H[2][0] = Fxz(p)/len; H[2][1] = Fyz(p)/len; H[2][2] = Fzz(p)/len;

      W = TNT::matmult( H, A );
	  
      ADM::Eigen eig(W);
      eig.getRealEigenvalues(E);
      eig.getV(V);
	  
      ADM::R3 v[3];
      for (int i=0; i<3; i++)
	{
	  v[i] = ADM::R3(V[0][i],V[1][i],V[2][i]);
	  v[i].normalize();
	}
      
      // Ordering eigenvalues
      double abs_dot[3];
      for (int i=0; i<3; i++)
	abs_dot[i] = std::fabs( dot( n, v[i] ) );
      
      int indices[3];
      if ( abs_dot[0]>abs_dot[1] && abs_dot[0]>abs_dot[2] )
	{
	  indices[0] = 0;
	  if (E[1]<E[2])
	    { indices[1] = 1; indices[2] = 2; }
	  else
	    { indices[1] = 2; indices[2] = 1; }
	}
      else if ( abs_dot[1]>abs_dot[0] && abs_dot[1]>abs_dot[2] )
      {
	indices[0] = 1;
	if (E[0]<E[2])
	  { indices[1] = 0; indices[2] = 2; }
	else
	  { indices[1] = 2; indices[2] = 0; }
      } else {
	indices[0] = 2;
	if (E[0]<E[1])
	  { indices[1] = 0; indices[2] = 1; }
	else
	  { indices[1] = 1; indices[2] = 0; }
      }
      
      // Setting principal curvatures and directions
      k1 = E[ indices[2] ];
      k2 = E[ indices[1] ];
      
      d1 = v[ indices[2] ];
      d2 = v[ indices[1] ];
    }
  
  ADM::R3 gradient(const ADM::R3& p) 
    { return ADM::R3(Fx(p),Fy(p),Fz(p)); }      
  
  /*------------------------------------------------------------*/

  virtual ADM::R3 position(const ADM::R3& p) = 0;

  virtual double F(const ADM::R3& p)   = 0;
  virtual double Fx(const ADM::R3& p)  = 0;
  virtual double Fy(const ADM::R3& p)  = 0;
  virtual double Fz(const ADM::R3& p)  = 0;
  virtual double Fxx(const ADM::R3& p) = 0;
  virtual double Fxy(const ADM::R3& p) = 0;
  virtual double Fxz(const ADM::R3& p) = 0;
  virtual double Fyy(const ADM::R3& p) = 0;
  virtual double Fyz(const ADM::R3& p) = 0;
  virtual double Fzz(const ADM::R3& p) = 0;
};

#endif
