#ifndef _ELLIPSOID_SURF_
#define _ELLIPSOID_SURF_ 1

#include "surfaces/implicit_genus0.h"

class Ellipsoid : public ImplicitSurfGenus0
{
 protected:
  double A, B, C;
  int   State;
  
 public:
  Ellipsoid(double a, double b, double c)
    : A(a), B(b), C(c), State(0) { }
  
  virtual ~Ellipsoid() { }
  
  /*------------------------------------------------------------*/
  
  void deform()
    {
      int i = State % 3;
      bool add = (State<3)? true : false;
      
      switch (i)
	{
	case 0:
	  if (add)
	    {
	      A += 0.05;
	      if (A>=2) State += (State==2)? 3 : 1;
	    } else {
	      A -= 0.05;
	      if (A<=1) State -= (State==3)? 3 : 1;
	    }
	  break;
	case 1:
	  if (add)
	    {
	      B += 0.05;
	      if (B>=2) State += (State==2)? 3 : 1;
	    } else {
	      B -= 0.05;
	      if (B<=1) State -= (State==3)? 3 : 1;
	    }
	  break;
	case 2:
	  if (add)
	    {
	      C += 0.05;
	      if (C>=2) State += (State==2)? 3 : 1;
	    } else {
	      C -= 0.05;
	      if (C<=1) State -= (State==3)? 3 : 1;
	    }
	  break;
	default: break;
	};
    }

  /*------------------------------------------------------------*/
  
  double F(const ADM::R3& p)   { return ( std::pow(p.x/A,2) + std::pow(p.y/B,2) + std::pow(p.z/C,2) - 1 ); }
  
  double Fx(const ADM::R3& p)  { return ( 2*p.x / std::pow(A,2) ); }
  
  double Fy(const ADM::R3& p)  { return ( 2*p.y / std::pow(B,2) ); }
  
  double Fz(const ADM::R3& p)  { return ( 2*p.z / std::pow(C,2) ); }
  
  double Fxx(const ADM::R3& p) { return ( 2 / std::pow(A,2) ); }
  
  double Fxy(const ADM::R3& p) { return 0; }
  
  double Fxz(const ADM::R3& p) { return 0; }
  
  double Fyy(const ADM::R3& p) { return ( 2 / std::pow(B,2) ); }
  
  double Fyz(const ADM::R3& p) { return 0; }
  
  double Fzz(const ADM::R3& p) { return ( 2 / std::pow(C,2) ); }

  /*------------------------------------------------------------*/

  ADM::R3 skeleton(const ADM::R3& p) { return ADM::R3(); }
};

#endif
