#ifndef _LP_3D_
#define _LP_3D_ 1

#include "surfaces/implicit_genus0.h"

class Lp3D : public ImplicitSurfGenus0
{
 protected:
  double P, R;
  bool up;
 
 public:
  Lp3D(double p, double r) 
    : P(p), R(r), up(true) 
    { if (P<2.) { TML::Error("I only accept P>=2"); } }

  virtual ~Lp3D() { }
  
  /*------------------------------------------------------------*/
  
  void deform()
    {
      if (up) P += 0.5;
      else    P -= 0.5;

      if (P >= 18.) up = false;
      if (P <= 2. ) up = true;

      //P +=0.5;
      std::cout << "p = " << P << endl;
    }

  /*------------------------------------------------------------*/
  
  double F(const ADM::R3& p)
    { 
      return ( std::pow(std::fabs(p.x),P) + 
	       std::pow(std::fabs(p.y),P) + 
	       std::pow(std::fabs(p.z),P) - 
	       std::pow(R,P) ); 
    }
  
  double Fx(const ADM::R3& p)
    { double val = P * std::pow( std::fabs(p.x), P-1 ); return (p.x>=0)? val : -val; }
  
  double Fy(const ADM::R3& p)
    { double val = P * std::pow( std::fabs(p.y), P-1 ); return (p.y>=0)? val : -val; }
  
  double Fz(const ADM::R3& p)
    { double val = P * std::pow( std::fabs(p.z), P-1 ); return (p.z>=0)? val : -val; }
  
  double Fxx(const ADM::R3& p)
    { return ( P * (P-1) * std::pow( std::fabs(p.x), P-2 ) ); }
  
  double Fxy(const ADM::R3& p)
    { return 0; }
  
  double Fxz(const ADM::R3& p)
    { return 0; }
  
  double Fyy(const ADM::R3& p)
    { return ( P * (P-1) * std::pow( std::fabs(p.y), P-2 ) ); }
  
  double Fyz(const ADM::R3& p)
    { return 0; }
  
  double Fzz(const ADM::R3& p)
    { return ( P * (P-1) * std::pow( std::fabs(p.z), P-2 ) ); }
  
  /*------------------------------------------------------------*/

  ADM::R3 skeleton(const ADM::R3& p)
    { return ADM::R3(); }
};

#endif
