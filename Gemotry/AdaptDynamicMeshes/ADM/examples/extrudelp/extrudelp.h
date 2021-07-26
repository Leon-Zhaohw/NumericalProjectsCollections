#ifndef _EXTRUDE_LP_SURF_
#define _EXTRUDE_LP_SURF_ 1

#include "surfaces/implicit_genus0.h"

class ExtrudeLP : public ImplicitSurfGenus0
{
 protected:
  double P, R;
 
 public:
  ExtrudeLP(double p, double r) 
    : P(p), R(r) 
    { if (P<2.) { TML::Error("I only accept P>=2"); } }

  virtual ~ExtrudeLP() { }
  
 /*------------------------------------------------------------*/
  
  void deform()
    {
      P +=1;
      std::cout << "p = " << P << std::endl;
    }

  /*------------------------------------------------------------*/
  
  double F(const ADM::R3& p)
    { return ( std::pow(std::fabs(p.y),P) + std::pow(std::fabs(p.z),P) - std::pow(R,P) ); }

  double Fx(const ADM::R3& p)
    { return 0; }

  double Fy(const ADM::R3& p)
    { double val = P * std::pow( std::fabs(p.y), P-1 ); return (p.y>=0)? val : -val; }

  double Fz(const ADM::R3& p)
    { double val = P * std::pow( std::fabs(p.z), P-1 ); return (p.z>=0)? val : -val; }

  double Fxx(const ADM::R3& p)
    { return 0; }

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
    { return ADM::R3(p.x,0,0); }
};

#endif
