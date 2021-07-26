#ifndef _TWIST_SURF_
#define _TWIST_SURF_ 1

#include "../extrudelp/extrudelp.h"

class Twist : public ImplicitSurfGenus0
{
 protected:
  ExtrudeLP   Elp;
  ADM::Matrix Hlp, Vlp, Dlp;
  double      Theta;
  bool        Add;
  
 public:
  Twist(double p, double r, double theta);
  
  virtual ~Twist();
  
  /*------------------------------------------------------------*/
  
  void deform()
    {
      if (Add) Theta += 0.1;
      else     Theta -= 0.1;
      
      if (Theta>M_PI || Theta<0.) Add = !Add;
      if (Theta<0.) Theta = 0.;
      cout << "Theta = " <<  Theta << endl;
    }

  /*------------------------------------------------------------*/
  
  double F(const ADM::R3& p);
  double Fx(const ADM::R3& p);
  double Fy(const ADM::R3& p);
  double Fz(const ADM::R3& p);
  double Fxx(const ADM::R3& p);
  double Fxy(const ADM::R3& p);
  double Fxz(const ADM::R3& p);
  double Fyy(const ADM::R3& p);
  double Fyz(const ADM::R3& p);
  double Fzz(const ADM::R3& p);

  /*------------------------------------------------------------*/

  ADM::R3 skeleton(const ADM::R3& p);
  
 /*------------------------------------------------------------*/

  inline double getPhase(double x) { return Theta*(x+1); }
  
  ADM::R3 T(const ADM::R3& p);
  ADM::R3 Tx(const ADM::R3& p);
  ADM::R3 Ty(const ADM::R3& p);
  ADM::R3 Tz(const ADM::R3& p);
  ADM::R3 Txx(const ADM::R3& p);
  ADM::R3 Txy(const ADM::R3& p);
  ADM::R3 Txz(const ADM::R3& p);
  ADM::R3 Tyy(const ADM::R3& p);
  ADM::R3 Tyz(const ADM::R3& p);
  ADM::R3 Tzz(const ADM::R3& p);
};

#endif
