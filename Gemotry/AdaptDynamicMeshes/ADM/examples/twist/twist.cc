#include "twist.h"

Twist::Twist(double p, double r, double theta) 
  : Elp(p,r), Hlp(3,3), Vlp(3,1), Dlp(3,1), 
    Theta(theta), Add(true) { }

Twist::~Twist() { } 

/*------------------------------------------------------------*/

double Twist::F(const ADM::R3& p)
{ return Elp.F( T(p) ); }

double Twist::Fx(const ADM::R3& p)
{
  ADM::R3 g = Elp.gradient( T(p) );
  return dot( g, Tx(p) );
}

double Twist::Fy(const ADM::R3& p)
{
  ADM::R3 g = Elp.gradient( T(p) );
  return dot( g, Ty(p) );
}

double Twist::Fz(const ADM::R3& p)
{
  ADM::R3 g = Elp.gradient( T(p) );
  return dot( g, Tz(p) );
}

double Twist::Fxx(const ADM::R3& p)
{
  ADM::R3 t = T(p);

  Hlp[0][0] = Elp.Fxx(t);   Hlp[0][1] = Elp.Fxy(t);   Hlp[0][2] = Elp.Fxz(t);
  Hlp[1][0] = Elp.Fxy(t);   Hlp[1][1] = Elp.Fyy(t);   Hlp[1][2] = Elp.Fyz(t);
  Hlp[2][0] = Elp.Fxz(t);   Hlp[2][1] = Elp.Fyz(t);   Hlp[2][2] = Elp.Fzz(t);
  
  ADM::R3 tx = Tx(p);
  Vlp[0][0] = tx.x;  Vlp[1][0] = tx.y;  Vlp[2][0] = tx.z;
  
  Dlp = TNT::matmult( Hlp, Vlp );
  ADM::R3 a = ADM::R3( Dlp[0][0], Dlp[1][0], Dlp[2][0] );
  
  ADM::R3 g = Elp.gradient(t);
  return ( dot( tx, a ) + dot( g, Txx(p) ) );
}

double Twist::Fxy(const ADM::R3& p)
{
  ADM::R3 t = T(p);
  
  Hlp[0][0] = Elp.Fxx(t);   Hlp[0][1] = Elp.Fxy(t);   Hlp[0][2] = Elp.Fxz(t);
  Hlp[1][0] = Elp.Fxy(t);   Hlp[1][1] = Elp.Fyy(t);   Hlp[1][2] = Elp.Fyz(t);
  Hlp[2][0] = Elp.Fxz(t);   Hlp[2][1] = Elp.Fyz(t);   Hlp[2][2] = Elp.Fzz(t);
  
  ADM::R3 ty = Ty(p);
  Vlp[0][0] = ty.x;  Vlp[1][0] = ty.y;  Vlp[2][0] = ty.z;
  
  Dlp = TNT::matmult( Hlp, Vlp );
  ADM::R3 a = ADM::R3( Dlp[0][0], Dlp[1][0], Dlp[2][0] );
  
  ADM::R3 tx = Tx(p);
  ADM::R3 g = Elp.gradient(t);
  return ( dot( tx, a ) + dot( g, Txy(p) ) );
}

double Twist::Fxz(const ADM::R3& p)
{
  ADM::R3 t = T(p);

  Hlp[0][0] = Elp.Fxx(t);   Hlp[0][1] = Elp.Fxy(t);   Hlp[0][2] = Elp.Fxz(t);
  Hlp[1][0] = Elp.Fxy(t);   Hlp[1][1] = Elp.Fyy(t);   Hlp[1][2] = Elp.Fyz(t);
  Hlp[2][0] = Elp.Fxz(t);   Hlp[2][1] = Elp.Fyz(t);   Hlp[2][2] = Elp.Fzz(t);
 
  ADM::R3 tz = Tz(p);
  Vlp[0][0] = tz.x;  Vlp[1][0] = tz.y;  Vlp[2][0] = tz.z;
  
  Dlp = TNT::matmult( Hlp, Vlp );
  ADM::R3 a = ADM::R3( Dlp[0][0], Dlp[1][0], Dlp[2][0] );
  
  ADM::R3 tx = Tx(p);
  ADM::R3 g = Elp.gradient(t);
  return ( dot( tx, a ) + dot( g, Txz(p) ) );
}

double Twist::Fyy(const ADM::R3& p)
{
  ADM::R3 t = T(p);

  Hlp[0][0] = Elp.Fxx(t);   Hlp[0][1] = Elp.Fxy(t);   Hlp[0][2] = Elp.Fxz(t);
  Hlp[1][0] = Elp.Fxy(t);   Hlp[1][1] = Elp.Fyy(t);   Hlp[1][2] = Elp.Fyz(t);
  Hlp[2][0] = Elp.Fxz(t);   Hlp[2][1] = Elp.Fyz(t);   Hlp[2][2] = Elp.Fzz(t);
  
  ADM::R3 ty = Ty(p);
  Vlp[0][0] = ty.x;  Vlp[1][0] = ty.y;  Vlp[2][0] = ty.z;
  
  Dlp = TNT::matmult( Hlp, Vlp );
  ADM::R3 a = ADM::R3( Dlp[0][0], Dlp[1][0], Dlp[2][0] );
  
  ADM::R3 g = Elp.gradient(t);
  return ( dot( ty, a ) + dot( g, Tyy(p) ) );
}

double Twist::Fyz(const ADM::R3& p)
{
  ADM::R3 t = T(p);

  Hlp[0][0] = Elp.Fxx(t);   Hlp[0][1] = Elp.Fxy(t);   Hlp[0][2] = Elp.Fxz(t);
  Hlp[1][0] = Elp.Fxy(t);   Hlp[1][1] = Elp.Fyy(t);   Hlp[1][2] = Elp.Fyz(t);
  Hlp[2][0] = Elp.Fxz(t);   Hlp[2][1] = Elp.Fyz(t);   Hlp[2][2] = Elp.Fzz(t);

  ADM::R3 tz = Tz(p);
  Vlp[0][0] = tz.x;  Vlp[1][0] = tz.y;  Vlp[2][0] = tz.z;
  
  Dlp = TNT::matmult( Hlp, Vlp );
  ADM::R3 a = ADM::R3( Dlp[0][0], Dlp[1][0], Dlp[2][0] );
  
  ADM::R3 ty = Ty(p);
  ADM::R3 g = Elp.gradient(t);
  return ( dot( ty, a ) + dot( g, Tyz(p) ) );
}

double Twist::Fzz(const ADM::R3& p)
{
  ADM::R3 t = T(p);

  Hlp[0][0] = Elp.Fxx(t);   Hlp[0][1] = Elp.Fxy(t);   Hlp[0][2] = Elp.Fxz(t);
  Hlp[1][0] = Elp.Fxy(t);   Hlp[1][1] = Elp.Fyy(t);   Hlp[1][2] = Elp.Fyz(t);
  Hlp[2][0] = Elp.Fxz(t);   Hlp[2][1] = Elp.Fyz(t);   Hlp[2][2] = Elp.Fzz(t);

  ADM::R3 tz = Tz(p);
  Vlp[0][0] = tz.x;  Vlp[1][0] = tz.y;  Vlp[2][0] = tz.z;
  
  Dlp = TNT::matmult( Hlp, Vlp );
  ADM::R3 a = ADM::R3( Dlp[0][0], Dlp[1][0], Dlp[2][0] );
  
  ADM::R3 g = Elp.gradient(t);
  return ( dot( tz, a ) + dot( g, Tzz(p) ) );
}

/*------------------------------------------------------------*/

ADM::R3 Twist::skeleton(const ADM::R3& p)
{ return Elp.skeleton( T(p) ); }

/*------------------------------------------------------------*/

ADM::R3 Twist::T(const ADM::R3& p)
{ 
  double phi = getPhase(p.x);
  double co=std::cos(phi), si=std::sin(phi);
  return ADM::R3(p.x, -p.y*si + p.z*co, p.y*co + p.z*si);
}

ADM::R3 Twist::Tx(const ADM::R3& p)
{
  double phi = getPhase(p.x);
  double co=std::cos(phi), si=std::sin(phi);
  return ADM::R3(1, -p.y*Theta*co - p.z*Theta*si, -p.y*Theta*si + p.z*Theta*co);
}

ADM::R3 Twist::Ty(const ADM::R3& p)
{
  double phi = getPhase(p.x);
  return ADM::R3(0, -std::sin(phi), std::cos(phi));
}

ADM::R3 Twist::Tz(const ADM::R3& p)
{
  double phi = getPhase(p.x);
  return ADM::R3(0, std::cos(phi), std::sin(phi));  
}

ADM::R3 Twist::Txx(const ADM::R3& p)
{
  double phi = getPhase(p.x);
  double co=std::cos(phi), si=std::sin(phi);
  return ( std::pow(Theta,2) * ADM::R3(0, p.y*si - p.z*co, -p.y*co - p.z*si) );
}

ADM::R3 Twist::Txy(const ADM::R3& p)
{
  double phi = getPhase(p.x);
  return ADM::R3(0, -Theta*std::cos(phi), -Theta*std::sin(phi));  
}

ADM::R3 Twist::Txz(const ADM::R3& p)
{
  double phi = getPhase(p.x);
  return ADM::R3(0, -Theta*std::sin(phi), Theta*std::cos(phi));  
}

ADM::R3 Twist::Tyy(const ADM::R3& p)
{ return ADM::R3(); }

ADM::R3 Twist::Tyz(const ADM::R3& p)
{ return ADM::R3(); }


ADM::R3 Twist::Tzz(const ADM::R3& p)
{ return ADM::R3(); }
