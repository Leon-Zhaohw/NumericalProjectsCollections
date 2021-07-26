#ifndef _IMPLICIT_SURF_GENUS0_
#define _IMPLICIT_SURF_GENUS0_ 1

#include "implicit.h"

class ImplicitSurfGenus0 : public ImplicitSurf
{
 public:
  ImplicitSurfGenus0() { }
  
  virtual ~ImplicitSurfGenus0() { }
      
  /*------------------------------------------------------------*/

  ADM::R3 position(const ADM::R3& p)
    {
      double val = F(p);
      if ( TML::zero( val, 0., 1e-6 ) ) return p;
      
      ADM::R3 c = skeleton(p);
      ADM::R3 d = p - c;
      double a=0., b=1., t;
	  
      if (val<0.)
	{
	  while ( F( c+d*b ) < 0. ) b++;
	  a = b-1;
	}
      t = (a+b)*0.5;
      
      while ( !TML::zero( val=F( c+d*t ), 0., 1e-3 ) )
	{
	  if ( val*F( c+d*b ) < 0. ) a = t;
	  else                       b = t;
	  t = (a+b)/2;
	}
      
      return ( c+d*t );
    }
  
  /*------------------------------------------------------------*/
  
  virtual ADM::R3 skeleton(const ADM::R3& p) = 0;
};

#endif
