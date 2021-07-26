#ifndef _BUMP_SINC_
#define _BUMP_SINC_ 1

#include "sinc.h"
#include "drawing.h"

int main(int argc, char** argv)
{
  BumpSinc surf( 0.5, 3*M_PI, 0.9, 0., 0. );
  
  ADM::StochasticSampling criteria( &surf );

  ADM::AdaptiveMesh mesh( "../base_mesh/plane.off", &criteria );
  mesh.set_stiffness( -0.1 );

  ADMViewer win( "Bump_sinc", &mesh, &surf, &criteria );
  //win.KeyboardOptions();
  win.show();
  Fl::run();
  
  return 1;
}

#endif
