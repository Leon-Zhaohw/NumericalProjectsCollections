#include "ellipsoid.h"
#include "drawing.h"

int main(int argc, char** argv)
{
  Ellipsoid surf( 1.,1.,1. );
  
  ADM::StochasticSampling criteria( &surf );

  ADM::AdaptiveMesh mesh( "../base_mesh/cube.off", &criteria );
  mesh.set_threshold( 0.02 ); //0.05 );

  ADMViewer win( "Ellipsoid", &mesh, &surf, &criteria );
  //win.KeyboardOptions();
  win.show();
  Fl::run();
  
  return 1;
}
