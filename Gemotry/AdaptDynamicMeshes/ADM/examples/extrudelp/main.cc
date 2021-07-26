#include "extrudelp.h"
#include "drawing.h"

int main(int argc, char** argv)
{
  ExtrudeLP surf( 2., 1. );

  ADM::StochasticSampling criteria( &surf );

  ADM::AdaptiveMesh mesh( "../base_mesh/cylinder.off", &criteria );

  ADMViewer win( "ExtrudeLp", &mesh, &surf, &criteria );
  //win.KeyboardOptions();
  win.show();
  Fl::run();
  
  return 1;
}
