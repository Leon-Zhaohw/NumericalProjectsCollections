#include "lp3d.h"
#include "multi_criteria.h"
#include "drawing.h"

int main(int argc, char** argv)
{
  Lp3D surf( 2., 1. );

  MultiCriteria criteria( &surf, 0.2 );

  ADM::AdaptiveMesh mesh( "../base_mesh/cube.off", &criteria );
  mesh.set_threshold( 0.002 );
  mesh.set_stiffness( -10 );

  ADMViewer win( "MultiSphere", &mesh, &surf, &criteria );
  //win.KeyboardOptions();
  win.show();
  Fl::run();

  return 1;
}
